/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2023 of Nicholas Herringer and Siva Dasetty.

The PINES module is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The PINES module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "colvar/Colvar.h"
#include "colvar/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionWithVirtualAtom.h"
#include "tools/NeighborList.h"
#include "tools/SwitchingFunction.h"
#include "tools/PDB.h"
#include "tools/Pbc.h"
#include "tools/Stopwatch.h"

// -- SD header file for PINES
#include "PINES.h"

#include <string>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <unordered_map>
#include <set>
#include <numeric>

using namespace std;

namespace PLMD
{
  namespace PINES
  {

    PLUMED_REGISTER_ACTION(PINES, "PINES")

    void PINES::registerKeywords(Keywords &keys)
    {
      Colvar::registerKeywords(keys);
      keys.add("numbered", "SWITCH", "The switching functions parameter. You must specify a Switching function for all PINES blocks."
                                     "Details of the various switching functions you can use are provided on \\ref switchingfunction.");
      keys.add("numbered", "BLOCK", "Each block of the PIV");
      keys.add("compulsory", "REF_FILE", "PDB file name that contains the information about system connectivity and labels.");
      keys.add("compulsory", "N_BLOCKS", "Number of blocks in PIV");
      keys.add("compulsory", "BUFFER", "Number of additional pairwise distances to include as a buffer for each PIV block");
      keys.add("compulsory", "PL_REFRESH", "Upper limit refresh rate for each block pair list");
      keys.add("compulsory", "SIZE", "Length of each PIV Block");
      keys.add("optional", "ATOMID", "AtomIDs");
      keys.add("optional", "RESID", "ResIDs");
      keys.add("optional", "NAME", "Atom Names");
      keys.add("optional", "EXCLUDE_PAIRS", "Excluded pairs");
      componentsAreNotOptional(keys);
      keys.addOutputComponent("ELEMENT", "default", "Elements of the PINES block"); 
      keys.reset_style("SWITCH", "compulsory");
    }

    bool PINES::atomMatchesFilters(int n, int g, AtomNumber ind, int resid, const std::string& atom_name) {
      bool id_check = !input_filters[n][g][0];
      bool res_check = !input_filters[n][g][1];
      bool name_check = !input_filters[n][g][2];
    
      if (input_filters[n][g][0]) {
        for (auto& atomID : ID_list[n][g]) {
          if (ind == atomID) { id_check = true; break; }
        }
      }
      if (input_filters[n][g][1]) {
        for (auto& rID : ResID_list[n][g]) {
          if (resid == rID) { res_check = true; break; }
        }
      }
      if (input_filters[n][g][2]) {
        for (auto& name : Name_list[n][g]) {
          if (atom_name == name) { name_check = true; break; }
        }
      }
      return id_check && res_check && name_check;
    }
    
    void PINES::buildMaxHeapVecBlock(int n, const PDB& mypdb, std::vector<std::pair<double, std::pair<AtomNumber, AtomNumber>>>& heap) {
      logMsg("Entering buildMaxHeapVecBlock", "buildMaxHeapVecBlock");
      logMsg("Block index: " + std::to_string(n), "buildMaxHeapVecBlock");
      logMsg("heap.size(): " + std::to_string(heap.size()), "buildMaxHeapVecBlock");
      logMsg("block_groups_atom_list[n][0].size(): " + std::to_string(block_groups_atom_list[n][0].size()), "buildMaxHeapVecBlock");
      logMsg("block_groups_atom_list[n][1].size(): " + std::to_string(block_groups_atom_list[n][1].size()), "buildMaxHeapVecBlock");

      bool isFirstStep = (getStep() == 0);

      std::vector<std::pair<AtomNumber, AtomNumber> > unique_pairs;
      heap.clear();
      logMsg("alpha0", "buildMaxHeapVecBlock");
      for (int i = 0; i < block_groups_atom_list[n][0].size(); i++) {
        logMsg("alpha1", "buildMaxHeapVecBlock");
        AtomNumber ind0 = block_groups_atom_list[n][0][i];
        logMsg("alpha2", "buildMaxHeapVecBlock");
        Vector Pos0 = isFirstStep ? mypdb.getPosition(ind0) : getPosition(atom_ind_hashmap[ind0.index()]);
        logMsg("alpha3", "buildMaxHeapVecBlock");
        for (int j = 0; j < block_groups_atom_list[n][1].size(); j++) {
          AtomNumber ind1 = block_groups_atom_list[n][1][j];

          if (ind1 == ind0) continue;

          auto test_pair = std::make_pair(ind0, ind1);
          auto reverse_pair = std::make_pair(ind1, ind0);
    
          if (std::find(Exclude_Pairs[n].begin(), Exclude_Pairs[n].end(), test_pair) != Exclude_Pairs[n].end() ||
              std::find(Exclude_Pairs[n].begin(), Exclude_Pairs[n].end(), reverse_pair) != Exclude_Pairs[n].end() ||
              std::find(unique_pairs.begin(), unique_pairs.end(), reverse_pair) != unique_pairs.end()) {
            continue;
          } else {
            unique_pairs.push_back(test_pair);
          }
    
          Vector Pos1 = isFirstStep ? mypdb.getPosition(ind1) : getPosition(atom_ind_hashmap[ind1.index()]);
          double mag = pbcDistance(Pos0, Pos1).modulo();

          heap.push_back({mag, {ind0, ind1}});
          std::push_heap(heap.begin(), heap.end(), MaxCompareDist());
          logMsg("alpha4", "buildMaxHeapVecBlock");
          if (heap.size() > tot_num_pairs[n]) {
            std::pop_heap(heap.begin(), heap.end(), MaxCompareDist());
            heap.pop_back();
          }
        }
      }
      logMsg("alpha5", "buildMaxHeapVecBlock");
      std::sort(heap.begin(), heap.end(), MinCompareDist());
      logMsg("alpha6", "buildMaxHeapVecBlock");
      int chk1 = block_lengths[n]-1;
      int chk2 = chk1 + Buffer_Pairs[n];
      delta_pd[n] = heap[chk1].first - heap[chk2].first;
      r_tolerance[n] = delta_pd[n]/4;
      listreduced[n].clear();
      std::set<AtomNumber, AtomNumberLess> uniqueAtoms;
      logMsg("alpha7", "buildMaxHeapVecBlock");
      for (const auto& pair : heap) {
        uniqueAtoms.insert(pair.second.first);  // Atom1
        uniqueAtoms.insert(pair.second.second); // Atom2
      }
      for (const auto& atom : uniqueAtoms) {
        listreduced[n].push_back(atom);
    }
      logMsg("alpha8", "buildMaxHeapVecBlock");
      logMsg("listreduced[" + std::to_string(n) + "].size(): " + std::to_string(listreduced[n].size()), "buildMaxHeapVecBlock");
      PL_atoms_ref_coords[n].resize(listreduced[n].size());
      logMsg("PL_atoms_ref_coords[" + std::to_string(n) + "].size(): " + std::to_string(PL_atoms_ref_coords[n].size()), "buildMaxHeapVecBlock");
      for (int i=0; i<listreduced[n].size(); i++)
      {
        PL_atoms_ref_coords[n][i].zero();
        logMsg("listreduced[" + std::to_string(n) + "][" + std::to_string(i) + "].index(): " + std::to_string(listreduced[n][i].index()), "buildMaxHeapVecBlock");
        Vector debug_pos = isFirstStep ? mypdb.getPosition(listreduced[n][i]) : getPosition(atom_ind_hashmap[listreduced[n][i].index()]);
        logMsg(debug_pos, "buildMaxHeapVecBlock");
        PL_atoms_ref_coords[n][i] = isFirstStep ? mypdb.getPosition(listreduced[n][i]) : getPosition(atom_ind_hashmap[listreduced[n][i].index()]);
      }
      logMsg("alpha9", "buildMaxHeapVecBlock");
      ann_deriv.resize(listreduced[n].size());
      for (int i = 0; i < ann_deriv.size(); i++)
      {
        ann_deriv[i].resize(total_PIV_length);
      }
      logMsg("alpha10", "buildMaxHeapVecBlock");
    }
    
    void PINES::updateBlockPairList(int n, std::vector<std::pair<double, std::pair<AtomNumber, AtomNumber>>>& heap) {
      // This is separate from staleness, which triggers a refresh instead of an update.
      // This is simply to update the pair distances in the MaxHeap/Pairlist with the new atom positions
      // and rearrange the ordering in the MaxHeap/Pairlisat

      for (int i = 0; i < tot_num_pairs[n]; i++)
      {
        AtomNumber ind0 = heap[i].second.first;
        AtomNumber ind1 = heap[i].second.second;
        heap[i].first = calculateDistance(ind0,ind1,mypdb);
      }
      std::nth_element(heap.begin(), heap.begin() + block_lengths[n], heap.end(), MinCompareDist());
      std::sort(heap.begin(), heap.begin() + block_lengths[n], MinCompareDist());
    }
    
    double PINES::calculateDistance(const AtomNumber& ind0, const AtomNumber& ind1, const PDB& mypdb) {
      bool isFirstStep = (getStep() == 0);
      Vector Pos0 = isFirstStep ? mypdb.getPosition(ind0) : getPosition(atom_ind_hashmap[ind0.index()]);
      Vector Pos1 = isFirstStep ? mypdb.getPosition(ind1) : getPosition(atom_ind_hashmap[ind1.index()]);
      return pbcDistance(Pos0, Pos1).modulo();
    }

    void PINES::resizeAllContainers(int N) {
      // Outer containers
      nstride.resize(N);
      steps_since_update.resize(N);
      block_params.resize(N);
      block_groups_atom_list.resize(N);
      block_lengths.resize(N);
      Buffer_Pairs.resize(N);
      tot_num_pairs.resize(N);
      Exclude_Pairs.resize(N);
      vecMaxHeapVecs.resize(N);
      PIV.resize(N);
      listall.resize(N);
      listreduced.resize(N);
      stale_tolerance.resize(N);
      r00.resize(N);
      sw.resize(N);
      sfs.resize(N);
      delta_pd.resize(N, 0.0);
      r_tolerance.resize(N, 0.0);
      PL_atoms_ref_coords.resize(N);
      input_filters.resize(N);
      ID_list.resize(N);
      ResID_list.resize(N);
      Name_list.resize(N);
      atom_ind_hashmap.clear();
    
      // Per-block inner structures
      for (int n = 0; n < N; n++) {
        block_groups_atom_list[n].resize(2);
        ID_list[n].resize(2);
        ResID_list[n].resize(2);
        Name_list[n].resize(2);
        input_filters[n].resize(2);
        PL_atoms_ref_coords[n].resize(0);  // will be filled per atom if needed
      }
    
      // Final 3D input_filters init
      for (int n = 0; n < N; n++) {
        for (int g = 0; g < 2; g++) {
          input_filters[n][g].resize(3, false);  // [ID, ResID, Name]
        }
      }
    }    

    void PINES::logMsg(const std::string& msg, const std::string& section) {
      log << "[" << plumed.getStep() << "] "
          << "[" << section << "] " << msg << std::endl;
    }

    void PINES::logMsg(const Vector& vec, const std::string& section) {
      log << "[" << section << "] Vector = (" 
          << vec[0] << ", " << vec[1] << ", " << vec[2] << ")" << std::endl;
    }

    PINES::PINES(const ActionOptions &ao) : PLUMED_COLVAR_INIT(ao),
                                            N_Blocks(1),
                                            total_PIV_length(1),
                                            steps_since_update(std::vector<int>(1)),
                                            nstride(std::vector<int>(1,10)),
                                            ref_file(std::string()),
                                            sfs(),
                                            sw(std::vector<string>(1)),
                                            r00(std::vector<double>(1)),
                                            PIV(std::vector<std::vector<double>>(1)),
                                            ann_deriv(std::vector<std::vector<Vector> >(1)),
                                            listall(std::vector<std::vector<AtomNumber> >()),
                                            listreduced(std::vector<std::vector<AtomNumber> >()),
                                            listreducedall(std::set<AtomNumber, AtomNumberLess>()),
                                            listreducedall_vec(std::vector<AtomNumber>()),
                                            atom_ind_hashmap(),
                                            stale_tolerance(std::vector<bool>(1,false)),
                                            mypdb(),
                                            block_params(std::vector<string>(1)),
                                            block_groups_atom_list(std::vector<std::vector<std::vector<AtomNumber> > >()),
                                            block_lengths(std::vector<int>(1,1)),
                                            Buffer_Pairs(std::vector<int>(1,0)),
                                            tot_num_pairs(std::vector<int>(1,1)),
                                            input_filters(
                                              std::vector<std::vector<std::vector<bool> > >(
                                                1,  // outermost dimension
                                                std::vector<std::vector<bool>>(
                                                  1,  // middle dimension
                                                  std::vector<bool>(1, false)  // innermost dimension with value `false`
                                                )
                                              )
                                            ),
                                            delta_pd(std::vector<double>()),
                                            r_tolerance(std::vector<double>()),
                                            PL_atoms_ref_coords(),
                                            Exclude_Pairs(std::vector<std::vector<std::pair<AtomNumber, AtomNumber> > >()),
                                            Name_list(std::vector<std::vector<std::vector<string> > >()),
                                            ID_list(std::vector<std::vector<std::vector<AtomNumber> > >()),
                                            ResID_list(std::vector<std::vector<std::vector<int> > >()),
                                            vecMaxHeapVecs()
    {
      log.open("pines_debug.log", std::ios::out);
      if (!log.is_open()) {
        error("Failed to open debug log file.");
      }

      logMsg("Constructor", "Finished initializer list");
      // Reference PDB file from which atom names, types, ids, and initial positions are determined
      parse("REF_FILE", ref_file);
      FILE *fp = fopen(ref_file.c_str(), "r");
      if (fp != NULL)
      {
        mypdb.readFromFilepointer(fp, plumed.getAtoms().usingNaturalUnits(), 0.1 / atoms.getUnits().getLength());
        fclose(fp);
      }
      else error("Error in reference PDB file");

      // Create variable to get number of blocks
      parse("N_BLOCKS", N_Blocks);
      resizeAllContainers(N_Blocks);

      // Check that the correct number of Blocks are specified
      for (unsigned n = 0; n < N_Blocks; n++)
      {
        if (!parseNumbered("BLOCK", n + 1, block_params[n]))
          break;
      }

      // Parse blocks for keywords
      for (int n = 0; n < N_Blocks; n++)
      {
        string block_length;
        std::vector<string> ex_pairs_n;
        std::vector<string> g1_ids;
        std::vector<string> g2_ids;
        std::vector<string> g1_resids;
        std::vector<string> g2_resids;
        std::vector<string> g1_names;
        std::vector<string> g2_names;

        std::vector<string> block_data = Tools::getWords(block_params[n]);


        std::vector<string> G1_data;
        Tools::parseVector(block_data, "G1", G1_data);

        std::vector<string> G2_data;
        Tools::parseVector(block_data, "G2", G2_data);

        Tools::parseVector(G1_data, "ATOMID", g1_ids);

        for (int i = 0; i < g1_ids.size(); i++)
        {
          AtomNumber g1i;
          g1i.setIndex(std::stoi(g1_ids[i]));
          ID_list[n][0].push_back(g1i);
        }

        Tools::parseVector(G2_data, "ATOMID", g2_ids);

        for (int i = 0; i < g2_ids.size(); i++)
        {
          AtomNumber g2i;
          g2i.setIndex(std::stoi(g2_ids[i]));
          ID_list[n][1].push_back(g2i);
        }

        Tools::parseVector(G1_data, "RESID", g1_resids);
        for (int i = 0; i < g1_resids.size(); i++)
        {
          ResID_list[n][0].push_back(std::stoi(g1_resids[i]));
        }
        Tools::parseVector(G2_data, "RESID", g2_resids);
        for (int i = 0; i < g2_resids.size(); i++)
        {
          ResID_list[n][1].push_back(std::stoi(g2_resids[i]));
        }
        Tools::parseVector(G1_data, "NAME", g1_names);
        for (int i = 0; i < g1_names.size(); i++)
        {
          Name_list[n][0].push_back(g1_names[i]);
        }
        Tools::parseVector(G2_data, "NAME", g2_names);
        for (int i = 0; i < g2_names.size(); i++)
        {
          Name_list[n][1].push_back(g2_names[i]);
        }

        Tools::parse(block_data, "SIZE", block_length);
        block_lengths[n] = std::stoi(block_length);

        Tools::parseVector(block_data, "EXCLUDE_PAIRS", ex_pairs_n);
        if (!ex_pairs_n.empty())
        {
          for (int i = 0; i < ex_pairs_n.size() - 1; i+=2)
          {
            AtomNumber atom1; 
            atom1.setIndex(std::stoi(ex_pairs_n[i]));
            AtomNumber atom2; 
            atom2.setIndex(std::stoi(ex_pairs_n[i+1]));

            std::pair<AtomNumber, AtomNumber> excluded_pair;
            excluded_pair = {atom1, atom2};
            Exclude_Pairs[n].push_back(excluded_pair);
          }
        }
        string buffer_pairs;
        string pl_refresh;
        Tools::parse(block_data, "BUFFER", buffer_pairs);
        Tools::parse(block_data, "PL_REFRESH", pl_refresh);

        if (!buffer_pairs.empty())
        {
          Buffer_Pairs[n] = std::stoi(buffer_pairs);
        }
        else
        {
          Buffer_Pairs[n] = 0;
        }

        if (!pl_refresh.empty())
        {
          nstride[n] = std::stoi(pl_refresh);
        }
        else
        {
          nstride[n] = 25;
        }
        tot_num_pairs[n] = block_lengths[n] + Buffer_Pairs[n];
      }

      total_PIV_length = std::accumulate(block_lengths.begin(), block_lengths.end(), 0);

      for (int n = 0; n < N_Blocks; n++)
      {
        // pseudo-code ish
        for (int g = 0; g < 2; g++)
        {
          if (ID_list[n][g].size() > 0)
          {
            input_filters[n][g][0] = true;
          }
          if (ResID_list[n][g].size() > 0)
          {
            input_filters[n][g][1] = true;
          }
          if (Name_list[n][g].size() > 0)
          {
            input_filters[n][g][2] = true;
          }
        }
      }

      for (int n=0; n < N_Blocks; n++) listall[n].clear();

      for (int i=0; i < mypdb.getAtomNumbers().size(); i++)
      {
        AtomNumber ind = mypdb.getAtomNumbers()[i];
        int resid = mypdb.getResidueNumber(ind);
        string atom_name = mypdb.getAtomName(ind);
        for (int n = 0; n < N_Blocks; n++)
        {
          bool atom_added = false;
          for (int g = 0; g < 2; g++)
          {
            if (atomMatchesFilters(n, g, ind, resid, atom_name)) {
              block_groups_atom_list[n][g].push_back(ind);
              atom_added = true;
            }
          }
          if (atom_added) listall[n].push_back(ind);
        }
      }

      r00.resize(N_Blocks);
      sw.resize(N_Blocks);
      sfs.resize(N_Blocks);

      for (unsigned n = 0; n < N_Blocks; n++) if (!parseNumbered("SWITCH", n + 1, sw[n])) break;
      
      std::string errors;
      for (unsigned n = 0; n < N_Blocks; n++)
      {
        sfs[n].set(sw[n], errors);
        std::string num;
        Tools::convert(n + 1, num);
        if (errors.length() != 0) error("problem reading SWITCH" + num + " keyword : " + errors);
        r00[n] = sfs[n].get_r0();
      }
      checkRead();
      logMsg("Constructor", "Parsed all keywords/values");
      log.flush();

      int total_count = 0;
      for (int n = 0; n < N_Blocks; n++)
      {
        for (int i = 0; i < block_lengths[n]; i++)
        {
          string comp = "ELEMENT-" + to_string(total_count);
          addComponentWithDerivatives(comp);
          componentIsNotPeriodic(comp);
          total_count += 1;
        }
      }
    }

    void PINES::prepare(){
      logMsg("Prepare", "Made it to prepare");
      bool heap_refreshed = false;
      bool stale_refresh_prep = false;
      for (int n = 0; n < N_Blocks; n++)
      {
        if (steps_since_update[n] == 0){
          buildMaxHeapVecBlock(n, mypdb, vecMaxHeapVecs[n]);
          heap_refreshed = true;
        }
        else if(steps_since_update[n] >= nstride[n] || stale_tolerance[n]){
          listreduced[n] = listall[n];
          stale_refresh_prep = true;
          steps_since_update[n] = -1;
          stale_tolerance[n] = false;
        }
        steps_since_update[n]+=1;
      }
      if (heap_refreshed || stale_refresh_prep){
        // Collate all atoms from all lists
        listreducedall.clear();
        for (const auto& blockList : listreduced) {        // Loop over each block's vector
          listreducedall.insert(blockList.begin(), blockList.end()); // Insert all atoms into set
        }
        atom_ind_hashmap.clear();
        listreducedall_vec = std::vector<AtomNumber>(listreducedall.begin(), listreducedall.end());
        for (int i=0; i < listreducedall_vec.size(); i++) atom_ind_hashmap[listreducedall_vec[i].index()] = i;
        requestAtoms(listreducedall_vec);
        ann_deriv.resize(listreducedall_vec.size());

        for (int i=0; i < ann_deriv.size(); i++) ann_deriv[i].resize(total_PIV_length);
      }
    }
      
    void PINES::calculate()
    {
      logMsg("Calculate", "Made it to calculate");
#pragma region VarsAndToleranceCheck

      Vector ref_xyz, step_xyz;
      AtomNumber aID;
      float delta_r;

      for(int n = 0; n < N_Blocks; n++){
        if(steps_since_update[n] > 0 && steps_since_update[n] < nstride[n]){
          // No point in checking if ssu is 0 because ssu is reference
          // No point in checking if ssu is nstride - 1 because next step will always trigger update anyway
          // Things I'm pretending exist:
          //   Hashmap -> PL_atoms_ref_coords = {atomID: [x, y, z]}
          //   Hashmap -> PL_pairs_dists = {(atomID1, atomID2): dist}
          //   Float -> r_tolerance = PBCdist(vecBL,vecL)/4
          //   My confidence in a bright future
          for(int i = 0; i < listreduced[n].size(); i++){
            aID = listreduced[n][i];
            ref_xyz = PL_atoms_ref_coords[n][i];
            step_xyz = getPosition(atom_ind_hashmap[aID.index()]);
            delta_r = pbcDistance(ref_xyz, step_xyz).modulo();
            if(delta_r >= r_tolerance[n]){
              stale_tolerance[n] = true;
            }
          }
        }
      }
#pragma endregion
#pragma region UpdatePLs

      // Build ann_deriv
      for (unsigned j = 0; j < ann_deriv.size(); j++)
      {
        for (unsigned i = 0; i < ann_deriv[j].size(); i++)
        {
          for (unsigned k = 0; k < 3; k++)
          {
            ann_deriv[j][i][k] = 0.;
          }
        }
      }

      for (unsigned n = 0; n < N_Blocks; n++)
      {
        PIV[n].resize(block_lengths[n]);
        for (unsigned i = 0; i < block_lengths[n]; i++)
        {
          PIV[n][i] = 0.;
        }
      }

      int PINES_element = 0;
      for (unsigned n = 0; n < N_Blocks; n++) {
        if (steps_since_update[n] != 1){
          updateBlockPairList(n, vecMaxHeapVecs[n]);
        }
        for (int i=0; i < block_lengths[n]; i++) {
          double ds_element = 0.;
          double dfunc = 0.;
          int local_aid0, local_aid1;
          PIV[n][i] = sfs[n].calculate(vecMaxHeapVecs[n][i].first, dfunc);
          local_aid0 = atom_ind_hashmap[vecMaxHeapVecs[n][i].second.first.index()];
          local_aid1 = atom_ind_hashmap[vecMaxHeapVecs[n][i].second.second.index()];
          Vector Pos0 = getPosition(local_aid0);
          Vector Pos1 = getPosition(local_aid1);
          Vector dr_dcoord = pbcDistance(Pos0, Pos1) / vecMaxHeapVecs[n][i].first;
          ds_element = dfunc * vecMaxHeapVecs[n][i].first;
          ann_deriv[local_aid0][PINES_element] = -ds_element * dr_dcoord;
          ann_deriv[local_aid1][PINES_element] = ds_element * dr_dcoord;
          PINES_element += 1;
        }
      }
#pragma endregion
#pragma region parallelizationAndPassing
      if (comm.initialized())
      {
        int count = 0;
        for (unsigned j = 0; j < N_Blocks; j++)
        {
          for (unsigned i = 0; i < PIV[j].size(); i++)
          {
            count += 1;
          }
        }

        comm.Barrier();

        for (unsigned j = 0; j < N_Blocks; j++)
        {
          for (unsigned k = 0; k < PIV[j].size(); k++)
          {
            comm.Sum(PIV[j][k]);
            PIV[j][k] /= comm.Get_size();
          }
        }

        if (!ann_deriv.empty())
        {
          for (unsigned i = 0; i < ann_deriv.size(); i++)
          {
            for (unsigned j = 0; j < ann_deriv[j].size(); j++)
            {
              for (unsigned k = 0; k < 3; k++)
              {
                comm.Sum(ann_deriv[i][j][k]);
                ann_deriv[i][j][k] /= comm.Get_size();
              }
            }
          }
        }
      }

      // Pass values and derivates to next stage
      unsigned total_count = 0;
      for (unsigned j = 0; j < N_Blocks; j++)
      {
        for (unsigned i = 0; i < block_lengths[j]; i++)
        {
          string comp = "ELEMENT-" + to_string(total_count);
          Value *valueNew = getPntrToComponent(comp);
          valueNew->set(PIV[j][i]);
          for (unsigned k = 0; k < ann_deriv.size(); k++)
          {
            setAtomsDerivatives(valueNew, k, ann_deriv[k][total_count]);
          }
          total_count += 1;
        }
      }
#pragma endregion
    }
  }
}
