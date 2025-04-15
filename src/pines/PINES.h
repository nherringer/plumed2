/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2021, Andrea Arsiccio

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_PINES_vec_PINES_h
#define __PLUMED_PINES_vec_PINES_h
#include <unordered_map>
#include <set>
#include <fstream>
#include "tools/AtomNumber.h"

using namespace std;

namespace PLMD {
namespace PINES {

struct AtomNumberLess {
  bool operator()(const AtomNumber& a, const AtomNumber& b) const {
    return a.index() < b.index();
  }
};    
// Ideally core/Colvar.h should be moved to this directory and Colvar should stay in namespace PLMD::Sasa
// With this trick, PLMD::Colvar is visible as PLMD::Sasa::Colvar
using PLMD::Colvar;

class PINES      : public Colvar
{
private:

  int N_Blocks;
  int total_PIV_length;
  std::vector<int> steps_since_update;
  std::vector<int> nstride;
  std::string ref_file;
  std::vector<SwitchingFunction> sfs;
  std::vector<string> sw;
  std::vector<double> r00;
  std::vector<std:: vector<double> > PIV;
  std::vector<std:: vector<Vector> > ann_deriv;

  std::vector<std:: vector<AtomNumber> > listall;
  std::vector<std:: vector<AtomNumber> > listreduced;
  std::set<AtomNumber, AtomNumberLess> listreducedall;
  std::vector<AtomNumber> listreducedall_vec;
  std::unordered_map<int,int> atom_ind_hashmap;

  std::vector<bool> stale_tolerance;
  PDB mypdb;
  std::vector<string> block_params;
  std::vector<std::vector<std::vector<AtomNumber> > > block_groups_atom_list;
  std::vector<int> block_lengths;
  std::vector<int> Buffer_Pairs;
  std::vector<int> tot_num_pairs;
  std::vector<std::vector<std::vector<bool> > > input_filters;
  std::vector<double> delta_pd;
  std::vector<double> r_tolerance;
  std::vector<std::vector<Vector> > PL_atoms_ref_coords;
  std::vector<std::vector<std::pair<AtomNumber,AtomNumber> > > Exclude_Pairs;
  std::vector<std::vector<std::vector<string> > > Name_list;
  std::vector<std::vector<std::vector<AtomNumber> > > ID_list;
  std::vector<std::vector<std::vector<int> > > ResID_list;
  std::vector<std::vector<std::pair<double, std::pair<AtomNumber,AtomNumber> > > > vecMaxHeapVecs;

  bool atomMatchesFilters(int n, int g, AtomNumber ind, int resid, const std::string& atom_name);
  void buildMaxHeapVecBlock(int n, const PDB& mypdb, std::vector<std::pair<double, std::pair<AtomNumber, AtomNumber>>>& heap);
  void updateBlockPairList(int n, std::vector<std::pair<double, std::pair<AtomNumber, AtomNumber>>>& heap);
  double calculateDistance(const AtomNumber& ind0, const AtomNumber& ind1, const PDB& mypdb);
  std::ofstream log;
  void logMsg(const std::string& msg, const std::string& section);
  void logMsg(const Vector& vec, const std::string& section);
  void resizeAllContainers(int N);

public:
  static void registerKeywords( Keywords& keys );                                                                       
  explicit PINES(const ActionOptions&); 
  //~PINES();                                                                                                               
  // active methods:
  struct MinCompareDist {
    bool operator()(const std::pair<double, std::pair<AtomNumber, AtomNumber>>& p1, const std::pair<double, std::pair<AtomNumber, AtomNumber>>& p2) {
      return p1.first < p2.first; // Min heap
    }
  };
  struct MaxCompareDist {
    bool operator()(const std::pair<double, std::pair<AtomNumber, AtomNumber>>& p1, const std::pair<double, std::pair<AtomNumber, AtomNumber>>& p2) {
      return p1.first > p2.first; // Max heap
    }
  };
                                                                                            
  virtual void calculate();
  void checkFieldsAllowed() {}                                                                                           
  // -- SD prepare to requestAtoms during simulation 
  void prepare() override;
};

}
}

#endif
