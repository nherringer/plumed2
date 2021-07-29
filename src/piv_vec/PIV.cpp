/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2017 of Pipolo Silvio and Fabio Pietrucci.

The piv module is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The piv module is distributed in the hope that it will be useful,
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

// -- SD header file for both ANN function and PIV
#include "PIV.h"

#include <string>
#include <cmath>
#include <iostream>
#include <stdio.h>

using namespace std;

namespace PLMD
{
namespace piv
{

//+PLUMEDOC PIVMOD_COLVAR PIV
/*
Calculates the PIV-distance.

PIV distance is the squared Cartesian distance between the PIV \cite gallet2013structural \cite pipolo2017navigating
associated to the configuration of the system during the dynamics and a reference configuration provided
as input (PDB file format).
PIV can be used together with \ref FUNCPATHMSD to define a path in the PIV space.

\par Examples

The following example calculates PIV-distances from three reference configurations in Ref1.pdb, Ref2.pdb and Ref3.pdb
and prints the results in a file named colvar.
Three atoms (PIVATOMS=3) with names (pdb file) A B and C are used to construct the PIV and all PIV blocks (AA, BB, CC, AB, AC, BC) are considered.
SFACTOR is a scaling factor that multiplies the contribution to the PIV-distance given by the single PIV block.
NLIST sets the use of neighbor lists for calculating atom-atom distances.
The SWITCH keyword specifies the parameters of the switching function that transforms atom-atom distances.
SORT=1 means that the PIV block elements are sorted (SORT=0 no sorting.)
Values for SORT, SFACTOR and the neighbor list parameters have to be specified for each block.
The order is the following: AA,BB,CC,AB,AC,BC. If ONLYDIRECT (ONLYCROSS) is used the order is AA,BB,CC (AB,AC,BC).
The sorting operation within each PIV block is performed using the counting sort algorithm, PRECISION specifies the size of the counting array.

\plumedfile
PIV ...
LABEL=Pivd1
PRECISION=1000
NLIST
REF_FILE=Ref1.pdb
PIVATOMS=3
ATOMTYPES=A,B,C
SFACTOR=0.3,0.5,1.0,0.2,0.2,0.2
SORT=1,1,1,1,1,1
SWITCH1={RATIONAL R_0=0.6 MM=12 NN=4}
SWITCH2={RATIONAL R_0=0.4 MM=10 NN=5}
SWITCH3={RATIONAL R_0=0.4 MM=10 NN=5}
SWITCH4={RATIONAL R_0=0.5 MM=12 NN=6}
SWITCH5={RATIONAL R_0=0.5 MM=12 NN=6}
SWITCH6={RATIONAL R_0=0.5 MM=12 NN=6}
NL_CUTOFF=0.8,0.6,0.6,0.7,0.7,0.7
NL_STRIDE=10,10,10,10,10,10
NL_SKIN=0.1,0.1,0.1,0.1,0.1,0.1
... PIV
PIV ...
LABEL=Pivd2
PRECISION=1000
NLIST
REF_FILE=Ref2.pdb
PIVATOMS=3
ATOMTYPES=A,B,C
SFACTOR=0.3,0.5,1.0,0.2,0.2,0.2
SORT=1,1,1,1,1,1
SWITCH1={RATIONAL R_0=0.6 MM=12 NN=4}
SWITCH2={RATIONAL R_0=0.4 MM=10 NN=5}
SWITCH3={RATIONAL R_0=0.4 MM=10 NN=5}
SWITCH4={RATIONAL R_0=0.5 MM=12 NN=6}
SWITCH5={RATIONAL R_0=0.5 MM=12 NN=6}
SWITCH6={RATIONAL R_0=0.5 MM=12 NN=6}
NL_CUTOFF=0.8,0.6,0.6,0.7,0.7,0.7
NL_STRIDE=10,10,10,10,10,10
NL_SKIN=0.1,0.1,0.1,0.1,0.1,0.1
... PIV
PIV ...
LABEL=Pivd3
PRECISION=1000
NLIST
REF_FILE=Ref3.pdb
PIVATOMS=3
ATOMTYPES=A,B,C
SFACTOR=0.3,0.5,1.0,0.2,0.2,0.2
SORT=1,1,1,1,1,1
SWITCH1={RATIONAL R_0=0.6 MM=12 NN=4}
SWITCH2={RATIONAL R_0=0.4 MM=10 NN=5}
SWITCH3={RATIONAL R_0=0.4 MM=10 NN=5}
SWITCH4={RATIONAL R_0=0.5 MM=12 NN=6}
SWITCH5={RATIONAL R_0=0.5 MM=12 NN=6}
SWITCH6={RATIONAL R_0=0.5 MM=12 NN=6}
NL_CUTOFF=0.8,0.6,0.6,0.7,0.7,0.7
NL_STRIDE=10,10,10,10,10,10
NL_SKIN=0.1,0.1,0.1,0.1,0.1,0.1
... PIV

PRINT ARG=Pivd1,Pivd2,Pivd3 FILE=colvar
\endplumedfile

WARNING:
Both the "CRYST" and "ATOM" lines of the PDB files must conform precisely to the official pdb format, including the width of each alphanumerical field:

\verbatim
CRYST1   31.028   36.957   23.143  89.93  92.31  89.99 P 1           1
ATOM      1  OW1 wate    1      15.630  19.750   1.520  1.00  0.00
\endverbatim

In each pdb frame, atoms must be numbered in the same order and with the same element symbol as in the input of the MD program.

The following example calculates the PIV-distances from two reference configurations Ref1.pdb and Ref2.pdb
and uses PIV-distances to define a Path Collective Variable (\ref FUNCPATHMSD) with only two references (Ref1.pdb and Ref2.pdb).
With the VOLUME keyword one scales the atom-atom distances by the cubic root of the ratio between the specified value and the box volume of the initial step of the trajectory file.

\plumedfile
PIV ...
LABEL=c1
PRECISION=1000
VOLUME=12.15
NLIST
REF_FILE=Ref1.pdb
PIVATOMS=2
ATOMTYPES=A,B
ONLYDIRECT
SFACTOR=1.0,0.2
SORT=1,1
SWITCH1={RATIONAL R_0=0.6 MM=12 NN=4}
SWITCH2={RATIONAL R_0=0.5 MM=10 NN=5}
NL_CUTOFF=1.2,1.2
NL_STRIDE=10,10
NL_SKIN=0.1,0.1
... PIV
PIV ...
LABEL=c2
PRECISION=1000
VOLUME=12.15
NLIST
REF_FILE=Ref2.pdb
PIVATOMS=2
ATOMTYPES=A,B
ONLYDIRECT
SFACTOR=1.0,0.2
SORT=1,1
SWITCH1={RATIONAL R_0=0.6 MM=12 NN=4}
SWITCH2={RATIONAL R_0=0.4 MM=10 NN=5}
NL_CUTOFF=1.2,1.2
NL_STRIDE=10,10
NL_SKIN=0.1,0.1
... PIV

p1: FUNCPATHMSD ARG=c1,c2 LAMBDA=0.180338
METAD ARG=p1.s,p1.z SIGMA=0.01,0.2 HEIGHT=0.8 PACE=500   LABEL=res
PRINT ARG=c1,c2,p1.s,p1.z,res.bias STRIDE=500  FILE=colvar FMT=%15.6f
\endplumedfile

When using PIV please cite \cite pipolo2017navigating .

(See also \ref PRINT)

*/
//+ENDPLUMEDOC

// -- SD Moved the following class declaration to header file.
//class PIV      : public Colvar
//{
//private:
//  bool pbc, serial, timer;
//  ForwardDecl<Stopwatch> stopwatch_fwd;
//  Stopwatch& stopwatch=*stopwatch_fwd;
//  // Added NL_const_size to fix solute-solvent elements as constant size
//  int updatePIV,NL_const_size;
//  size_t Nprec;
//  unsigned Natm,Nlist,NLsize;
//  double Fvol,Vol0,m_PIVdistance;
//  std::string ref_file;
//  NeighborList *nlall;
//  std::vector<SwitchingFunction> sfs;
//  std::vector<std:: vector<double> > rPIV;
//  std::vector<double> scaling,r00;
//  std::vector<double> nl_skin;
//  std::vector<double> fmass;
//  std::vector<bool> dosort;
//  std::vector<Vector> compos;
//  std::vector<string> sw;
//  std::vector<NeighborList *> nl;
//  std::vector<NeighborList *> nlcom;
//  std::vector<Vector> m_deriv;
//  // ann_deriv is the 3D array (dv(r)/dxyz) passed to the plumed core --NH
//  std::vector<std:: vector<Vector> > ann_deriv;
//  // dr_dxyz_array is the 3D array (dr/dxyz) used to build ann_deriv and ANN_sum_array --NH
//  std::vector<std:: vector<Vector> > dr_dxyz_array;
//  // ds_array is the 1D array (dv(r)/dr) of the switching function --NH
//  std::vector<double> ds_array;
//  // ANN_sum_array is the 1D array (sum dv_d/dv_n) written to an output file for use by the ANN code --NH
//  //std::vector<double> ANN_sum_array;
//  // ANN piv derivatives array written to output file for use by ANN code --SD
//  std::vector<std::vector<double>> ANN_piv_deriv;
//  // The PIV_Pair vectors record the atom IDs for the PIV elements that are passed to the VAE --NH
//  std::vector<int> PIV_Pair0;
//  std::vector<int> PIV_Pair1;
//  Tensor m_virial;
//  // adding a flag (cart2piv) for post-processing a trajectory in cartesian coordinates to a PIV representation
//  bool Svol,cross,direct,doneigh,test,CompDer,com,cart2piv;
//  int writestride;
//public:
//  static void registerKeywords( Keywords& keys );
//  explicit PIV(const ActionOptions&);
//  ~PIV();
//  // active methods:
//  virtual void calculate();
//  void checkFieldsAllowed() {}
//  // SD ANN SUM DERIVATIVE
//  std::vector<vector<double>> get_ann_sum_derivative( );
//};

PLUMED_REGISTER_ACTION(PIV,"PIV")

void PIV::registerKeywords( Keywords& keys )
{
  Colvar::registerKeywords( keys );
  keys.add("numbered","SWITCH","The switching functions parameter."
           "You should specify a Switching function for all PIV blocks."
           "Details of the various switching "
           "functions you can use are provided on \\ref switchingfunction.");
  keys.add("compulsory","PRECISION","the precision for approximating reals with integers in sorting.");
  keys.add("compulsory","REF_FILE","PDB file name that contains the \\f$i\\f$th reference structure.");
  keys.add("compulsory","PIVATOMS","Number of atoms to use for PIV.");
  keys.add("compulsory","SORT","Whether to sort or not the PIV block.");
  keys.add("compulsory","ATOMTYPES","The atom types to use for PIV.");
  keys.add("optional","SFACTOR","Scale the PIV-distance by such block-specific factor");
  keys.add("optional","VOLUME","Scale atom-atom distances by the cubic root of the cell volume. The input volume is used to scale the R_0 value of the switching function. ");
  keys.add("optional","UPDATEPIV","Frequency (in steps) at which the PIV is updated.");
  keys.addFlag("TEST",false,"Print the actual and reference PIV and exit");
  keys.addFlag("COM",false,"Use centers of mass of groups of atoms instead of atoms as specified in the Pdb file");
  keys.addFlag("ONLYCROSS",false,"Use only cross-terms (A-B, A-C, B-C, ...) in PIV");
  keys.addFlag("ONLYDIRECT",false,"Use only direct-terms (A-A, B-B, C-C, ...) in PIV");
  keys.addFlag("DERIVATIVES",false,"Activate the calculation of the PIV for every class (needed for numerical derivatives).");
  keys.addFlag("NLIST",false,"Use a neighbor list for distance calculations.");
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.addFlag("TIMER",false,"Perform timing analysis on heavy loops.");
  keys.addFlag("PIVREP",false,"Post process a trajectory from cartesian coordinates to a PIV representation.");
  keys.add("optional","NL_CONSTANT_SIZE","Fix the number of elements in all blocks to be constant. Blocks which have a total number of possible elements less than the chosen constant size will not be affected.");
  keys.add("optional","NL_CUTOFF","Neighbor lists cutoff.");
  keys.add("optional","NL_STRIDE","Update neighbor lists every NL_STRIDE steps.");
  keys.add("optional","NL_SKIN","The maximum atom displacement tolerated for the neighbor lists update.");
<<<<<<< HEAD
  // -- SD Flag for writing PIV values in a single file when using plumed driver.
  keys.addFlag("WRITEPIVTRAJ",false,"Flag to enable or disable writing PIV_representation when using plumed driver.");
  // -- SD Variables to control frequency of writing PIV values and ANN PIV derivatives during simulation.
=======
  keys.addFlag("WRITEPIVTRAJ",false,"Flag to enable or disable writing PIV_representation when using plumed driver.");
>>>>>>> 608e0a56b108e5745bff0001c31a5c0459882838
  keys.add("optional","WRITEPIVSTRIDE","STRIDE to write PIV_representation.");
  keys.add("optional","WRITEANNSTRIDE","STRUDE to write ANN_derivative files.");
  componentsAreNotOptional(keys);
  // Changing "COMPONENTS" to "default" and slightly modifying the name. Added components for ANN_SUM_DERIV
  keys.addOutputComponent("ELEMENT", "default", "Elements of the PIV block. The position in the N choose 2 interactions (i) and the neighbor in the neighbor list (j) is given as PIV-i-j.");
  //keys.addOutputComponent("ANNSUMDERIV", "default", "2D array of PIV element partial derivatives (used with ANN module).");
  keys.reset_style("SWITCH","compulsory");
}

PIV::PIV(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true),
  serial(false),
  timer(false),
  NL_const_size(0),
  updatePIV(1),
  Nprec(1000),
  Natm(1),
  Nlist(1),
  NLsize(1),
  Fvol(1.),
  Vol0(0.),
  m_PIVdistance(0.),
  rPIV(std:: vector<std:: vector<double> >(Nlist)),
  scaling(std:: vector<double>(Nlist)),
  r00(std:: vector<double>(Nlist)),
  nl_skin(std:: vector<double>(Nlist)),
  fmass(std:: vector<double>(Nlist)),
  dosort(std:: vector<bool>(Nlist)),
  compos(std:: vector<Vector>(NLsize)),
  sw(std:: vector<string>(Nlist)),
  nl(std:: vector<NeighborList *>(Nlist)),
  nlcom(std:: vector<NeighborList *>(NLsize)),
  m_deriv(std:: vector<Vector>(1)),
  dr_dxyz_array(std:: vector<std:: vector<Vector> >(1)),
  ds_array(std:: vector<double>(1)),
  //ANN_sum_array(std:: vector<double>(1)),
  ANN_piv_deriv(std:: vector<std:: vector<double>>(Nlist)),
  ann_deriv(std:: vector<std:: vector<Vector> >(1)),
  PIV_Pair0(std:: vector<int>(1)),
  PIV_Pair1(std:: vector<int>(1)),
  Svol(false),
  cross(true),
  direct(true),
  doneigh(false),
  test(false),
  CompDer(false),
  com(false),
<<<<<<< HEAD
  // SD -- local variables corresponding to user defined flags.
  writepivtraj(false),
  writepivstride(1),
  writeannstride(1),
  cart2piv(false),
  // SD -- used in prepare function.
  invalidateList(true),
  firsttime(true)
=======
  writepivtraj(false),
  writepivstride(1),
  writeannstride(1),
  cart2piv(false)
>>>>>>> 608e0a56b108e5745bff0001c31a5c0459882838
{
  log << "Starting PIV Constructor\n";

  // Precision on the real-to-integer transformation for the sorting
  parse("PRECISION",Nprec);
  if(Nprec<2) error("Precision must be => 2");

  // PBC
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
  if(pbc) {
    log << "Using Periodic Boundary Conditions\n";
  } else  {
    log << "Isolated System (NO PBC)\n";
  }

  // SERIAL/PARALLEL
  parseFlag("SERIAL",serial);
  if(serial) {
    log << "Serial PIV construction\n";
  } else     {
    log << "Parallel PIV construction\n";
  }

  // Derivatives
  parseFlag("DERIVATIVES",CompDer);
  if(CompDer) log << "Computing Derivatives\n";

  // Timing
  parseFlag("TIMER",timer);
  if(timer) {
    log << "Timing analysis\n";
    stopwatch.start();
    stopwatch.pause();
  }

  // Test
  parseFlag("TEST",test);

  // PIV Representation
  parseFlag("PIVREP",cart2piv);

  // Constant Neighbor List Size
  if(keywords.exists("NL_CONSTANT_SIZE")) {
    parse("NL_CONSTANT_SIZE",NL_const_size);
  }

  // UPDATEPIV
  if(keywords.exists("UPDATEPIV")) {
    parse("UPDATEPIV",updatePIV);
  }

  // Test
  parseFlag("COM",com);
  if(com) log << "Building PIV using COMs\n";

  // Volume Scaling
  parse("VOLUME",Vol0);
  if (Vol0>0) {
    Svol=true;
  }

  // PIV direct and cross blocks
  bool oc=false,od=false;
  parseFlag("ONLYCROSS",oc);
  parseFlag("ONLYDIRECT",od);
  if (oc&&od) {
    error("ONLYCROSS and ONLYDIRECT are incompatible options!");
  }
  if(oc) {
    direct=false;
    log << "Using only CROSS-PIV blocks\n";
  }
  if(od) {
    cross=false;
    log << "Using only DIRECT-PIV blocks\n";
  }

  // Atoms for PIV
  parse("PIVATOMS",Natm);
  std:: vector<string> atype(Natm);
  parseVector("ATOMTYPES",atype);
  //if(atype.size()!=getNumberOfArguments() && atype.size()!=0) error("not enough values for ATOMTYPES");

  // Reference PDB file
  parse("REF_FILE",ref_file);
  PDB mypdb;
  FILE* fp=fopen(ref_file.c_str(),"r");
  if (fp!=NULL) {
    log<<"Opening PDB file with reference frame: "<<ref_file.c_str()<<"\n";
    mypdb.readFromFilepointer(fp,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength());
    fclose (fp);
  } else {
    error("Error in reference PDB file");
  }

  // Build COM/Atom lists of AtomNumbers (this might be done in PBC.cpp)
  // Atomlist or Plist used to build pair lists
  std:: vector<std:: vector<AtomNumber> > Plist(Natm);
  // Atomlist used to build list of atoms for each COM
  std:: vector<std:: vector<AtomNumber> > comatm(1);
  // NLsize is the number of atoms in the pdb cell
  NLsize=mypdb.getAtomNumbers().size();
  // In the following P stands for Point (either an Atom or a COM)
  unsigned resnum=0;
  // Presind (array size: number of residues) contains the contains the residue number
  //   this is because the residue numbers may not always be ordered from 1 to resnum
  std:: vector<unsigned> Presind;
  // Build Presind
  for (unsigned i=0; i<mypdb.getAtomNumbers().size(); i++) {
    unsigned rind=mypdb.getResidueNumber(mypdb.getAtomNumbers()[i]);
    bool oldres=false;
    for (unsigned j=0; j<Presind.size(); j++) {
      if(rind==Presind[j]) {
        oldres=true;
      }
    }
    if(!oldres) {
      Presind.push_back(rind);
    }
  }
  resnum=Presind.size();

  // Pind0 is the atom/COM used in Nlists (for COM Pind0 is the first atom in the pdb belonging to that COM)
  unsigned Pind0size;
  if(com) {
    Pind0size=resnum;
  } else {
    Pind0size=NLsize;
  }
  std:: vector<unsigned> Pind0(Pind0size);
<<<<<<< HEAD
  // SD -- following resize of COM arrays to NLsize don't make sense to me. Should it be resnum? We don't use COM anyway.
=======
  // SD -- following resize of COM arrays to NLsize don't make sense. Should it be resnum?
>>>>>>> 608e0a56b108e5745bff0001c31a5c0459882838
  // If COM resize important arrays
  comatm.resize(NLsize);
  if(com) {
    nlcom.resize(NLsize);
    compos.resize(NLsize);
    fmass.resize(NLsize,0.);
  }
<<<<<<< HEAD
  // SD -- following total atoms don't make sense to me.
=======
  // SD -- following total atoms don't make sense.
>>>>>>> 608e0a56b108e5745bff0001c31a5c0459882838
  log << "Total COM/Atoms: " << Natm*resnum << " \n";
  // Build lists of Atoms/COMs for NLists
  //   comatm filled also for non_COM calculation for analysis purposes
  unsigned countIndex = 0;
  for (unsigned j=0; j<Natm; j++) {
    unsigned oind;
    for (unsigned i=0; i<Pind0.size(); i++) {
      Pind0[i]=0;
    }
    for (unsigned i=0; i<mypdb.getAtomNumbers().size(); i++) {
      // Residue/Atom AtomNumber: used to build NL for COMS/Atoms pairs.
      AtomNumber anum=mypdb.getAtomNumbers()[i];
      // ResidueName/Atomname associated to atom
      string rname=mypdb.getResidueName(anum);
      string aname=mypdb.getAtomName(anum);
      // Index associated to residue/atom: used to separate COM-lists
      unsigned rind=mypdb.getResidueNumber(anum);
      unsigned aind=anum.index();
      // This builds lists for NL
      string Pname;
      unsigned Pind;
      if(com) {
        Pname=rname;
        for(unsigned l=0; l<resnum; l++) {
          if(rind==Presind[l]) {
            Pind=l;
          }
        }
      } else {
        Pname=aname;
        Pind=aind;
      }
      if(Pname==atype[j]) {
        if(Pind0[Pind]==0) {
          // adding the atomnumber to the atom/COM list for pairs
<<<<<<< HEAD
          // SD local variable of type AtomNumber. Its value is set using countIndex.
          AtomNumber ati;
          ati.setIndex(countIndex);
          Plist[j].push_back(ati); //anum) -- SD (previously, it is same as atom number in PDB file).;
=======
          AtomNumber ati;
          ati.setIndex(countIndex);
          Plist[j].push_back(ati); //num); //anum); //(ati);
>>>>>>> 608e0a56b108e5745bff0001c31a5c0459882838
          Pind0[Pind]=aind+1;
          oind=Pind;
          countIndex += 1;
        }
        // adding the atomnumber to list of atoms for every COM/Atoms
        comatm[Pind0[Pind]-1].push_back(anum);
      }
    }
    // Output Lists
    log << "  Groups of type  " << j << ": " << Plist[j].size() << " \n";
    string gname;
    unsigned gsize;
    if(com) {
      gname=mypdb.getResidueName(comatm[Pind0[oind]-1][0]);
      gsize=comatm[Pind0[oind]-1].size();
    } else {
      gname=mypdb.getAtomName(comatm[Pind0[oind]-1][0]);
      gsize=1;
    }
    // SDlog
    //log.printf("    %6s %3s %13s %10i %6s\n", "type  ", gname.c_str(),"   containing ",gsize," atoms");
  }

  // SD This is to build the list with the atoms required for PIV.
  std:: vector<AtomNumber> listall;
  for (unsigned j=0; j<Natm; j++) {
    for (unsigned i=0; i<mypdb.getAtomNumbers().size(); i++) {
<<<<<<< HEAD
      // SD --- including only user defined atom types;
=======
      // SD --- including only user defined atom types
>>>>>>> 608e0a56b108e5745bff0001c31a5c0459882838
      AtomNumber at_num = mypdb.getAtomNumbers()[i];                                                                        
      // ResidueName/Atomname associated to atom                                                                        
      string at_name = mypdb.getAtomName(at_num);                                                                             
      if(at_name == atype[j]) {
<<<<<<< HEAD
        // -- SD listall should contain the actual atom numbers in the PDB file.
=======
>>>>>>> 608e0a56b108e5745bff0001c31a5c0459882838
        listall.push_back(at_num);
      }                                                                                                               
    }                                                                                                                 
  }    

<<<<<<< HEAD
  // SD previously, listall has all the atoms in the system.
  //for (unsigned i=0; i<mypdb.getAtomNumbers().size(); i++) {                                                          
  //  listall.push_back(mypdb.getAtomNumbers()[i]);                                                                     
  //}
=======
  //printf("countIndex: %d\n", countIndex); 

  //for (unsigned i=0; i<mypdb.getAtomNumbers().size(); i++) {                                                          
  //  listall.push_back(mypdb.getAtomNumbers()[i]);                                                                     
  //}

  //for (unsigned l=0; l<listall.size(); l++){
  //  printf("l: %d \t list: %d\t total: %d\n", l, listall[l], mypdb.getAtomNumbers().size());
  //}
  //exit();

>>>>>>> 608e0a56b108e5745bff0001c31a5c0459882838

  // PIV blocks and Neighbour Lists
  Nlist=0;
  // Direct adds the A-A ad B-B blocks (N)
  if(direct) {
    Nlist=Nlist+unsigned(Natm);
  }
  // Cross adds the A-B blocks (N*(N-1)/2)
  if(cross) {
    Nlist=Nlist+unsigned(double(Natm*(Natm-1))/2.);
  }
  // Resize vectors according to Nlist
  rPIV.resize(Nlist);

  // PIV scaled option
  scaling.resize(Nlist);
  for(unsigned j=0; j<Nlist; j++) {
    scaling[j]=1.;
  }

  if(keywords.exists("SFACTOR")) {
    parseVector("SFACTOR",scaling);
    //if(scaling.size()!=getNumberOfArguments() && scaling.size()!=0) error("not enough values for SFACTOR");
  }

  // Added STRIDE to write PIV representation and ANN sum derivatives -- SD
  if(keywords.exists("WRITEPIVTRAJ")){
      parseFlag("WRITEPIVTRAJ",writepivtraj);
  }
  if(keywords.exists("WRITEPIVSTRIDE")) { 
    parse("WRITEPIVSTRIDE",writepivstride);
  }
  if(keywords.exists("WRITEANNSTRIDE")){
    parse("WRITEANNSTRIDE",writeannstride);
  }

  // Neighbour Lists option
  parseFlag("NLIST",doneigh);
  nl.resize(Nlist);
  nl_skin.resize(Nlist);
  if(doneigh) {
    std:: vector<double> nl_cut(Nlist,0.);
    std:: vector<int> nl_st(Nlist,0);
    parseVector("NL_CUTOFF",nl_cut);
    //if(nl_cut.size()!=getNumberOfArguments() && nl_cut.size()!=0) error("not enough values for NL_CUTOFF");
    parseVector("NL_STRIDE",nl_st);
    //if(nl_st.size()!=getNumberOfArguments() && nl_st.size()!=0) error("not enough values for NL_STRIDE");
    parseVector("NL_SKIN",nl_skin);
    //if(nl_skin.size()!=getNumberOfArguments() && nl_skin.size()!=0) error("not enough values for NL_SKIN");
    for (unsigned j=0; j<Nlist; j++) {
      if(nl_cut[j]<=0.0) error("NL_CUTOFF should be explicitly specified and positive");
      if(nl_st[j]<=0) error("NL_STRIDE should be explicitly specified and positive");
      if(nl_skin[j]<=0.) error("NL_SKIN should be explicitly specified and positive");
      nl_cut[j]=nl_cut[j]+nl_skin[j];
    }
    log << "Creating Neighbor Lists \n";
<<<<<<< HEAD
    // SD -- nlall is a neighbor list created using list all. nl_cut[0] and nl_st[0] are probably not needed.
=======
    // SD -- Why nl_cut[0] and nl_st[0]? Why listall -- all atoms?
>>>>>>> 608e0a56b108e5745bff0001c31a5c0459882838
    // WARNING: is nl_cut meaningful here?
    nlall= new NeighborList(listall,true,pbc,getPbc(),comm,nl_cut[0],nl_st[0]);
    if(com) {
      //Build lists of Atoms for every COM
      for (unsigned i=0; i<compos.size(); i++) {
        // WARNING: is nl_cut meaningful here?
        nlcom[i]= new NeighborList(comatm[i],true,pbc,getPbc(),comm,nl_cut[0],nl_st[0]);
      }
    }
    unsigned ncnt=0;
    // Direct blocks AA, BB, CC, ...
    if(direct) {
      for (unsigned j=0; j<Natm; j++) {
        nl[ncnt]= new NeighborList(Plist[j],true,pbc,getPbc(),comm,nl_cut[j],nl_st[j]);
        ncnt+=1;
      }
    }

<<<<<<< HEAD
=======
    // SD -- Shouldn't the following change to account for NL constant size for each block?
>>>>>>> 608e0a56b108e5745bff0001c31a5c0459882838
    // Cross blocks AB, AC, BC, ...
    if(cross) {

      // -- SD No changes here. nl depends on Plist for each j. Plist for each j = [0, Num of atoms of type j]
      // -- SD example: if j=0 corresponds to C1, Plist[0] = [0] because there is only one C1; if is OW, it is [0, 1, 2, ... Nwaters]
      for (unsigned j=0; j<Natm; j++) {
        for (unsigned i=j+1; i<Natm; i++) {
          nl[ncnt]= new NeighborList(Plist[i],Plist[j],true,false,pbc,getPbc(),comm,nl_cut[ncnt],nl_st[ncnt]);
          //printf("\n\n\n nl[%d]: %d \n", ncnt, nl[ncnt]->getFullAtomList().size());
          ncnt+=1;
        }
      }

<<<<<<< HEAD
=======
    //printf("\n\n##########################################\n\n");
    //int tempcnt=0
    //for(int j = 0; j < Natm; j ++) {                                                                                    
    //  for(int i= j+1; i < Natm; i++) {                                                                                  
    //      if(i == Natm - 1) {                                                                                           
    //        for(int n = 0; n < NL_const_size; n++) {                                                                    
    //          string comp = "ELEMENT-" + to_string(total_count);                                                        
    //          addComponentWithDerivatives(comp);                                                                        
    //          componentIsNotPeriodic(comp);                                                                             
    //          total_count += 1;                                                                                         
    //        }                                                                                                           
    //      } else {                                                                                                      
    //        string comp = "ELEMENT-" + to_string(total_count);                                                          
    //        addComponentWithDerivatives(comp);                                                                          
    //        componentIsNotPeriodic(comp);                                                                               
    //        total_count +=1;                                                                                            
    //      }                                                                                                             
    //    printf("nl[%d]: %f\t]");
    //  }
    //} 


>>>>>>> 608e0a56b108e5745bff0001c31a5c0459882838
    }
  } else {
    log << "WARNING: Neighbor List not activated this has not been tested!!  \n";
    nlall= new NeighborList(listall,true,pbc,getPbc(),comm);
    for (unsigned j=0; j<Nlist; j++) {
      nl[j]= new NeighborList(Plist[j],Plist[j],true,true,pbc,getPbc(),comm);
    }
  }
  // Output Nlist
  log << "Total Nlists: " << Nlist << " \n";
  for (unsigned j=0; j<Nlist; j++) {
    log << "  list " << j+1 << "   size " << nl[j]->size() << " \n";
  }
  // Calculate COM masses once and for all from lists
  if(com) {
    for(unsigned j=0; j<compos.size(); j++) {
      double commass=0.;
      for(unsigned i=0; i<nlcom[j]->getFullAtomList().size(); i++) {
        unsigned andx=nlcom[j]->getFullAtomList()[i].index();
        commass+=mypdb.getOccupancy()[andx];
      }
      for(unsigned i=0; i<nlcom[j]->getFullAtomList().size(); i++) {
        unsigned andx=nlcom[j]->getFullAtomList()[i].index();
        if(commass>0.) {
          fmass[andx]=mypdb.getOccupancy()[andx]/commass;
        } else {
          fmass[andx]=1.;
        }
      }
    }
  }

  // Sorting
  dosort.resize(Nlist);
  std:: vector<int> ynsort(Nlist);
  parseVector("SORT",ynsort);
  if(cart2piv) {
    for (unsigned i=0; i<Nlist; i++) {
      if(ynsort[i]==0) {
        dosort[i]=false;
      } else {
        dosort[i]=true;
      }
    }
  } else {
    for (unsigned i=0; i<Nlist; i++) {
      if(ynsort[i]==0||CompDer) {
        dosort[i]=false;
      } else {
        dosort[i]=true;
      }
    }
  }
  //build box vectors and correct for pbc
  log << "Building the box from PDB data ... \n";
  Tensor Box=mypdb.getBoxVec();
  log << "  Done! A,B,C vectors in Cartesian space:  \n";
  log.printf("  A:  %12.6f%12.6f%12.6f\n", Box[0][0],Box[0][1],Box[0][2]);
  log.printf("  B:  %12.6f%12.6f%12.6f\n", Box[1][0],Box[1][1],Box[1][2]);
  log.printf("  C:  %12.6f%12.6f%12.6f\n", Box[2][0],Box[2][1],Box[2][2]);
  log << "Changing the PBC according to the new box \n";
  Pbc mypbc;
  mypbc.setBox(Box);
  log << "The box volume is " << mypbc.getBox().determinant() << " \n";

  //Compute scaling factor
  if(Svol) {
    Fvol=cbrt(Vol0/mypbc.getBox().determinant());
    log << "Scaling atom distances by  " << Fvol << " \n";
  } else {
    log << "Using unscaled atom distances \n";
  }

  r00.resize(Nlist);
  sw.resize(Nlist);
  for (unsigned j=0; j<Nlist; j++) {
    if( !parseNumbered( "SWITCH", j+1, sw[j] ) ) break;
  }
  if(CompDer) {
    // Set switching function parameters here only if computing derivatives
    //   now set at the beginning of the dynamics to solve the r0 issue
    log << "Switching Function Parameters \n";
    sfs.resize(Nlist);
    std::string errors;
    for (unsigned j=0; j<Nlist; j++) {
      if(Svol) {
        double r0;
        vector<string> data=Tools::getWords(sw[j]);
        data.erase(data.begin());
        Tools::parse(data,"R_0",r0);
        std::string old_r0; Tools::convert(r0,old_r0);
        r0*=Fvol;
        std::string new_r0; Tools::convert(r0,new_r0);
        std::size_t pos = sw[j].find("R_0");
        sw[j].replace(pos+4,old_r0.size(),new_r0);
      }
      sfs[j].set(sw[j],errors);
      std::string num;
      Tools::convert(j+1, num);
      if( errors.length()!=0 ) error("problem reading SWITCH" + num + " keyword : " + errors );
      r00[j]=sfs[j].get_r0();
      log << "  Swf: " << j << "  r0=" << (sfs[j].description()).c_str() << " \n";
    }
  }

  // build COMs from positions if requested
  if(com) {
    for(unsigned j=0; j<compos.size(); j++) {
      compos[j][0]=0.;
      compos[j][1]=0.;
      compos[j][2]=0.;
      for(unsigned i=0; i<nlcom[j]->getFullAtomList().size(); i++) {
        unsigned andx=nlcom[j]->getFullAtomList()[i].index();
        compos[j]+=fmass[andx]*mypdb.getPositions()[andx];
      }
    }
  }
  // build the rPIV distances (transformation and sorting is done afterwards)
  if(CompDer) {
    log << "  PIV  |  block   |     Size      |     Zeros     |     Ones      |" << " \n";
  }
  for(unsigned j=0; j<Nlist; j++) {
    for(unsigned i=0; i<nl[j]->size(); i++) {
      unsigned i0=(nl[j]->getClosePairAtomNumber(i).first).index();
      unsigned i1=(nl[j]->getClosePairAtomNumber(i).second).index();

      //printf("i0: %d \t i1: %d \n", i0, i1);

      //calculate/get COM position of centers i0 and i1
      Vector Pos0,Pos1;
      if(com) {
        //if(pbc) makeWhole();
        Pos0=compos[i0];
        Pos1=compos[i1];
      } else {
        Pos0=mypdb.getPositions()[i0];
        Pos1=mypdb.getPositions()[i1];
      }
      Vector ddist;
      if(pbc) {
        ddist=mypbc.distance(Pos0,Pos1);
      } else {
        ddist=delta(Pos0,Pos1);
      }
      double df=0.;
      // Transformation and sorting done at the first timestep to solve the r0 definition issue
      if(CompDer) {
        rPIV[j].push_back(sfs[j].calculate(ddist.modulo()*Fvol, df));
      } else {
        rPIV[j].push_back(ddist.modulo()*Fvol);
      }
    }
    if(CompDer) {
      if(dosort[j]) {
        std::sort(rPIV[j].begin(),rPIV[j].end());
      }
      int lmt0=0;
      int lmt1=0;
      for(unsigned i=0; i<rPIV[j].size(); i++) {
        if(int(rPIV[j][i]*double(Nprec-1))==0) {
          lmt0+=1;
        }
        if(int(rPIV[j][i]*double(Nprec-1))==1) {
          lmt1+=1;
        }
      }
      // SDlog
      //log.printf("       |%10i|%15zu|%15i|%15i|\n", j, rPIV[j].size(), lmt0, lmt1);
    }
  }

  checkRead();
  // Create components of PIV
  if(cart2piv) {
    // cPIV hasn't been created yet so it can't be used in the loop.
    // The loop is set up generally as N(N-1)/2. 
    // The if/else statement accounts for there being >1 elements in the solute-solvent blocks, 
    // and expects that the interaction with solvent is the last block of the each solute atom's interactions
    // Total count keeps a running tally of the elements in the entire PIV so that there will be an equal number of components.
    unsigned total_count=0;
    for(int j = 0; j < Natm; j ++) {
      for(int i= j+1; i < Natm; i++) {
          if(i == Natm - 1) {
            for(int n = 0; n < NL_const_size; n++) {
              string comp = "ELEMENT-" + to_string(total_count);
              addComponentWithDerivatives(comp); 
              componentIsNotPeriodic(comp);
              total_count += 1;
            }
          } else {
            string comp = "ELEMENT-" + to_string(total_count);
            addComponentWithDerivatives(comp); 
            componentIsNotPeriodic(comp);
            total_count +=1;
          }
      }
    }
    //printf("step:%d \t natm: %d \t Nlconstsize: %d \t total count: %d\n",getStep(), Natm, NL_const_size, total_count);
    //for(int j = 0; j < total_count; j++) {
    //  for(int i = 0; i < total_count; i++) {
    //    string comp = "ANNSUMDERIV-" + to_string(j) + "-" + to_string(i);
    //          addComponent(comp); 
    //          componentIsNotPeriodic(comp);
    //  }
    //}
    requestAtoms(nlall->getFullAtomList());

<<<<<<< HEAD
    //printf("\n\n\n Num of Atoms: %d \t total_count: %d \n\n\n", getNumberOfAtoms(), total_count);
=======
    //for (unsigned j=0; j<Natm; j++) {                                                                                 
    //  for (unsigned i=j+1; i<Natm; i++) {                                                                             
    //    nl[cnt] = nl[ncnt]->getReducedAtomList();                                                                 
    //    ncnt+=1;                                                                                                      
    //  }                                                                                                               
    //}

    printf("\n\n\n Num of Atoms: %d \t total_count: %d \n\n\n", getNumberOfAtoms(), total_count);
>>>>>>> 608e0a56b108e5745bff0001c31a5c0459882838
    //printf("\n\n\n Num of Atoms in List: %d\n\n\n", listall.size());

    ann_deriv.resize(getNumberOfAtoms());
    for (int i = 0; i < getNumberOfAtoms(); i++) {
      ann_deriv[i].resize(total_count);
    }
    ds_array.resize(total_count);
  } else {
    addValueWithDerivatives();
    requestAtoms(nlall->getFullAtomList());
    setNotPeriodic();
    // getValue()->setPeridodicity(false);
    // set size of derivative vector
    m_deriv.resize(getNumberOfAtoms());
  }
}

// The following deallocates pointers
PIV::~PIV()
{
  for (unsigned j=0; j<Nlist; j++) {
    delete nl[j];
  }
  if(com) {
    for (unsigned j=0; j<NLsize; j++) {
      delete nlcom[j];
    }
  }
  delete nlall;
}


// SD request atoms in every frame.
void PIV::prepare() {
  if(nlall->getStride()>0) {
    if(firsttime || (getStep()%nlall->getStride()==0)) {
      requestAtoms(nlall->getFullAtomList());
      invalidateList=true;
      firsttime=false;
    } else {
      requestAtoms(nlall->getReducedAtomList());
      invalidateList=false;
      if(getExchangeStep()) error("Neighbor lists should be updated on exchange steps - choose a NL_STRIDE which divides the exchange stride!");
    }
    if(getExchangeStep()) firsttime=true;
  }
}



// SD function to return ANN PIV derivatives to ANN code
std::vector<vector<double>> PIV::get_ann_sum_derivative( ) {

  std::vector<vector<double>> ann_piv_deriv_arr = ANN_piv_deriv; 

  return ann_piv_deriv_arr;

}

void PIV::calculate()
{

  // Local variables
  // The following are probably needed as static arrays
  static int prev_stp=-1;
  static int init_stp=1;
  static std:: vector<std:: vector<Vector> > prev_pos(Nlist);
  static std:: vector<std:: vector<double> > cPIV(Nlist);
  static std:: vector<std:: vector<int> > Atom0(Nlist);
  static std:: vector<std:: vector<int> > Atom1(Nlist);
  std:: vector<std:: vector<int> > A0(Nprec);
  std:: vector<std:: vector<int> > A1(Nprec);
  size_t stride=1;
  unsigned rank=0;

  if(!serial) {
    stride=comm.Get_size();
    rank=comm.Get_rank();
  } else {
    stride=1;
    rank=0;
  }

  // Transform (and sort) the rPIV before starting the dynamics
  if (((prev_stp==-1) || (init_stp==1)) &&!CompDer) {
    if(prev_stp!=-1) {init_stp=0;}
    // Calculate the volume scaling factor
    if(Svol) {
      Fvol=cbrt(Vol0/getBox().determinant());
    }
    //Set switching function parameters
    log << "\n";
    log << "REFERENCE PDB # " << prev_stp+2 << " \n";
    // Set switching function parameters here only if computing derivatives
    //   now set at the beginning of the dynamics to solve the r0 issue
    log << "Switching Function Parameters \n";
    sfs.resize(Nlist);
    std::string errors;
    for (unsigned j=0; j<Nlist; j++) {
      if(Svol) {
        double r0;
        vector<string> data=Tools::getWords(sw[j]);
        data.erase(data.begin());
        Tools::parse(data,"R_0",r0);
        std::string old_r0; Tools::convert(r0,old_r0);
        r0*=Fvol;
        std::string new_r0; Tools::convert(r0,new_r0);
        std::size_t pos = sw[j].find("R_0");
        sw[j].replace(pos+4,old_r0.size(),new_r0);
      }
      sfs[j].set(sw[j],errors);
      std::string num;
      Tools::convert(j+1, num);
      if( errors.length()!=0 ) error("problem reading SWITCH" + num + " keyword : " + errors );
      r00[j]=sfs[j].get_r0();
      log << "  Swf: " << j << "  r0=" << (sfs[j].description()).c_str() << " \n";
    }
    //Transform and sort
    log << "Building Reference PIV Vector \n";
    log << "  PIV  |  block   |     Size      |     Zeros     |     Ones      |" << " \n";
    double df=0.;
    for (unsigned j=0; j<Nlist; j++) {
      for (unsigned i=0; i<rPIV[j].size(); i++) {
        rPIV[j][i]=sfs[j].calculate(rPIV[j][i], df);
      }
      if(dosort[j]) {
        std::sort(rPIV[j].begin(),rPIV[j].end());
      }
      int lmt0=0;
      int lmt1=0;
      for(unsigned i=0; i<rPIV[j].size(); i++) {
        if(int(rPIV[j][i]*double(Nprec-1))==0) {
          lmt0+=1;
        }
        if(int(rPIV[j][i]*double(Nprec-1))==1) {
          lmt1+=1;
        }
      }
      // SDlog
      //log.printf("       |%10i|%15zu|%15i|%15i|\n", j, rPIV[j].size(), lmt0, lmt1);
    }
    log << "\n";
  }
  // Do the sorting only once per timestep to avoid building the PIV N times for N rPIV PDB structures!
  if ((getStep()>prev_stp&&getStep()%updatePIV==0)||CompDer) {
    if (CompDer) log << " Step " << getStep() << "  Computing Derivatives NON-SORTED PIV \n";
    //
    // build COMs from positions if requested
    if(com) {
      if(pbc) makeWhole();
      for(unsigned j=0; j<compos.size(); j++) {
        compos[j][0]=0.;
        compos[j][1]=0.;
        compos[j][2]=0.;
        for(unsigned i=0; i<nlcom[j]->getFullAtomList().size(); i++) {
          unsigned andx=nlcom[j]->getFullAtomList()[i].index();
          compos[j]+=fmass[andx]*getPosition(andx);
        }
      }
    }
    // update neighbor lists when an atom moves out of the Neighbor list skin
    if (doneigh) {
      bool doupdate=false;
      // For the first step build previous positions = actual positions
      if (prev_stp==-1) {
        bool docom=com;
        for (unsigned j=0; j<Nlist; j++) {
          for (unsigned i=0; i<nl[j]->getFullAtomList().size(); i++) {
            Vector Pos;
            if(docom) {
              Pos=compos[i];
            } else {
              Pos=getPosition(nl[j]->getFullAtomList()[i].index());
            }
            prev_pos[j].push_back(Pos);
          }
        }
        doupdate=true;
      }
      // Decide whether to update lists based on atom displacement, every stride
      std:: vector<std:: vector<Vector> > tmp_pos(Nlist);
      if (getStep() % nlall->getStride() ==0) {
        bool docom=com;
        for (unsigned j=0; j<Nlist; j++) {
          for (unsigned i=0; i<nl[j]->getFullAtomList().size(); i++) {
            Vector Pos;
            if(docom) {
              Pos=compos[i];
            } else {
              Pos=getPosition(nl[j]->getFullAtomList()[i].index());
            }
            tmp_pos[j].push_back(Pos);
            if (pbcDistance(tmp_pos[j][i],prev_pos[j][i]).modulo()>=nl_skin[j]) {
              doupdate=true;
            }
          }
        }
      }
      // Update Nlists if needed
      if (doupdate==true) {
        for (unsigned j=0; j<Nlist; j++) {
          for (unsigned i=0; i<nl[j]->getFullAtomList().size(); i++) {
            prev_pos[j][i]=tmp_pos[j][i];
          }
          nl[j]->update(prev_pos[j]);
          log << " Step " << getStep() << "  Neighbour lists updated " << nl[j]->size() << " \n";
        }
      }
    }
    // Calculate the volume scaling factor
    if(Svol) {
      Fvol=cbrt(Vol0/getBox().determinant());
    }
    Vector ddist;
    // Global to local variables
    bool doserial=serial;
    // Build "Nlist" PIV blocks
    for(unsigned j=0; j<Nlist; j++) {
      if(dosort[j]) {
        // from global to local variables to speedup the for loop with if statements
        bool docom=com;
        bool dopbc=pbc;
        // Vectors collecting occupancies: OrdVec one rank, OrdVecAll all ranks
        std:: vector<int> OrdVec(Nprec,0);
        cPIV[j].resize(0);
        Atom0[j].resize(0);
        Atom1[j].resize(0);
        // Building distances for the PIV vector at time t
        if(timer) stopwatch.start("1 Build cPIV");
        for(unsigned i=rank; i<nl[j]->size(); i+=stride) {
          unsigned i0=(nl[j]->getClosePairAtomNumber(i).first).index();
          unsigned i1=(nl[j]->getClosePairAtomNumber(i).second).index();
          //printf("i0: %d \t i1: %d \n", i0, i1);
          Vector Pos0,Pos1;
          if(docom) {
            Pos0=compos[i0];
            Pos1=compos[i1];
          } else {
            Pos0=getPosition(i0);
            Pos1=getPosition(i1);
          }
          if(dopbc) {
            ddist=pbcDistance(Pos0,Pos1);
          } else {
            ddist=delta(Pos0,Pos1);
          }
          double df=0.;
          //Integer sorting ... faster!
          //Transforming distances with the Switching function + real to integer transformation
          int Vint=int(sfs[j].calculate(ddist.modulo()*Fvol, df)*double(Nprec-1)+0.5);
          // SD for debugging.
          // double temp = 0;
          // printf(" Frame: %d \t j: %d \t i: %d \t r: %f \t: %f \n", getStep(), i0, i1,  ddist.modulo()*Fvol, sfs[j].calculate(ddist.modulo()*Fvol, temp));
          //
          // Enables low precision with standard PIV sizes.
          if(cart2piv) {
            if(Vint == 0) {
              Vint = 1;
            }
          }
<<<<<<< HEAD
=======

          //printf("i0: %d \t i1: %d \t Vint: %d \t ddistmodulo:%f \t df: %f \n", i0, i1, Vint, ddist.modulo(), df);
>>>>>>> 608e0a56b108e5745bff0001c31a5c0459882838

          //Integer transformed distance values as index of the Ordering Vector OrdVec
          OrdVec[Vint]+=1;
          //Keeps track of atom indices for force and virial calculations
          A0[Vint].push_back(i0);
          A1[Vint].push_back(i1);
        }
        if(timer) stopwatch.stop("1 Build cPIV");
        if(timer) stopwatch.start("2 Sort cPIV");
        if(!doserial && comm.initialized()) {
          // Vectors keeping track of the dimension and the starting-position of the rank-specific pair vector in the big pair vector.
          std:: vector<int> Vdim(stride,0);
          std:: vector<int> Vpos(stride,0);
          // Vectors collecting occupancies: OrdVec one rank, OrdVecAll all ranks
          std:: vector<int> OrdVecAll(stride*Nprec);
          // Big vectors containing all Atom indexes for every occupancy (Atom0O(Nprec,n) and Atom1O(Nprec,n) matrices in one vector)
          std:: vector<int> Atom0F;
          std:: vector<int> Atom1F;
          // Vector used to reconstruct arrays
          std:: vector<unsigned> k(stride,0);
          // Zeros might be many, this slows down a lot due to MPI communication
          // Avoid passing the zeros (i=1) for atom indices
          for(unsigned i=1; i<Nprec; i++) {
            // Building long vectors with all atom indexes for occupancies ordered from i=1 to i=Nprec-1
            // Can this be avoided ???
            Atom0F.insert(Atom0F.end(),A0[i].begin(),A0[i].end());
            Atom1F.insert(Atom1F.end(),A1[i].begin(),A1[i].end());
            A0[i].resize(0);
            A1[i].resize(0);
          }
          // Resize partial arrays to fill up for the next PIV block
          A0[0].resize(0);
          A1[0].resize(0);
          A0[Nprec-1].resize(0);
          A1[Nprec-1].resize(0);
          // Avoid passing the zeros (i=1) for atom indices
          OrdVec[0]=0;
          OrdVec[Nprec-1]=0;

          // Wait for all ranks before communication of Vectors
          comm.Barrier();

          // pass the array sizes before passing the arrays
          int dim=Atom0F.size();
          // Vdim and Vpos keep track of the dimension and the starting-position of the rank-specific pair vector in the big pair vector.
          comm.Allgather(&dim,1,&Vdim[0],1);

          // TO BE IMPROVED: the following may be done by the rank 0 (now every rank does it)
          int Fdim=0;
          for(unsigned i=1; i<stride; i++) {
            Vpos[i]=Vpos[i-1]+Vdim[i-1];
            Fdim+=Vdim[i];
          }
          Fdim+=Vdim[0];
          // build big vectors for atom pairs on all ranks for all ranks
          std:: vector<int> Atom0FAll(Fdim);
          std:: vector<int> Atom1FAll(Fdim);
          // TO BE IMPROVED: Allgathers may be substituted by gathers by proc 0
          //   Moreover vectors are gathered head-to-tail and assembled later-on in a serial step.
          // Gather the full Ordering Vector (occupancies). This is what we need to build the PIV
          comm.Allgather(&OrdVec[0],Nprec,&OrdVecAll[0],Nprec);
          // Gather the vectors of atom pairs to keep track of the idexes for the forces
          comm.Allgatherv(&Atom0F[0],Atom0F.size(),&Atom0FAll[0],&Vdim[0],&Vpos[0]);
          comm.Allgatherv(&Atom1F[0],Atom1F.size(),&Atom1FAll[0],&Vdim[0],&Vpos[0]);

          // Reconstruct the full vectors from collections of Allgathered parts (this is a serial step)
          // This is the tricky serial step, to assemble together PIV and atom-pair info from head-tail big vectors
          // Loop before on l and then on i would be better but the allgather should be modified
          // Loop on blocks
          //for(unsigned m=0;m<Nlist;m++) 
          // Loop on Ordering Vector size excluding zeros (i=1)
          if(timer) stopwatch.stop("2 Sort cPIV");
          if(timer) stopwatch.start("3 Reconstruct cPIV");
          for(unsigned i=1; i<Nprec; i++) {
            // Loop on the ranks
            for(unsigned l=0; l<stride; l++) {
              // Loop on the number of head-to-tail pieces
              for(unsigned m=0; m<OrdVecAll[i+l*Nprec]; m++) {
                // cPIV is the current PIV at time t
                cPIV[j].push_back(double(i)/double(Nprec-1));
                Atom0[j].push_back(Atom0FAll[k[l]+Vpos[l]]);
                Atom1[j].push_back(Atom1FAll[k[l]+Vpos[l]]);
                k[l]+=1;
              }
            }
          }
          if(timer) stopwatch.stop("3 Reconstruct cPIV");
        } else {
          for(unsigned i=1; i<Nprec; i++) {
            for(unsigned m=0; m<OrdVec[i]; m++) {
              cPIV[j].push_back(double(i)/double(Nprec-1));
              Atom0[j].push_back(A0[i][m]);
              Atom1[j].push_back(A1[i][m]);
            }
          }
        }
      }
    }
  }
  Vector distance;
  double dfunc=0.;
  // Calculate volume scaling factor
  if(Svol) {
    Fvol=cbrt(Vol0/getBox().determinant());
  }

  // This test may be run by specifying the TEST keyword as input, it pritnts rPIV and cPIV and quits
  if(test) {
    unsigned limit=0;
    for(unsigned j=0; j<Nlist; j++) {
      if(dosort[j]) {
        limit = cPIV[j].size();
      } else {
        limit = rPIV[j].size();
      }
      // SDlog
      //log.printf("PIV Block:  %6i %12s %6i \n", j, "      Size:", limit);
      //log.printf("%6s%6s%12s%12s%36s\n","     i","     j", "    c-PIV   ","    r-PIV   ","   i-j distance vector       ");
      for(unsigned i=0; i<limit; i++) {
        unsigned i0=0;
        unsigned i1=0;
        if(dosort[j]) {
          i0=Atom0[j][i];
          i1=Atom1[j][i];
        } else {
          i0=(nl[j]->getClosePairAtomNumber(i).first).index();
          i1=(nl[j]->getClosePairAtomNumber(i).second).index();
        }
        Vector Pos0,Pos1;
        if(com) {
          Pos0=compos[i0];
          Pos1=compos[i1];
        } else {
          Pos0=getPosition(i0);
          Pos1=getPosition(i1);
        }
        if(pbc) {
          distance=pbcDistance(Pos0,Pos1);
        } else {
          distance=delta(Pos0,Pos1);
        }
        dfunc=0.;
        double cP,rP;
        if(dosort[j]) {
          cP = cPIV[j][i];
          rP = rPIV[j][rPIV[j].size()-cPIV[j].size()+i];
        } else {
          double dm=distance.modulo();
          cP = sfs[j].calculate(dm*Fvol, dfunc);
          rP = rPIV[j][i];
        }
        // SDlog
        //log.printf("%6i%6i%12.6f%12.6f%12.6f%12.6f%12.6f\n",i0,i1,cP,rP,distance[0],distance[1],distance[2]);
      }
    }
    log.printf("This was a test, now exit \n");
    exit();
  }
  // Write out file of PIV representation for each frame of the trajectory
  if(cart2piv) {
    // open a file in append mode.

    //int STRIDE = 1000;                                                        
    FILE *piv_rep_file = NULL;
    if (getStep() % writepivstride == 0) {                                                                                      
      string piv_rep_fileName = "PIV_representation_" + to_string(getStep()) + ".dat";                                    
      piv_rep_file = fopen(piv_rep_fileName.c_str(), "w+"); 
    }
   
    FILE *piv_rep_file_traj = NULL;
    if (writepivtraj) {
      string piv_rep_fileName_traj = "PIV_representation_traj.dat";
      piv_rep_file_traj = fopen(piv_rep_fileName_traj.c_str(), "a");
    }

<<<<<<< HEAD
    // SD countLoopLimit for debugging.
=======
>>>>>>> 608e0a56b108e5745bff0001c31a5c0459882838
    int countLoopLimit = 0;
    for(unsigned j=0; j<Nlist; j++) {
      bool dosorting=dosort[j];
      unsigned limit=0;
      // Set limit to size of PIV block. Solute-solute blocks will be 1 (for tetracosane) 
      // and much larger for solute-solvent blocks (likely hundreds)
      limit = cPIV[j].size();
      //printf("limit: %d\t", limit);
      // Allow for non-constant block sizes if desired
      if(NL_const_size > 0) {
        // Solute-solvent blocks have more neighbors than necessary so that padding is not necessary.
        // The solute-solvent blocks are already sorted so the last elements of the block (size NL_const_size) are the desired interactions to include
        // i.e. the closest solute-solvent interactions. This is irrelevant for the solute-solute interactions. 
        int start_val=0;
        // This sets the start value to be NL_const_size away from the end of the sorted block to choose the desired interactions.
        if(limit > NL_const_size) {
          start_val = limit - NL_const_size;
        }
        if (writepivtraj) {
<<<<<<< HEAD
          for(unsigned i=start_val; i<limit; i++) {
            fprintf(piv_rep_file_traj, "%8.6f\t", cPIV[j][i]);
          }
        }
        if ( getStep() % writepivstride == 0) {
          for(unsigned i=start_val; i<limit; i++) {
=======
          for(unsigned i=start_val; i<limit; i++) {
            fprintf(piv_rep_file_traj, "%8.6f\t", cPIV[j][i]);
          }
        }
        if ( getStep() % writepivstride == 0) {
          for(unsigned i=start_val; i<limit; i++) {
>>>>>>> 608e0a56b108e5745bff0001c31a5c0459882838
            fprintf(piv_rep_file, "%8.6f\t", cPIV[j][i]);
            countLoopLimit += 1;
          }
        }
      } else {
        // Prints out in the same PIV block element format as TEST
        if (writepivtraj) {
          for(unsigned i=0; i<limit; i++) {
            fprintf(piv_rep_file_traj, "%8.6f\t", cPIV[j][i]);
          } 
        }
        if ( getStep() % writepivstride == 0) {
          for(unsigned i=0; i<limit; i++) {
            fprintf(piv_rep_file, "%8.6f\t", cPIV[j][i]);
            countLoopLimit += 1;
          }
        }
      }
    }
<<<<<<< HEAD
    // SD for debugging
=======
>>>>>>> 608e0a56b108e5745bff0001c31a5c0459882838
    //printf("\nEnd limit: %d\n", countLoopLimit);
    if (writepivtraj) {
      fprintf(piv_rep_file_traj, "\n#END OF FRAME\n");
      fclose(piv_rep_file_traj);
    }
    if ( getStep() % writepivstride == 0) {
      fprintf(piv_rep_file, "\n#END OF FRAME: %d \n", getStep());
      fclose(piv_rep_file);
    }
  }

  if(timer) stopwatch.start("4 Build For Derivatives");
  // non-global variables Nder and Scalevol defined to speedup if structures in cycles
  bool Nder=CompDer;
  bool Scalevol=Svol;
  // Build derivatives for PIVREP differently because indexing is slightly different (only taking some of the solute-solvent interactions)
  if(cart2piv) {
    for(unsigned j=0; j<ann_deriv.size(); j++) {
      for(unsigned i=0; i<ann_deriv[j].size(); i++) {
        for(unsigned k=0; k<3; k++) {ann_deriv[j][i][k]=0.;}
      }
    }
    for(unsigned j=0; j<3; j++) {
      for(unsigned k=0; k<3; k++) {
        m_virial[j][k]=0.;
      }
    }
    // resize vectors to the appropriate sizes and set starting values to zero --NH
    //ds_array.resize(ann_deriv[0].size());
    //ANN_sum_array.resize(ds_array.size());
    ANN_piv_deriv.resize(ds_array.size());
    //for(unsigned j=0; j<ANN_sum_array.size(); j++) {
    //  ANN_sum_array[j] = 0.;
    //}
    for(unsigned j=0; j<ANN_piv_deriv.size(); j++) {
        ANN_piv_deriv[j].resize(ds_array.size());
    }
    dr_dxyz_array.resize(ann_deriv.size());
    PIV_Pair0.resize(ds_array.size());
    PIV_Pair1.resize(ds_array.size());
    for(unsigned j=0; j<dr_dxyz_array.size(); j++) {
      dr_dxyz_array[j].resize(ds_array.size());
    }
    for(unsigned j=0; j<dr_dxyz_array.size(); j++) {
      for(unsigned i=0; i<dr_dxyz_array[j].size(); i++) {
        for(unsigned k=0; k<3; k++) {
          dr_dxyz_array[j][i][k] = 0.;
        }
      }
    }
    //unsigned countLoop=0;
    unsigned PIV_element=0;
    for(unsigned j=0; j<Nlist; j++) {
      // Same loop as used/described above
      unsigned limit=0;
      limit = cPIV[j].size();
      int start_val=0;
      if(limit > NL_const_size) {
        start_val = limit - NL_const_size;
      }
      for(unsigned i=start_val; i<limit; i++) {
        unsigned i0=0;
        unsigned i1=0;
        // i0 and i1 take scalar values, so Atom0[j][i]/Atom1[j][i] must not be xyz coordinates. Atom0 and Atom1 seem
        // to be lists that index atoms within the neighbor list
        i0=Atom0[j][i];
        // Record the atom IDs for the PIV elements of interest --NH
        PIV_Pair0[PIV_element] = i0;
        i1=Atom1[j][i];
        PIV_Pair1[PIV_element] = i1;
        // Pos0 and Pos1 seem to be 1x3 vectors that hold the xyz coordinates of the indexed atoms
        Vector Pos0,Pos1;
        Pos0=getPosition(i0);
        Pos1=getPosition(i1);
        // distance is also a 1x3 vector of the xyz distances between the two atoms after consideration of the pbc
        distance=pbcDistance(Pos0,Pos1);
        dfunc=0.;
        // dm seems to be scalar value that is the magnitude of the distance between the atoms
        double dm=distance.modulo();
        // sfs[j] is the parameters for the switching function, which can be chosen to be different for different blocks
        // In this case, all blocks use the same switching function so all sfs[j] are the same function.
        // Used with .calculate(dm*Fvol, dfunc), the PIV element value is returned and the derivative stored in dfunc
        double tPIV = sfs[j].calculate(dm*Fvol, dfunc);
        //double tmp=0.;
        double ds_element=0.;
        // In our case, Fvol is 1 and so is scaling[j], so tmp is really 2*s(r)*derivative_of_s(r)
        //tmp = 2.*scaling[j]*tPIV*Fvol*Fvol*dfunc;
        // Create the ds_array one element at a time --NH
        ds_element = scaling[j]*Fvol*Fvol*dfunc*dm;
        ds_array[PIV_element] = ds_element;
<<<<<<< HEAD


=======
>>>>>>> 608e0a56b108e5745bff0001c31a5c0459882838
        //printf("ds_array: %8.10f\n", ds_array[PIV_element]);
        //countLoop += 1;
        // Create 1x3 vector of (dr/dx,dr/dy,dr/dz) --NH
        Vector dr_dcoord = distance/dm;
        // SD for debugging
        //double temp = 0;
        //printf(" Frame: %d \t j: %d \t i: %d \t r: %f \t PIV: %f \n", getStep(), i0, i1,  dm*Fvol, sfs[j].calculate(dm*Fvol, temp));
        //printf(" Frame: %d \t j: %d \t i: %d \t PIVELEMENT: %d \t dPELE/dx: %f \t dPELE/dy: %f \t dPELE/dz: %f \n", getStep(), i0, i1, PIV_element, ds_element*dr_dcoord[0], ds_element*dr_dcoord[1], ds_element*dr_dcoord[2]);
        //printf(" Frame: %d \t j: %d \t i: %d \t PIVELEMENT: %d \t dr_dcoordx: %f \t dr_dcoordy: %f \t dr_dcoordz: %f \n", getStep(), i0, i1, PIV_element, dr_dcoord[0], dr_dcoord[1], dr_dcoord[2]);
        
        //Vector tmpder = ds*distance;
        // the xyz components of the distance between atoms is scaled by tmp and added or subtracted to reflect
        // that distance is calculated as Pos1 - Pos0
        // Calculate ann_deriv values for the current PIV element in the loop --NH
        ann_deriv[i0][PIV_element] = -ds_element*dr_dcoord;
        ann_deriv[i1][PIV_element] =  ds_element*dr_dcoord;
        // Record dr/dxyz values for ANN_sum_array calculation later in code --NH 
        dr_dxyz_array[i0][PIV_element] = -dr_dcoord;
        dr_dxyz_array[i1][PIV_element] =  dr_dcoord;
        // This m_virial is likely not correct but has been included in case it is necessary to test the code --NH
        m_virial    -= ds_element*Tensor(distance,distance); // Question
        PIV_element += 1;
        //fprintf(atom0_file, "%8u\n", i0);
        //fprintf(atom1_file, "%8u\n", i1);
      }
    }
    // The file "dri_drj_values.dat" was used for debugging and has since been commented out --NH
    //printf("Count loop: %d \t ds_array size: %d \n", countLoop, ds_array.size());
    
    double dri_drj = 0.;
    //FILE *dri_drj_file = NULL;
    //dri_drj_file = fopen("dri_drj_values.dat", "a");

    //This loops over the two PIV element sets (dv_d and dv_n) --NH
    //for(unsigned j=0; j<ANN_sum_array.size(); j++) 
    for(unsigned j=0; j<ANN_piv_deriv.size(); j++){
      unsigned i0_j = PIV_Pair0[j];
      unsigned i1_j = PIV_Pair1[j];
      //for(unsigned i=0; i<ANN_sum_array.size(); i++) 
        for(unsigned i=0; i<ANN_piv_deriv[j].size(); i++){
        double dri_drjalpha=0.;
        double dri_drjbeta=0.;
        for(unsigned k=0; k<3; k++) {
          // This is where the vectors are summed into a scalar --NH
          // We will likely need to address the possibility of a zero term in the denominator of either term --NH
          if ( dr_dxyz_array[i0_j][j][k] != 0.0 ){ 
            dri_drjalpha +=(double) dr_dxyz_array[i0_j][i][k] / (double) dr_dxyz_array[i0_j][j][k];
          }
          if ( dr_dxyz_array[i1_j][j][k] != 0.0) {
            dri_drjbeta +=(double) dr_dxyz_array[i1_j][i][k] / (double) dr_dxyz_array[i1_j][j][k];
          }
        }  
        dri_drj = dri_drjalpha + dri_drjbeta;
        //fprintf(dri_drj_file, "%8.6f\n", dri_drj);
        //printf("ds_array: %8.10f\n", ds_array[j]);
        // Calculate ANN_sum_array from sub-arrays --NH 
        //ANN_sum_array[j] += ds_array[i]*dri_drj/ds_array[j];
        ANN_piv_deriv[j][i] = ds_array[i]*dri_drj/ds_array[j];
        dri_drj=0.;
      }
    }
    //fprintf(dri_drj_file, "END OF FRAME\n");
    //fclose(dri_drj_file);
    // Output the values of ANN_sum_array to a file to be read-in by the ANN module --NH
    //FILE *ANN_sum_file = NULL;
    //ANN_sum_file = fopen("ANN_deriv_sum.dat", "a");
    //for(unsigned j=0; j<ANN_sum_array.size(); j++) {
    //  fprintf(ANN_sum_file, "%8.6f\n", ANN_sum_array[j]);
    //}
    //fprintf(ANN_sum_file, "END OF FRAME\n");
    //fclose(ANN_sum_file);
    //int STRIDE = 1000;
    if (getStep() % writeannstride == 0) {
      FILE *ANN_deriv_file = NULL;
      string ANN_deriv_fileName = "ANN_deriv_file_" + to_string(getStep()) + ".dat";
      ANN_deriv_file = fopen(ANN_deriv_fileName.c_str(), "w+"); // Question: Should this be w+; a works for trajectories.
      for(unsigned j=0; j<ANN_piv_deriv.size(); j++){
        for(unsigned i=0; i<ANN_piv_deriv[j].size(); i++){
          fprintf(ANN_deriv_file, "%8.6f\t", ANN_piv_deriv[j][i]);
        }
        fprintf(ANN_deriv_file, "\n");
      }
      fprintf(ANN_deriv_file, "END OF FRAME: %d \n", getStep());
      fclose(ANN_deriv_file);
    }
    // SD for debugging get_ann_sum_derivative function.
    //get_ann_sum_derivative( ); //ANN_piv_deriv);
    //for (int j = 0; j < ANN_piv_deriv.size(); j++) {
    //  for(int i = 0; i < ANN_piv_deriv[j].size(); i++) {
    //    int comp_count;
    //    comp_count = int(ANN_piv_deriv.size())*j + i;
    //    string comp = "ANNSUMDERIV-" + to_string(j) + "-" + to_string(i);
    //    Value* valueNew=getPntrToComponent(comp);
    //    valueNew -> set(ANN_piv_deriv[j][i]);
    //  }
    //}
    
    //fprintf(atom0_file, "END OF FRAME\n");
    //fclose(atom0_file);
    //fprintf(atom1_file, "END OF FRAME\n");
    //fclose(atom1_file);
    
    log.printf("cPIV size: %10d\n", cPIV.size());
    if (!serial && comm.initialized() ) {
      int count = 0;
      for(unsigned j=0; j<Nlist; j++) {
          for(unsigned i=0; i<cPIV[j].size(); i++) {
              count += 1;
          }
      }
      
      comm.Barrier();
      // SD -- This probably works because cPIV[j] size is variable for each j.
      for (unsigned j=0; j< Nlist; j++) {
        for (unsigned k=0; k<cPIV[j].size(); k++) {
          comm.Sum(cPIV[j][k]);
          cPIV[j][k] /= comm.Get_size();
        }
      }
      //comm.Sum(&cPIV[0][0], count);

      // SD -- This probably works because comm.Sum cannot handle 3D vectors.
      if(!ann_deriv.empty()) {
        for (unsigned i=0;  i<ann_deriv.size(); i++) {
          for (unsigned j=0; j<ann_deriv[j].size(); j++) {
            for (unsigned k=0; k<3; k++) {
              comm.Sum(ann_deriv[i][j][k]);
              ann_deriv[i][j][k] /= comm.Get_size();
            }
          }
        }
        //comm.Sum(&ann_deriv[0][0][0], 3*ann_deriv.size()*ann_deriv[0].size());
      }
      // SD -- this is probably not needed.
      comm.Sum(&m_virial[0][0],9);
    }
  } else {
    if(getStep()%updatePIV==0) {
      // set to zero PIVdistance, derivatives and virial when they are calculated
      for(unsigned j=0; j<m_deriv.size(); j++) {
        for(unsigned k=0; k<3; k++) {m_deriv[j][k]=0.;}
      }
      for(unsigned j=0; j<3; j++) {
        for(unsigned k=0; k<3; k++) {
          m_virial[j][k]=0.;
        }
      }
      m_PIVdistance=0.;
      // Re-compute atomic distances for derivatives and compute PIV-PIV distance
      for(unsigned j=0; j<Nlist; j++) {
        unsigned limit=0;
        // dosorting definition is to speedup if structure in cycles with non-global variables
        bool dosorting=dosort[j];
        bool docom=com;
        bool dopbc=pbc;
        if(dosorting) {
          limit = cPIV[j].size();
        } else {
          limit = rPIV[j].size();
        }
        for(unsigned i=rank; i<limit; i+=stride) {
          unsigned i0=0;
          unsigned i1=0;
          if(dosorting) {
            i0=Atom0[j][i];
            i1=Atom1[j][i];
          } else {
            i0=(nl[j]->getClosePairAtomNumber(i).first).index();
            i1=(nl[j]->getClosePairAtomNumber(i).second).index();
          }
          Vector Pos0,Pos1;
          if(docom) {
            Pos0=compos[i0];
            Pos1=compos[i1];
          } else {
            Pos0=getPosition(i0);
            Pos1=getPosition(i1);
          }
          if(dopbc) {
            distance=pbcDistance(Pos0,Pos1);
          } else {
            distance=delta(Pos0,Pos1);
          }
          dfunc=0.;
          // this is needed for dfunc and dervatives
          double dm=distance.modulo();
          double tPIV = sfs[j].calculate(dm*Fvol, dfunc);
          // PIV distance
          double coord=0.;
          if(!dosorting||Nder) {
            coord = tPIV - rPIV[j][i];
          } else {
            coord = cPIV[j][i] - rPIV[j][rPIV[j].size()-cPIV[j].size()+i];
          }
          // Calculate derivatives, virial, and variable=sum_j (scaling[j] *(cPIV-rPIV)_j^2)
          // WARNING: dfunc=dswf/(Fvol*dm)  (this may change in future Plumed versions)
          double tmp=0.;
          tmp = 2.*scaling[j]*coord*Fvol*Fvol*dfunc;
          Vector tmpder = tmp*distance;
          // 0.5*(x_i-x_k)*f_ik         (force on atom k due to atom i)
          if(docom) {
            Vector dist;
            for(unsigned k=0; k<nlcom[i0]->getFullAtomList().size(); k++) {
              unsigned x0=nlcom[i0]->getFullAtomList()[k].index();
              m_deriv[x0] -= tmpder*fmass[x0];
              for(unsigned l=0; l<3; l++) {
                dist[l]=0.;
              }
              Vector P0=getPosition(x0);
              for(unsigned l=0; l<nlcom[i0]->getFullAtomList().size(); l++) {
                unsigned x1=nlcom[i0]->getFullAtomList()[l].index();
                Vector P1=getPosition(x1);
                if(dopbc) {
                  dist+=pbcDistance(P0,P1);
                } else {
                  dist+=delta(P0,P1);
                }
              }
              for(unsigned l=0; l<nlcom[i1]->getFullAtomList().size(); l++) {
                unsigned x1=nlcom[i1]->getFullAtomList()[l].index();
                Vector P1=getPosition(x1);
                if(dopbc) {
                  dist+=pbcDistance(P0,P1);
                } else {
                  dist+=delta(P0,P1);
                }
              }
              m_virial    -= 0.25*fmass[x0]*Tensor(dist,tmpder);
            }
            for(unsigned k=0; k<nlcom[i1]->getFullAtomList().size(); k++) {
              unsigned x1=nlcom[i1]->getFullAtomList()[k].index();
              m_deriv[x1] += tmpder*fmass[x1];
              for(unsigned l=0; l<3; l++) {
                dist[l]=0.;
              }
              Vector P1=getPosition(x1);
              for(unsigned l=0; l<nlcom[i1]->getFullAtomList().size(); l++) {
                unsigned x0=nlcom[i1]->getFullAtomList()[l].index();
                Vector P0=getPosition(x0);
                if(dopbc) {
                  dist+=pbcDistance(P1,P0);
                } else {
                  dist+=delta(P1,P0);
                }
              }
              for(unsigned l=0; l<nlcom[i0]->getFullAtomList().size(); l++) {
                unsigned x0=nlcom[i0]->getFullAtomList()[l].index();
                Vector P0=getPosition(x0);
                if(dopbc) {
                  dist+=pbcDistance(P1,P0);
                } else {
                  dist+=delta(P1,P0);
                }
              }
              m_virial    += 0.25*fmass[x1]*Tensor(dist,tmpder);
            }
          } else {
            m_deriv[i0] -= tmpder;
            m_deriv[i1] += tmpder;
            m_virial    -= tmp*Tensor(distance,distance);
          }
          if(Scalevol) {
            m_virial+=1./3.*tmp*dm*dm*Tensor::identity();
          }
          m_PIVdistance    += scaling[j]*coord*coord;
        }
      }

      if (!serial && comm.initialized()) {
        comm.Barrier();
        comm.Sum(&m_PIVdistance,1);
        if(!m_deriv.empty()) comm.Sum(&m_deriv[0][0],3*m_deriv.size());
        comm.Sum(&m_virial[0][0],9);
      }
    }
  }
  prev_stp=getStep();

  //Timing
  if(timer) stopwatch.stop("4 Build For Derivatives");
  if(timer) {
    log.printf("Timings for action %s with label %s \n", getName().c_str(), getLabel().c_str() );
    log<<stopwatch;
  }
  //log.printf("m_deriv size: %10d\n", m_deriv.size());
  if(cart2piv) {
    //FILE *ann_deriv_file = NULL;
    //ann_deriv_file = fopen("PIV2ANN_Derivatives.dat", "a");
    unsigned total_count=0;
    for (int j = 0; j < Nlist; j++) {
      unsigned limit=0;
      limit = cPIV[j].size();
      int start_val=0;
      if(limit > NL_const_size) {
        start_val = limit - NL_const_size;
      }
      for (int i = start_val; i < limit; i++) {
        string comp = "ELEMENT-" + to_string(total_count);
        Value* valueNew=getPntrToComponent(comp);
        valueNew -> set(cPIV[j][i]);
        // Pass the 3D array to the plumed core --NH
        // A 2D array is passed for each PIV element (component) --NH
        for(unsigned k=0; k<ann_deriv.size(); k++) {
          // Question - Can we exclude hydrogen atoms so that we can optimize as some the derivatives will be zeros?
          setAtomsDerivatives(valueNew, k, ann_deriv[k][total_count]);
        }
        //setAtomsDerivatives(valueNew, total_count, m_deriv[total_count]);
        total_count += 1;
      }
    }
    //fprintf(ann_deriv_file, "END OF FRAME\n");
    //fclose(ann_deriv_file);
    
    //setBoxDerivatives  (m_virial); Question -- Probably not required.
    
    //FILE *atom0_file = NULL;
    //atom0_file = fopen("Atom0.dat", "a");
    //FILE *atom1_file = NULL;
    //atom1_file = fopen("Atom1.dat", "a");
    //for(int j=0; j < Nlist; j++) {
    //  unsigned limit=0;
    //  limit = cPIV[j].size();
    //  int start_val=0;
    //  if(limit > NL_const_size) {
    //    start_val = limit - NL_const_size;
    //  }
    //  for(unsigned i=start_val; i<limit; i++) {
    //    fprintf(atom0_file, "%8u\n", Atom0[j][i]);
    //    fprintf(atom1_file, "%8u\n", Atom1[j][i]);
    //  }
    //}
    //fprintf(atom0_file, "END OF FRAME\n");
    //fclose(atom0_file);
    //fprintf(atom1_file, "END OF FRAME\n");
    //fclose(atom1_file);
  } else {
    // Update derivatives, virial, and variable (PIV-distance^2)
    for(unsigned i=0; i<m_deriv.size(); ++i) setAtomsDerivatives(i,m_deriv[i]);
    setValue           (m_PIVdistance);
    setBoxDerivatives  (m_virial);
  }
}
//Close Namespaces at the very beginning
}
}
