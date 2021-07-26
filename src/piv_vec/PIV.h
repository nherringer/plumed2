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
#ifndef __PLUMED_piv_vec_PIV_h
#define __PLUMED_piv_vec_PIV_h


//#include "colvar/Colvar.h"
//#include "colvar/ActionRegister.h"
//#include "core/PlumedMain.h"
//#include "core/ActionWithVirtualAtom.h"
//#include "tools/NeighborList.h"
//#include "tools/SwitchingFunction.h"
//#include "tools/PDB.h"
//#include "tools/Pbc.h"
//#include "tools/Stopwatch.h"
//
//#include <string>
//#include <cmath>
//#include <iostream>
//#include <stdio.h>

using namespace std;

namespace PLMD {
namespace piv {
// Ideally core/Colvar.h should be moved to this directory and Colvar should stay in namespace PLMD::Sasa
// With this trick, PLMD::Colvar is visible as PLMD::Sasa::Colvar
using PLMD::Colvar;

class PIV      : public Colvar
{
private:
  bool pbc, serial, timer;
  ForwardDecl<Stopwatch> stopwatch_fwd;
  Stopwatch& stopwatch=*stopwatch_fwd;
  // Added NL_const_size to fix solute-solvent elements as constant size
  int updatePIV,NL_const_size;
  size_t Nprec;
  unsigned Natm,Nlist,NLsize;
  double Fvol,Vol0,m_PIVdistance;
  std::string ref_file;
  NeighborList *nlall;
  std::vector<SwitchingFunction> sfs;
  std::vector<std:: vector<double> > rPIV;
  std::vector<double> scaling,r00;
  std::vector<double> nl_skin;
  std::vector<double> fmass;
  std::vector<bool> dosort;
  std::vector<Vector> compos;
  std::vector<string> sw;
  std::vector<NeighborList *> nl;
  std::vector<NeighborList *> nlcom;
  std::vector<Vector> m_deriv;
  // ann_deriv is the 3D array (dv(r)/dxyz) passed to the plumed core --NH
  std::vector<std:: vector<Vector> > ann_deriv;
  // dr_dxyz_array is the 3D array (dr/dxyz) used to build ann_deriv and ANN_sum_array --NH
  std::vector<std:: vector<Vector> > dr_dxyz_array;
  // ds_array is the 1D array (dv(r)/dr) of the switching function --NH
  std::vector<double> ds_array;
  // ANN_sum_array is the 1D array (sum dv_d/dv_n) written to an output file for use by the ANN code --NH
  //std::vector<double> ANN_sum_array;
  // ANN piv derivatives array written to output file for use by ANN code --SD
  std::vector<std::vector<double>> ANN_piv_deriv;
  // The PIV_Pair vectors record the atom IDs for the PIV elements that are passed to the VAE --NH
  std::vector<int> PIV_Pair0;
  std::vector<int> PIV_Pair1;
  Tensor m_virial;
  // adding a flag (cart2piv) for post-processing a trajectory in cartesian coordinates to a PIV representation
  bool Svol,cross,direct,doneigh,test,CompDer,com,cart2piv;
  bool writepivtraj;
  int writepivstride, writeannstride;
public:
  static void registerKeywords( Keywords& keys );                                                                       
  explicit PIV(const ActionOptions&); //ao):
//      PLUMED_COLVAR_INIT(ao),                                                                                               
//  pbc(true),                                                                                                            
//  serial(false),                                                                                                        
//  timer(false),                                                                                                                      
//  NL_const_size(0),                                                                                                     
//  updatePIV(1),                                                                                                         
//  Nprec(1000),                                                                                                          
//  Natm(1),                                                                                                              
//  Nlist(1),                                                                                                             
//  NLsize(1),                                                                                                            
//  Fvol(1.),                                                                                                             
//  Vol0(0.),                                                                                                             
//  m_PIVdistance(0.),                                                                                                    
//  rPIV(std:: vector<std:: vector<double> >(Nlist)),                                                                     
//  scaling(std:: vector<double>(Nlist)),                                                                                 
//  r00(std:: vector<double>(Nlist)),                                                                                     
//  nl_skin(std:: vector<double>(Nlist)),                                                                                 
//  fmass(std:: vector<double>(Nlist)),                                                                                   
//  dosort(std:: vector<bool>(Nlist)),                                                                                    
//  compos(std:: vector<Vector>(NLsize)),                                                                                 
//  sw(std:: vector<string>(Nlist)),                                                                                      
//  nl(std:: vector<NeighborList *>(Nlist)),                                                                              
//  nlcom(std:: vector<NeighborList *>(NLsize)),                                                                          
//  m_deriv(std:: vector<Vector>(1)),                                                                                     
//  dr_dxyz_array(std:: vector<std:: vector<Vector> >(1)),                                                                
//  ds_array(std:: vector<double>(1)),                                                                                    
//  //ANN_sum_array(std:: vector<double>(1)),                                                                             
//  ANN_piv_deriv(std:: vector<std:: vector<double>>(Nlist)),                                                             
//  ann_deriv(std:: vector<std:: vector<Vector> >(1)),                                                                    
//  PIV_Pair0(std:: vector<int>(1)),                                                                                      
//  PIV_Pair1(std:: vector<int>(1)),                                                                                      
//  Svol(false),                                                                                                          
//  cross(true),                                                                                                          
//  direct(true),                                                                                                         
//  doneigh(false),                                                                                                       
//  test(false),                                                                                                          
//  CompDer(false),                                                                                                       
//  com(false),                                                                                                           
//  writestride(1),                                                                                                       
//  cart2piv(false) ;                                                                                   
  ~PIV();                                                                                                               
  // active methods:                                                                                                    
  virtual void calculate();
  void checkFieldsAllowed() {}                                                                                           
  // SD ANN SUM DERIVATIVE                                                                                              
  std::vector<vector<double>> get_ann_sum_derivative(); //{  // {; //; // vector<vector<double>>&ann_piv_deriv_arr );

   //std::vector<vector<double>> ann_piv_deriv_arr = ANN_piv_deriv;                                                        
                                                                                                                        
    //for(unsigned j=0; j<ann_piv_deriv_arr.size(); j++){                                                                   
    //  for(unsigned i=0; i<ann_piv_deriv_arr[j].size(); i++){                                                              
    //     printf("%8.6f\t", ann_piv_deriv_arr[j][i]);                                                                    
    //  }                                                                                                                 
    //  printf("\n");                                                                                                     
    //}                                                                                                                   
    //ann_piv_deriv_arr = ANN_piv_deriv;                                                                                  
                                                                                                                        
 //   return ANN_piv_deriv;  
 //}

};

}
}

#endif
