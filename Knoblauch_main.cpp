//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//                    Knoblauch
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
#include<iostream>
#include<fstream>
#include<string>
#include"Version/Name_Info.h"
#include"Version/Git_Version.h"
#include"Grid1D/Block_Uniform_1D.h"
#include"Grid2D/Block_2D.h"
#include"Grid2D/Write_Solution_VTK.h"
#include"Grid2D/Write_Output_VTK.h"
#include"Grid2D/Make_Cartesian_Block_2D.h"
#include"Grid2D/Make_Cylinder_Block_2D.h"
//#include"Discretizations2D/Solution_Block_2D.h"
//#include"Discretizations2D/First_Order_FV_Block_2D.h"
//#include"Discretizations2D/Second_Order_FV_Cartesian_Block_2D.h"
//#include"Discretizations2D/Explicit_Euler.h"
//#include"Discretizations2D/Predictor_Corrector.h"
//#include"Discretizations2D/PI_Predictor_Corrector.h"
//#include"Euler2D/Euler2D_cState.h"
//#include"Gaussian2D/Gaussian2D_cState.h"
//#include"Particule2D/Particle_cState.h"
//#include"Particule2D/Particle2_cState.h"
//#include"Initial_Conditions/Shockbox.h"
//#include"Initial_Conditions/Vacuum.h"
//#include"Initial_Conditions/Particle_Families.h"
//#include"Discretizations2D/Write_Solution_VTK.h"
//#include"Discretizations2D/HLLE.h"
//#include"Discretizations2D/Saurel.h"
//#include"Discretizations2D/Constant_Extrapolation.h"
//#include"Discretizations2D/Particle_Beams.h"
//#include"Physics/Calorically_Perfect_Gas.h"
//#include"Physics/Common_Calorically_Perfect_Gases.h"
#include"Input/Input_Parameter_Tree.h"
#include"Grid2D/Grid_Tree/Refine_Block2D_I.h"
#include"Grid2D/Grid_Tree/Refine_Block2D_J.h"
#include"Grid2D/Grid_Tree/Refine_Block2D_IJ.h"
#include"Grid2D/Grid_Tree/Grid_Tree_Node_Refine_I.h"
#include"Grid2D/Solution_Block_2D.h"
#include"Grid2D/Grid_Tree/Tree_Node.h"
#include<fenv.h>

using namespace knoblauch;

///////////////////////////////////////////////////////////////
// main
///////////////////////////////////////////////////////////////
int main() {

  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  std::cout << code_header() << std::endl;

  Solution_Block_2D grid = make_cylinder_block_2D(1.0, 2.0, 11, 11);
  //Solution_Block_2D grid = make_cartesian_block_2D( -1.0,1.0,11,-1.0,1.0,11);

  //grid.give_value_test();

  error_norms non_ref = grid.compute_errors();
  std::cout<<"Before refined:"<<std::endl;
  std::cout<<"L1    = "<<non_ref.L1<<std::endl;
  std::cout<<"L2    = "<<non_ref.L2<<std::endl;
  std::cout<<"L_inf = "<<non_ref.L_inf<<std::endl;
  std::cout<<"number_of_cells="<<non_ref.num_cells<<std::endl<<std::endl;
  std::cout<<"////////////////////////////"<<std::endl<<std::endl;;
  
  grid.give_value_test();
  write_output_vtk(grid, "grid_o.vts");
  // write_solution_vtk(grid, "solution_grid.vts");

  //Block2D_Refined_I  refined_i  = Refine_Block2D_In_I(grid);
  //Block2D_Refined_J  refined_j  = Refine_Block2D_In_J(grid);
  Block2D_Refined_IJ refined_ij = Refine_Block2D_In_IJ(grid);

  //write_output_vtk(*refined_i.i_pos,"ipos.vts");
  //write_output_vtk(*refined_i.i_neg,"ineg.vts");

  //write_output_vtk(*refined_j.j_pos,"jpos.vts");
  //write_output_vtk(*refined_j.j_neg,"jneg.vts");

  write_output_vtk(*refined_ij.i_neg_j_neg,"injn.vts");
  write_output_vtk(*refined_ij.i_pos_j_neg,"ipjn.vts");
  write_output_vtk(*refined_ij.i_neg_j_pos,"injp.vts");
  write_output_vtk(*refined_ij.i_pos_j_pos,"ipjp.vts");

  Tree_Node tree(grid);
  auto tree_ptr = std::make_shared<Tree_Node>(grid);
  // tree_ptr->Give_Variables_Test();
  // tree_ptr->Auto_Refine_Test();
  // tree_ptr->Auto_Refine_Test();
  // tree_ptr->Auto_Refine_Test();
  // tree_ptr->Auto_Refine_Test();
  // tree_ptr->Auto_Refine_Test();
  //tree_ptr->Auto_Refine_Test();
  //tree_ptr->Auto_Refine_Test();
  //tree_ptr->Auto_Refine_Test();
  //tree_ptr->Auto_Refine_Test();


  error_norms refined = tree_ptr->Compute_Errors();
  std::cout<<"After refined:"<<std::endl;
  std::cout<<"L1    = "<<refined.L1<<std::endl;
  std::cout<<"L2    = "<<refined.L2<<std::endl;
  std::cout<<"L_inf = "<<refined.L_inf<<std::endl;
  std::cout<<"number_of_cells="<<refined.num_cells<<std::endl;
  int num_blocks = tree_ptr->Write_All_VTK(0);

  //write .vtm file
  ofstream new_grid("whole_grid.vtm");
  new_grid<<"<?xml version="<<"\"1.0\""<<"?>\n";
  new_grid<<"<VTKFile type="<<"\"vtkMultiBlockDataSet\""<<" version="<<"\"0.1\""<<" byte_order="<<"\"LittleEndian\"\n";
  new_grid<<"compressor="<<"\"vtkZLibDataCompressor\""<<">\n";
  new_grid<<"<vtkMultiBlockDataSet>\n";
  for(int i=0; i<num_blocks; ++i) {
    new_grid<<"<DataSet group="<<"\""<<i<<"\""<<" dataset="<<"\""<<0<<"\""<<" file="<<"\"grid"<<i<<".vts\"/>\n";
  }
  new_grid<<"</vtkMultiBlockDataSet>\n";
  new_grid<<"</VTKFile>";
  




//  Solution_Block_2D<Euler2D_cState> solution_block(grid);
//
//  Set_Initial_Condition_Shockbox(solution_block);
//
//  write_solution_vtk(solution_block,"initial.vts");
//
//  using ode_type = Second_Order_FV_Cartesian_Block_2D<Euler2D_cState, HLLE, Constant_Extrapolation>;
//
//  ode_type fv_discretization(solution_block);
//
//  Predictor_Corrector<ode_type> time_marcher;
//
//  time_marcher.time_march(fv_discretization,0.00075,0.5);
//
//  write_solution_vtk(fv_discretization.solution_block,"final.vts");


  return 0;
}
