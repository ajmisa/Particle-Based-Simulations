#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "constants.h"
#include "structs.h"

// Set the parameters of this simulation
void set_parameters(struct Parameters *p_parameters)

{
  // Methane Parameters (CH4) - Stored in type[0]
  p_parameters->type[0].mass = 16.04; // Mass of methane particle [amu]
  p_parameters->type[0].epsilon = 148 * kb;     // LJ epsilon for methane [eV]
  p_parameters->type[0].sigma = 3.73;                       // LJ sigma for methane [A]

  // Ethane Parameters (CH3 groups) - Stored in type[1]
  p_parameters->type[1].mass = 30.07;         // Mass of ethane molecule [amu]
  p_parameters->type[1].epsilon = 98 * kb;      // LJ epsilon for ethane [eV]
  p_parameters->type[1].sigma = 3.75;                       // LJ sigma for ethane [A]

  // Use Lorentz-Berthelot mixing rules to calculate mixed interactions
  p_parameters->sigma_mixed = (p_parameters->type[0].sigma + p_parameters->type[1].sigma) / 2.0;
  p_parameters->epsilon_mixed = sqrt(p_parameters->type[0].epsilon * p_parameters->type[1].epsilon);

  // Here i make sure that the box size is big enough. No reason for it, it just works
  p_parameters->L = (struct Vec3D){14.938*p_parameters->sigma_mixed, 14.938*p_parameters->sigma_mixed, 14.938*p_parameters->sigma_mixed};//box size


  // lectyre slides suggestted that the r_cut and r_shell could be multiplied with sigma_mixed to get reasonable values.
  // I did so and good values came from the simulation so i think for now it is good!
  p_parameters->r_cut = 2.5*p_parameters->sigma_mixed;            //cut-off distance used for neigbor list
  p_parameters->r_shell = 0.4*p_parameters->sigma_mixed;          //shell thickness for neighbor list


  p_parameters->num_dt_pdb = 200;                           //number of time steps in between pdb outputs
  strcpy(p_parameters->filename_pdb, "trajectories");       //filename (without extension) for pdb file
  p_parameters->rescale_output = 1;                         //factor used to rescale output lengthscale (Most visualisation programs identify bonds based on distances of order 1)
  p_parameters->load_restart = 0;                           //if equal 1 restart file is loaded
  strcpy(p_parameters->restart_in_filename, "restart.dat"); //filename for loaded restart file
  p_parameters->num_dt_restart = 1000;                      // number of time steps between saves
  strcpy(p_parameters->restart_out_filename, "restart.dat");//filename for saved restart file


  if (p_parameters->r_cut > p_parameters->L.x / 2.0)
    fprintf(stderr, "Warning! r_cut > Lx/2");
  if (p_parameters->r_cut > p_parameters->L.y / 2.0)
    fprintf(stderr, "Warning! r_cut > Ly/2");
  if (p_parameters->r_cut > p_parameters->L.z / 2.0)
    fprintf(stderr, "Warning! r_cut > Lz/2");
}
