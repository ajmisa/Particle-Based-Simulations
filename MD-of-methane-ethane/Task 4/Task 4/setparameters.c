#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "constants.h"
#include "structs.h"

double mass[2];

// Set the parameters of this simulation
void set_parameters(struct Parameters *p_parameters)
{
  // System parameters
  p_parameters->Temp = 298;                                // Temperature in K
  p_parameters->kT = kb * p_parameters->Temp;              // Thermal energy (kB * T)

  // Set particle type-specific properties for Methane (index 0) and Ethane (index 1)
  p_parameters->particle_types[0].mass = 16.04;            // Methane mass (amu)
  p_parameters->particle_types[0].sigma = 3.73;            // Methane Lennard-Jones sigma (A)
  p_parameters->particle_types[0].epsilon = 148 * kb;      // Methane Lennard-Jones epsilon (eV)

  p_parameters->particle_types[1].mass = 30.07;            // Ethane mass (amu)
  p_parameters->particle_types[1].sigma = 3.75;            // Ethane Lennard-Jones sigma (A)
  p_parameters->particle_types[1].epsilon = 98 * kb;       // Ethane Lennard-Jones epsilon (eV)

  // The parameters below control core functionalities of the code, but many values will need to be changed
  p_parameters->num_part = 2000;                            //number of particles
  p_parameters->num_dt_steps = 2000;                        //number of time steps
  p_parameters->exclude_12_nb = 1;                          // 1-2 connected atoms exluded from non-bonded interactions 
  p_parameters->exclude_13_nb = 1;                          // 1-3 connected atoms exluded from non-bonded interactions    
  p_parameters->dt = 0.01;                                  //integration time step
  
  // Here i make sure that the box size is big enough. No reason for it, it just works
  p_parameters->L = (struct Vec3D){14.938 * 3.75, 14.938 * 3.75, 14.938 * 3.75};//box size


  // lectyre slides suggestted that the r_cut and r_shell could be multiplied with sigma to get reasonable values.
  // I did so and good values came from the simulation so i think for now it is good!
  p_parameters->r_cut = 2.5 * 3.75;             //cut-off distance used for neigbor list
  p_parameters->r_shell = 0.4 * 3.75;          //shell thickness for neighbor list


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
