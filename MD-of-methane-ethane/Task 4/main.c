/******************************************************************************/ 
/*                                                                            */
/*  A Molecular Dynamics simulation of Lennard-Jones particles                */
/*                                                                            */
/*	This code is part of the course "Particle-based Simulations"              */
/*  taught at Eindhoven University of Technology.                             */
/*  No part of this code may be reproduced without permission of the author:  */
/*  Dr. Ir. E.A.J.F. Peters                                                   */
/*                                                                            */
/*  Dr. Ir. J.T. Padding:    version 1.1, 30/1/2013                           */
/*  Jeroen Hofman:           version 1.2, 28/7/2015                           */
/*  Dr. Ir. E.A.J.F. Peters: version 6.0, 17/9/2024    			              */
/******************************************************************************/ 

/** 
 * For the 2024 PBS assignment, the code needs to be extended:
 * 
 * - Implement a Berendsen thermostat in dynamics.c 
 * - Implement bonds in initialise_bonds in file initialise.c 
 * - Initialize vectors.type so particles get the proper type 
 * - Implement the bonded and non-bonded force in forces.c (Make forces type-dependent) 
 * - Change the particle position initialization to account for bond lengths and angles 
 * - Implement the needed on-the-fly data analysis
 */ 

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include "constants.h" 
#include "structs.h" 
#include "setparameters.h" 
#include "initialise.h" 
#include "nbrlist.h" 
#include "forces.h" 
#include "dynamics.h" 
#include "memory.h" 
#include "fileoutput.h" 

/** 
 * @brief Main MD simulation code. After initialization, 
 * a velocity-Verlet scheme is executed for a specified number of time steps. 
 * 
 * This simulation runs for a binary mixture of methane and ethane. 
 * 
 * Tasks implemented in this assignment:
 * 
 * 1. Implement the Berendsen thermostat to maintain the system temperature (NVT ensemble). 
 * 2. Modify the code to include bond-stretch interactions for ethane molecules. 
 * 3. Add support for multiple particle types (methane and ethane) and handle non-bonded forces accordingly.
 * 4. Perform on-the-fly analysis for temperature, energy, and molecular dynamics properties.
 * 
 * @return int 0 on success, non-zero on failure. 
 */ 
int main(void) 
{ 
    struct Vectors vectors; 
    struct Parameters parameters; 
    struct Nbrlist nbrlist; 
    size_t step; 
    double Ekin, Epot, time; 
    const char *writenergies_output = "energies.csv";

    // Step 1: Set the simulation parameters from input files
    set_parameters(&parameters); 

    // Step 2: Allocate memory for particles, forces, and neighbor lists
    alloc_memory(&parameters, &vectors, &nbrlist); 

    // Check if a restart is required
    if (parameters.load_restart == 1) 
    { 
        load_restart(&parameters, &vectors); 
        initialise_structure(&parameters, &vectors, &nbrlist); 
        step = 0; 
        time = 0.0; 
    } 
    else 
    {   
    // TODO: Initialize particle types (methane and ethane) in vectors.type array
    // TODO: Implement the bonds between ethane's CH3 groups in initialise_bonds (initialise.c)
        initialise(&parameters, &vectors, &nbrlist, &step, &time); 

    }

    // Step 3: Build the neighbor list for non-bonded interactions
    build_nbrlist(&parameters, &vectors, &nbrlist); 

    // Step 4: Calculate initial forces (non-bonded and bonded if implemented)
    Epot = calculate_forces(&parameters, &nbrlist, &vectors); 

    // Output initial particle positions in PDB format
    record_trajectories_pdb(1, &parameters, &vectors, time); 

    // Removes the old energies.csv file so that it does not cause problems
    removeFile("energies.csv");
    
    // Main MD loop using velocity-Verlet integration
    while (step < parameters.num_dt_steps) 
    { 
        step++; 
        time += parameters.dt; 
        
        // Update velocities (half-step)
        Ekin = update_velocities_half_dt(&parameters, &nbrlist, &vectors); 

        // TODO: Implement and apply the Berendsen thermostat to maintain temperature (dynamics.c)
        thermostat(&parameters, &vectors, Ekin); 

        // Update positions
        update_positions(&parameters, &nbrlist, &vectors); 

        // Apply boundary conditions
        boundary_conditions(&parameters, &vectors); 

        // Rebuild neighbor list if needed
        update_nbrlist(&parameters, &vectors, &nbrlist); 

        // Calculate forces for the current configuration (bonded forces if implemented)
        Epot = calculate_forces(&parameters, &nbrlist, &vectors); 

        // Final velocity update (half-step)
        Ekin = update_velocities_half_dt(&parameters, &nbrlist, &vectors); 

        int N = parameters.num_part;
        double N_free = 3 * (N-1); // degrees of freedom
    
        // Calculate temperature based on E = 1.5kBT formula
        double T_meas = (Ekin) / (N_free * kB);
        //double T_current = (2.0 * Ekin) / (3.0 * N);
        
        //put calculated energy values into a csv file
        write_energy_to_file(writenergies_output, (long unsigned)step, Ekin, Epot, Ekin + Epot, T_meas);

        // Output system state every 'num_dt_pdb' steps
        if (step % parameters.num_dt_pdb == 0) 
            record_trajectories_pdb(0, &parameters, &vectors, time); 

        // Save restart file every 'num_dt_restart' steps
        if (step % parameters.num_dt_restart == 0) 
            save_restart(&parameters, &vectors); 

        // TODO: Implement on-the-fly analysis of temperature, energy, and RDF

        // Print to the screen to monitor the progress of the simulation
        //printf("Step %lu, Time %f, Epot %f, Ekin %f, Etot %f\n", (long unsigned)step, time, Ekin, Epot, Epot + Ekin);
    } 

    // Save final state
    save_restart(&parameters, &vectors); 

    // Step 5: Free memory and clean up
    free_memory(&vectors, &nbrlist); 

    return 0; 
}
