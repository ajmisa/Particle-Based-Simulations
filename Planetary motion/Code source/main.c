/**********************************************************************************/ 
/*                                                                                */
/*  Particle-Based Simulations: Assignemt 1 (2024)                                */
/*                                                                                */
/*  A simulation of the Solar system: Energy Conservation and Numerical Stability */
/*                                                                                */
/*	This code is part of the course "Particle-based Simulations"                  */
/*  No part of this code may be reproduced without permission of the authors:     */
/*  Adam Misa, Anita Jose, Luke de Jong                                           */
/*                                                                                */
/*  version 1.1, 19/09/2024                                                       */
/**********************************************************************************/ 

#include <stdio.h>          /* standard c input/output library */
#include <stdlib.h>         /* standard c library */
#include <math.h>           /* standard c math library */
#include <string.h>         /* String library */
#include <stdbool.h>        /* Boolean library*/
#include "constants.h" 
#include "structs.h" 
#include "computationfunctions.h" 
#include "algorithmfunctions.h" 
#include "filereadwrite.h" 

int main() {

    struct CelestialBody bodies[NUM_BODIES];  // Array of celestial bodies in 2021
    struct CelestialBody reference_bodies[NUM_BODIES]; // Array of celestial bodies in 2022

    // Read data from the file, initial positions and velocities
    read_data_file("solar_system_13sept2021.dat", bodies);

    // Read the reference data for comparison
    read_data_file("solar_system_13sept2022.dat", reference_bodies);

    double dt_d;  // Time step in days
    double dt_s;  // Time step in seconds
    int num_steps;  // Number of time steps
    int num_days;  // Number of year for simulation

    // Output file 
    const char *output_filename_velocity_verlet = "positions_over_time_velocity_verlet.txt";
    const char *energy_filename_velocity_verlet = "energy_output_velocity_verlet.txt";
    const char *deviation_filename_velocity_verlet = "positional_deviation_velocity_verlet.txt";  
    const char *output_filename_symplectic_euler = "positions_over_time_symplectic_euler.txt";
    const char *energy_filename_symplectic_euler = "energy_output_symplectic_euler.txt";
    const char *deviation_filename_symplectic_euler = "positional_deviation_symplectic_euler.txt"; 
    const char *output_filename_explicit_euler = "positions_over_time_explicit_euler.txt";
    const char *energy_filename_explicit_euler = "energy_output_explicit_euler.txt";
    const char *deviation_filename_explicit_euler = "positional_deviation_explicit_euler.txt"; 

    //variable for choosing the method
    int method;

    printf("Select the integration method:\n");
    printf("0: Velocity-Verlet\n");
    printf("1: Symplectic Euler\n");
    printf("2: Explicit Euler\n");
    printf("Enter choice (0-2): ");
    scanf("%d", &method);

    // Ensure the chosen method is within the valid range
    if (method < 0 || method > 2) {
        printf("Invalid method selection. Try again by choosing between 0-2.\n");
        return 1;
    }

    printf("Select the time step in days\n");
    printf("Note: you can also put a fraction of the day, let's say 0.5\n");
    scanf("%lf", &dt_d);
    
    dt_s = dt_d * 86400;

    printf("Select the number of steps. In this case the number of days\n");
    printf("that you want the simulation to run.\n");
    scanf("%i", &num_days);

    num_steps = num_days / dt_d;

    // Ask the user if they want to write all bodies or just the specific ones
    char user_input;
    bool write_all_bodies = false;  // Initialize as false (default to writing only specific bodies)

    printf("Do you want to write all celestial bodies to the file? (y/n): ");
    scanf(" %c", &user_input);

    if (user_input == 'y' || user_input == 'Y') {
        write_all_bodies = true;  // User wants to write all bodies
    }

    int output_frequency = 10;  // Output positions every 10 time steps

    // Compute the initial positions and velocities
    double initial_kinetic_energy, initial_potential_energy;
    double initial_energy = compute_forces_and_potential(bodies, NUM_BODIES, &initial_kinetic_energy, &initial_potential_energy);
    printf("Total Initial Energy: %e J\n", initial_energy); // Initial total Energy (J)



    // Simulation loop
    for (int step = 0; step < num_steps; step++) {
        switch (method){
            case 0:
                // Update positions and velocities using velocity-Verlet
                velocity_verlet(bodies, NUM_BODIES, dt_s);

                // Track energy change at each step and write to the file
                if (step % output_frequency == 0) {
                    double kinetic_energy, potential_energy;
                    double total_energy = compute_forces_and_potential(bodies, NUM_BODIES, &kinetic_energy, &potential_energy);
                    write_energy_to_file(energy_filename_velocity_verlet, step, kinetic_energy, potential_energy, total_energy);
                }
                // Output positions every 'output_frequency' steps
                output_positions(output_filename_velocity_verlet, bodies, NUM_BODIES, step, num_steps, output_frequency, write_all_bodies);
                break;
            case 1:
                // Update positions and velocities using velocity-Verlet
                symplectic_euler(bodies, NUM_BODIES, dt_s);

                // Track energy change at each step and write to the file
                if (step % output_frequency == 0) {
                    double kinetic_energy, potential_energy;
                    double total_energy = compute_forces_and_potential(bodies, NUM_BODIES, &kinetic_energy, &potential_energy);
                    write_energy_to_file(energy_filename_symplectic_euler, step, kinetic_energy, potential_energy, total_energy);
                }
                // Output positions every 'output_frequency' steps
                output_positions(output_filename_symplectic_euler, bodies, NUM_BODIES, step, num_steps, output_frequency, write_all_bodies);
                break; 
            case 2:
                // Update positions and velocities using velocity-Verlet
                explicit_euler(bodies, NUM_BODIES, dt_s);

                // Track energy change at each step and write to the file
                if (step % output_frequency == 0) {
                    double kinetic_energy, potential_energy;
                    double total_energy = compute_forces_and_potential(bodies, NUM_BODIES, &kinetic_energy, &potential_energy);
                    write_energy_to_file(energy_filename_explicit_euler, step, kinetic_energy, potential_energy, total_energy);
                }
                // Output positions every 'output_frequency' steps
                output_positions(output_filename_explicit_euler, bodies, NUM_BODIES, step, num_steps, output_frequency, write_all_bodies);
                break;             
        }
        // Check if we've reached the 365th day to compute positional deviation
        if (step == (int)(364 / dt_d)) { 
            printf("Calculating positional deviation at 365 days.\n");
            switch (method) {
                case 0:
                    calculate_positional_deviation(bodies, reference_bodies, NUM_BODIES, deviation_filename_velocity_verlet);
                    break;
                case 1:
                    calculate_positional_deviation(bodies, reference_bodies, NUM_BODIES, deviation_filename_symplectic_euler);
                    break;
                case 2:
                    calculate_positional_deviation(bodies, reference_bodies, NUM_BODIES, deviation_filename_explicit_euler);
                    break;
            }
        }
    }

    // Final energy after the simulation
    double final_kinetic_energy, final_potential_energy;
    double final_energy = compute_forces_and_potential(bodies, NUM_BODIES, &final_kinetic_energy, &final_potential_energy);
    printf("Total Final Energy: %e J\n", final_energy);

    // Energy change after a year of simulation
    double energy_change = (final_energy - initial_energy) / initial_energy * 100;
    printf("Energy Change after a year of simulation: %.10f%%\n", energy_change);

    printf("The simulation finished successfully\n");

    return 0;
}