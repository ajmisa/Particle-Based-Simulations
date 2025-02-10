#include <stdio.h>
#include <stdlib.h>
#include <math.h>  
#include "constants.h"
#include "structs.h"
#include "computationfunctions.h"

// Function to calculate positional deviation and write to a file
void calculate_positional_deviation(struct CelestialBody bodies[], struct CelestialBody reference_bodies[], int num_bodies, const char *deviation_filename) {
    // Open the file for writing
    FILE *deviation_file = fopen(deviation_filename, "w");
    if (deviation_file == NULL) {
        perror("Error opening deviation file");
        exit(EXIT_FAILURE);
    }

    // Write a header for the file
    fprintf(deviation_file, "%15s %20s %20s %20s\n", "Name", "X Deviation (%)", "Y Deviation (%)", "Z Deviation (%)");

    // Loop through each body and calculate the positional deviation
    for (int i = 0; i < num_bodies; i++) {
        double dx = (bodies[i].r.x - reference_bodies[i].r.x) / reference_bodies[i].r.x * 100;
        double dy = (bodies[i].r.y - reference_bodies[i].r.y) / reference_bodies[i].r.y * 100;
        double dz = (bodies[i].r.z - reference_bodies[i].r.z) / reference_bodies[i].r.z * 100;

        // Write the deviations to the file
        fprintf(deviation_file, "%15s %20.10f %20.10f %20.10f\n", bodies[i].name, dx, dy, dz);
    }

    // Close the file
    fclose(deviation_file);
}

// Function to compute forces and total potential energy
double compute_forces_and_potential(struct CelestialBody bodies[], int num_bodies, double *total_kinetic_energy, double *total_potential_energy) {
    // Reset the energies
    *total_kinetic_energy = 0.0;
    *total_potential_energy = 0.0;

    // Reset accelerations to zero before calculating new forces
    for (int i = 0; i < num_bodies; i++) {
        bodies[i].a.x = 0.0;
        bodies[i].a.y = 0.0;
        bodies[i].a.z = 0.0;
    }

    // Loop through all unique pairs of bodies to calculate forces and potential energy
    for (int i = 0; i < num_bodies; i++) {
        // Kinetic energy calculation
        double kinetic_energy_i = 0.5 * bodies[i].mass * 
                                  (bodies[i].v.x * bodies[i].v.x +
                                   bodies[i].v.y * bodies[i].v.y +
                                   bodies[i].v.z * bodies[i].v.z);
        // Calculation of total Kinetic Energy
        *total_kinetic_energy += kinetic_energy_i;

        for (int j = i + 1; j < num_bodies; j++) {
            if (bodies[i].mass == 1) {
                continue;  // Massless bodies don't contribute forces
            }

            // Calculation of position difference between body i and j
            struct Vec3D r_ij;
            r_ij.x = bodies[j].r.x - bodies[i].r.x;
            r_ij.y = bodies[j].r.y - bodies[i].r.y;
            r_ij.z = bodies[j].r.z - bodies[i].r.z;

            // Calculation of the distance between body i and j
            double distance = sqrt(r_ij.x * r_ij.x + r_ij.y * r_ij.y + r_ij.z * r_ij.z);

            if (distance < 1e-10) {
                continue; // If the distance is too small skip the calculation
            }

            // Calculate the force/acceleration of the celestial bodies
            double force_magnitude = G * bodies[i].mass * bodies[j].mass / (distance * distance * distance);
            bodies[i].a.x += (force_magnitude * r_ij.x) / bodies[i].mass;
            bodies[i].a.y += (force_magnitude * r_ij.y) / bodies[i].mass;
            bodies[i].a.z += (force_magnitude * r_ij.z) / bodies[i].mass;

            // The massive bodies contribute forces
            if (bodies[j].mass != 1) {  
                bodies[j].a.x -= (force_magnitude * r_ij.x) / bodies[j].mass;
                bodies[j].a.y -= (force_magnitude * r_ij.y) / bodies[j].mass;
                bodies[j].a.z -= (force_magnitude * r_ij.z) / bodies[j].mass;
            }

            // Calculation of potential energy bewteen the pairs
            double potential_energy_ij = -G * bodies[i].mass * bodies[j].mass / distance;

            // Calculation of the total potential potential energy
            *total_potential_energy += potential_energy_ij;
        }

    }

    // Return the total Energy of the system. In addition, it stores the total kinetic energy and potential energy.
    return *total_kinetic_energy + *total_potential_energy;
}