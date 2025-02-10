#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h> 
#include "constants.h"
#include "structs.h"
#include "filereadwrite.h"

// Function to read data from file and set mass to 1 if it is zero
void read_data_file(const char *filename, struct CelestialBody bodies[]) {
    
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    // Skip the first line (titles)
    char header[200];
    fgets(header, sizeof(header), file);

    // Read the data file line by line
    for (int i = 0; i < NUM_BODIES; i++) {
        fscanf(file, "%d %s %lf %lf %lf %lf %lf %lf %lf",
               &bodies[i].id, bodies[i].name, &bodies[i].mass,
               &bodies[i].r.x, &bodies[i].r.y, &bodies[i].r.z,
               &bodies[i].v.x, &bodies[i].v.y, &bodies[i].v.z);

        // If mass is zero, set it to 1
        if (bodies[i].mass == 0.0) {
            bodies[i].mass = 1;
        }
    }

    fclose(file);
}

// Output positions of bodies to a file every few time steps 
void output_positions(const char *filename, struct CelestialBody bodies[], int num_bodies, int time_step, int num_steps, int frequency, bool write_all_bodies) {
    
    if (time_step % frequency == 0 || time_step == num_steps - 1) {
        FILE *file = fopen(filename, "a");
        if (file == NULL) {
            perror("Error opening file");
            exit(EXIT_FAILURE);
        }

        // Print a header for each time step
        fprintf(file, "Time Step: %d\n", time_step);

        // Print a header once or as needed
        if (time_step == 0) { 
            fprintf(file, "%s %12s %20s %20s %20s %20s %20s %20s %20s\n",
                    "#     id", "           name", "           mass (kg)", "               X (m)", "               Y (m)", "               Z (m)", "            VX (m/s)", "            VY (m/s)", "            VZ (m/s)");
        }

        for (int i = 0; i < num_bodies; i++) {
            
        // If user wants all bodies to be written, skip name checks
        if (write_all_bodies) {
            fprintf(file, "%8d %15s %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n",
                    bodies[i].id, bodies[i].name, bodies[i].mass,
                    bodies[i].r.x, bodies[i].r.y, bodies[i].r.z,
                    bodies[i].v.x, bodies[i].v.y, bodies[i].v.z);
        } 
        // Otherwise, only write specific bodies (Sun, Jupiter, Saturn, Uranus, Neptune and Pluto)
        else if (strcmp(bodies[i].name, "Sun") == 0 ||
                 strcmp(bodies[i].name, "Jupiter") == 0 ||
                 strcmp(bodies[i].name, "Saturn") == 0 ||
                 strcmp(bodies[i].name, "Uranus") == 0 ||
                 strcmp(bodies[i].name, "Neptune") == 0 ||
                 strcmp(bodies[i].name, "Pluto") == 0) {
            
            fprintf(file, "%8d %15s %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e %20.12e\n",
                    bodies[i].id, bodies[i].name, bodies[i].mass,
                    bodies[i].r.x, bodies[i].r.y, bodies[i].r.z,
                    bodies[i].v.x, bodies[i].v.y, bodies[i].v.z);
        }
        }
        fprintf(file, "\n");  // Separate data blocks with a newline

        fclose(file);
    }
}

// Function to write energy data to a file
void write_energy_to_file(const char *filename, int time_step, double kinetic_energy, double potential_energy, double total_energy) {
    FILE *file = fopen(filename, "a");
    if (file == NULL) {
        perror("Error opening energy file");
        exit(EXIT_FAILURE);
    }
    
    if (time_step ==0) {
    fprintf(file, "%20s %20s %20s %20s\n",
            "       Time Step", "    Kinetic Energy J", "    Potential Energy J", "      Total Energy J");
    }

    fprintf(file, "%20i %20.12e %20.12e %20.12e\n",
                    time_step, kinetic_energy, potential_energy, total_energy);

    fclose(file);
}