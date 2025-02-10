#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include "structs.h"
#include "computationfunctions.h"
#include "algorithmfunctions.h"

// Function to update positions and velocities using velocity-Verlet algorithm directly
void velocity_verlet(struct CelestialBody bodies[], int num_bodies, double dt_s) {

    // First, compute half-step velocities and full-step positions
    for (int i = 0; i < num_bodies; i++) {

        // Calculate v_i(t + dt_s/2)
        bodies[i].half_velocity.x = bodies[i].v.x + (bodies[i].a.x * dt_s / 2.0);
        bodies[i].half_velocity.y = bodies[i].v.y + (bodies[i].a.y * dt_s / 2.0);
        bodies[i].half_velocity.z = bodies[i].v.z + (bodies[i].a.z * dt_s / 2.0);

        // Update positions using half-step velocities
        bodies[i].r.x += bodies[i].half_velocity.x * dt_s;
        bodies[i].r.y += bodies[i].half_velocity.y * dt_s;
        bodies[i].r.z += bodies[i].half_velocity.z * dt_s;
    }

    // Now recompute forces to get F(t + Δt)
    double total_kinetic_energy, total_potential_energy;
    compute_forces_and_potential(bodies, num_bodies, &total_kinetic_energy, &total_potential_energy);  // This updates the acceleration (force) at t + Δt

    // Finally, compute full-step velocities
    for (int i = 0; i < num_bodies; i++) {
        bodies[i].v.x = bodies[i].half_velocity.x + (bodies[i].a.x * dt_s / 2.0);
        bodies[i].v.y = bodies[i].half_velocity.y + (bodies[i].a.y * dt_s / 2.0);
        bodies[i].v.z = bodies[i].half_velocity.z + (bodies[i].a.z * dt_s / 2.0);
    }
}

void symplectic_euler(struct CelestialBody bodies[], int num_bodies, double dt_s) {
    
    double total_kinetic_energy, total_potential_energy;
    compute_forces_and_potential(bodies, num_bodies, &total_kinetic_energy, &total_potential_energy);
   
    // Update velocities
    for (int i = 0; i < num_bodies; i++) {
        bodies[i].v.x += bodies[i].a.x * dt_s;
        bodies[i].v.y += bodies[i].a.y * dt_s;
        bodies[i].v.z += bodies[i].a.z * dt_s;
        // Update positions
        bodies[i].r.x += bodies[i].v.x * dt_s;
        bodies[i].r.y += bodies[i].v.y * dt_s;
        bodies[i].r.z += bodies[i].v.z * dt_s;
    }
}

void explicit_euler(struct CelestialBody bodies[], int num_bodies, double dt_s) {
    // Recompute forces
    double total_kinetic_energy, total_potential_energy;
    compute_forces_and_potential(bodies, num_bodies, &total_kinetic_energy, &total_potential_energy);
    
    // Update positions
    for (int i = 0; i < num_bodies; i++) {
        bodies[i].r.x += bodies[i].v.x * dt_s;
        bodies[i].r.y += bodies[i].v.y * dt_s;
        bodies[i].r.z += bodies[i].v.z * dt_s;
    }

    // Update velocities
    for (int i = 0; i < num_bodies; i++) {
        bodies[i].v.x += bodies[i].a.x * dt_s;
        bodies[i].v.y += bodies[i].a.y * dt_s;
        bodies[i].v.z += bodies[i].a.z * dt_s;
    }
}