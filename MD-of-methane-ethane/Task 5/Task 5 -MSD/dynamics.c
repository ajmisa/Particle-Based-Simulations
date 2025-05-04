#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "structs.h"

// This function updates particle positions using their velocities.
// The positions are advanced by one full time step (dt), and displacement vectors
// (dr) for one time step are updated. The displacement since the last neighbor 
// list creation (stored in p_nbrlist->dr) is also updated for each particle.
void update_positions(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors)
{
    struct Vec3D dr_loc;
    struct Vec3D *r = p_vectors->r;     // Particle positions
    struct Vec3D *dr = p_vectors->dr;   // Displacement in one timestep
    struct Vec3D *v = p_vectors->v;     // Particle velocities
    struct DeltaR *dr_nbrlist = p_nbrlist->dr;  // Displacement since last neighbor list creation
    size_t num_part = p_parameters->num_part;
    double dt = p_parameters->dt;

    // Loop over all particles to update their positions
    for (size_t i = 0; i < num_part; i++)
    {
        dr_loc.x = v[i].x * dt;  // Compute displacement in x-direction for one timestep
        dr_loc.y = v[i].y * dt;  // Compute displacement in y-direction for one timestep
        dr_loc.z = v[i].z * dt;  // Compute displacement in z-direction for one timestep

        dr[i] = dr_loc;          // Store the displacement for this timestep
        r[i].x += dr_loc.x;      // Update position in x-direction
        r[i].y += dr_loc.y;      // Update position in y-direction
        r[i].z += dr_loc.z;      // Update position in z-direction

        // Update the displacement since last neighbor list creation
        dr_nbrlist[i].x += dr_loc.x;
        dr_nbrlist[i].y += dr_loc.y;
        dr_nbrlist[i].z += dr_loc.z;
        dr_nbrlist[i].sq = (dr_nbrlist[i].x) * (dr_nbrlist[i].x) + 
                           (dr_nbrlist[i].y) * (dr_nbrlist[i].y) + 
                           (dr_nbrlist[i].z) * (dr_nbrlist[i].z);  // Square of total displacement since last neighbor list creation
    }
}

// This function updates particle velocities by half a time step using the current forces.
// The updated velocities are used in the velocity-Verlet integration scheme.
// The function also calculates and returns the kinetic energy of the system.
double update_velocities_half_dt(struct Parameters *p_parameters, struct Nbrlist *p_nbrlist, struct Vectors *p_vectors)
{
    double Ekin = 0.0;  // Initialize kinetic energy
    struct Vec3D *v = p_vectors->v;  // Particle velocities
    struct Vec3D *f = p_vectors->f;  // Forces acting on particles
    int *type = p_vectors->type;     // Particle types (0 = methane, 1 = ethane)

    // Loop over all particles and update their velocities
    for (size_t i = 0; i < p_parameters->num_part; i++)
    {   
        // Select the mass based on particle type
        double mass = p_parameters->particle_types[type[i]].mass;

        // Factor for velocity update, considering the particle's mass
        double factor = 0.5 / mass * p_parameters->dt;

        v[i].x += factor * f[i].x;  // Update velocity in x-direction
        v[i].y += factor * f[i].y;  // Update velocity in y-direction
        v[i].z += factor * f[i].z;  // Update velocity in z-direction
        Ekin += mass * (v[i].x * v[i].x + v[i].y * v[i].y + v[i].z * v[i].z);  // Accumulate kinetic energy, incl. the mass of the different type
    }

    // Final kinetic energy calculation
    Ekin = 0.5 * Ekin;
    return Ekin;  // Return the system's kinetic energy
}

// This function applies periodic boundary conditions to ensure particles stay inside the simulation box.
// If a particle moves beyond the box, it is wrapped around to the opposite side.
void boundary_conditions(struct Parameters *p_parameters, struct Vectors *p_vectors)
{
    struct Vec3D invL;  // Inverse of the box size
    struct Vec3D *r = p_vectors->r;  // Particle positions
    struct Vec3D L = p_parameters->L;  // Box dimensions
    size_t num_part = p_parameters->num_part;  // Number of particles

    invL.x = 1.0 / L.x;
    invL.y = 1.0 / L.y;
    invL.z = 1.0 / L.z;

    // Loop over all particles and apply periodic boundary conditions
    for (size_t i = 0; i < num_part; i++)
    {
        r[i].x -= L.x * floor(r[i].x * invL.x);  // Apply periodic boundary in x-direction
        r[i].y -= L.y * floor(r[i].y * invL.y);  // Apply periodic boundary in y-direction
        r[i].z -= L.z * floor(r[i].z * invL.z);  // Apply periodic boundary in z-direction
    }
}

// This function applies a thermostat to maintain the system's temperature.
// This function applies a thermostat to maintain the system's temperature.
void thermostat(struct Parameters *p_parameters, struct Vectors *p_vectors, double Ekin)
// TODO: Change velocities by thermostatting
{ 
    // Set desired temperature - set here as kT in the parameters 
    double T_target = 298;
    double N = p_parameters->num_part ; 
    double N_free = 3 * (N-1); // degrees of freedom
    struct Vec3D *v = p_vectors->v;  // Particle velocities
    
    // Calculate temperature based on E = 1.5kBT formula
    double T_meas = (Ekin) / (N_free * kB);

    double dt = p_parameters->dt; // time step 
    double tau = 0.1;  // tau-time constant for brendson thermostat
    
    // Compute the rescaling factor for velocities (lambda)
    double lambda = sqrt(1.0 + (dt / tau) * ((T_target / T_meas) - 1));
    
    // Apply the rescaling to the velocities
    for (size_t i = 0; i < p_parameters->num_part; i++) {

        v[i].x *= lambda;
        v[i].y *= lambda;
        v[i].z *= lambda;
    }

}

void apply_periodic_boundary_conditions(double *d, double box_length) {
    if (*d > 0.5 * box_length) {
        *d -= box_length;
    } else if (*d < -0.5 * box_length) {
        *d += box_length;
    }
}



// Function to calculate MSD and diffusion coefficient
void calculate_msd_and_diffusion(struct Parameters *p_parameters, struct Vectors *p_vectors, struct Vectors *initial_positions, double *msd, double *diffusion_coefficient) {
    double total_msd_methane = 0.0;
    double total_msd_ethane = 0.0;
    int count_methane = 0;
    int count_ethane = 0;
    
    // Loop over all particles to calculate the MSD
    for (size_t i = 0; i < p_parameters->num_part; i++) {
        double dx = p_vectors->r[i].x - initial_positions->r[i].x;
        double dy = p_vectors->r[i].y - initial_positions->r[i].y;
        double dz = p_vectors->r[i].z - initial_positions->r[i].z;

        // Apply periodic boundary conditions
        apply_periodic_boundary_conditions(&dx, p_parameters->L.x);
        apply_periodic_boundary_conditions(&dy, p_parameters->L.y);
        apply_periodic_boundary_conditions(&dz, p_parameters->L.z);

        double squared_displacement = dx * dx + dy * dy + dz * dz;
        
        // Sum MSD based on type
        if (p_vectors->type[i] == 0) {  // Methane
            total_msd_methane += squared_displacement;
            count_methane++;
        } else if (p_vectors->type[i] == 1) {  // Ethane
            total_msd_ethane += squared_displacement;
            count_ethane++;
        }

    }
    // Average MSD for each type
    double average_msd_methane = total_msd_methane / count_methane;
    double average_msd_ethane = total_msd_ethane / count_ethane;
    
    // Calculate the self-diffusion coefficient D = MSD / (6 * t)
    double time = p_parameters->dt * p_parameters->num_dt_steps;  // Total simulation time
    diffusion_coefficient[0] = average_msd_methane / (6 * time);  // Methane
    diffusion_coefficient[1] = average_msd_ethane / (6 * time);   // Ethane
    
    // Store MSD results
    msd[0] = average_msd_methane;
    msd[1] = average_msd_ethane;

}

