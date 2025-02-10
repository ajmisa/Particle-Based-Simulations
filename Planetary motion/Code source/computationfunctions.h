#ifndef COMPUTATIONFUNCTIONS_H_
#define COMPUTATIONFUNCTIONS_H_

/**
 * @brief Calculation of the deviation of positions in comparison with the provided data in 2022. 
 * 
 * @param bodies Array for saving the characteristics of a celestial body.
 * @param reference_bodies  Array of celestial bodies' positions and velocities in 2022. 
 * @param num_bodies  Number of celestial bodies
 * @param deviation_filename The name of the ouput filename to store the deviation of the position in comparison with the reference in 2022.
 */
void calculate_positional_deviation(struct CelestialBody bodies[], struct CelestialBody reference_bodies[], int num_bodies, const char *deviation_filename);

/**
 * @brief Calculation of forces/accelarations. In addition, it calculates the total kinetic energy and potential energy.
 * 
 * @param bodies Array for saving the characteristics of a celestial bodies.
 * @param num_bodies  Number of celestial bodies
 * @param total_kinetic_energy The total kinetic energy 
 * @param total_potential_energy The total potential energy
 */
double compute_forces_and_potential(struct CelestialBody bodies[], int num_bodies, double *total_kinetic_energy, double *total_potential_energy);

#endif