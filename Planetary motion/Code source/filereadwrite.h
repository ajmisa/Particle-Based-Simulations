#ifndef FILESREADWRITE_H_
#define FILESREADWRITE_H_

/**
 * @brief Read the file (in the given format).
 * 
 * @param filename  File's name.
 * @param bodies Array for saving the characteristics of a celestial body.
 */
void read_data_file(const char *filename, struct CelestialBody bodies[]);

/**
 * @brief Outputs the positions and velocities o fthe celestial bodies. 
 * 
 * @param filename  File's name.
 * @param bodies Array for saving the characteristics of a celestial body.
 * @param num_bodies  Number of celestial bodies
 * @param time_step Time step of the simulation in seconds. 
 * @param frequency Frequency of the time steps that we want write the positions and velocties of the celestial bodies. 
 * @param write_all_bodies Boolean value to check if we want to write the positions of all celestial bodies.
 */
void output_positions(const char *filename, struct CelestialBody bodies[], int num_bodies, int time_step, int num_steps, int frequency, bool write_all_bodies);

/**
 * @brief Calculation of the deviation of positions in comparison with the provided data in 2022. 
 * 
 * @param filename  File's name.
 * @param time_step Time step of the simulation in seconds. 
 * @param kinetic_energy The total kinetic energy of the system. 
 * @param potential_energy The total potential energy of the system. 
 * @param total_energy The total energy of the system (K+U) 
 */
void write_energy_to_file(const char *filename, int time_step, double kinetic_energy, double potential_energy, double total_energy);

#endif