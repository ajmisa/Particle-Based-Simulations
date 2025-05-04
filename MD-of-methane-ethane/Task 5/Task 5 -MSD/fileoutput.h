#ifndef FILEOUTPUT_H_
#define FILEOUTPUT_H_

/**
 * @brief Output particle positions to a pdb file for visualization in tools such as OVITO.
 * @param reset 1: write new file or overwrite existing file, 0: append data
 * @param p_parameters used members: filename_pdb
 * @param p_vectors used members: r
 * @param time time stamp
 */
void record_trajectories_pdb(int reset, struct Parameters * p_parameters, struct Vectors * p_vectors, double time);

/**
 * @brief Output particle positions to xyz file
 * @param reset 1: write new file or overwrite existing file, 0: append data
 * @param p_parameters used members: filename_xyz
 * @param p_vectors used members: r
 * @param time time stamp
 */
void record_trajectories_xyz(int reset, struct Parameters * p_parameters, struct Vectors * p_vectors, double time);

/**
 * @brief Save a restart file
 * 
 * @param p_parameters 
 * @param p_vectors 
 */
void save_restart(struct Parameters *p_parameters, struct Vectors *p_vectors);

/**
 * @brief Load a restart file
 * 
 * @param p_parameters 
 * @param p_vectors 
 */
void load_restart(struct Parameters *p_parameters, struct Vectors *p_vectors);

/**
 * @brief Calculation of the deviation of positions in comparison with the provided data in 2022. 
 * 
 * @param filename  File's name.
 * @param time_step Time step of the simulation in seconds. 
 * @param kinetic_energy The total kinetic energy of the system. 
 * @param potential_energy The total potential energy of the system. 
 * @param total_energy The total energy of the system (K+U) 
 */
void write_energy_to_file(const char *filename, int time_step, double kinetic_energy, double potential_energy, double total_energy, double temp);


/**
 * @brief removes old files so they dont cause problems when appending new information
 */
void removeFile(const char *fileName);

/**
 * @brief Write the msd of ethane and methane into csv file
 */
void write_msd_data(const char *msd_output, unsigned long step, double msd[], double diffusion_coefficient[]);

#endif /* FILEOUTPUT_H_ */