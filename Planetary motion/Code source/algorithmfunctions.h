#ifndef ALGORITHMFUNCTIONS_H_
#define ALGORITHMFUNCTIONS_H_

/**
 * @brief Velocity-Verlet algorithm
 * 
 * @param bodies Array for saving the characteristics of a celestial bodies.
 * @param num_bodies  Number of celestial bodies
 * @param dt_s Time step in seconds
 */
void velocity_verlet(struct CelestialBody bodies[], int num_bodies, double dt_s);

/**
 * @brief Symplectic Euler algorithm
 * 
 * @param bodies Array for saving the characteristics of a celestial bodies.
 * @param num_bodies  Number of celestial bodies
 * @param dt_s Time step in seconds
 */
void symplectic_euler(struct CelestialBody bodies[], int num_bodies, double dt_s);

/**
 * @brief Explicit Euler algorithm
 * 
 * @param bodies Array for saving the characteristics of a celestial bodies.
 * @param num_bodies  Number of celestial bodies
 * @param dt_s Time step in second
 */
void explicit_euler(struct CelestialBody bodies[], int num_bodies, double dt_s);

#endif