#ifndef STRUCTURES_H_
#define STRUCTURES_H_

struct Vec3D
/* Structure to store x, y, z coordinates of a 3D vector */
{
    double x, y, z;
};

struct CelestialBody
/* Structure to store a celestial body characteristics: id, name, mass, velocity, position, acceleration and half velocity (velocty-Verlet algorithm)*/
{
    int id;                         /* id of celestial body */
    char name [30];                 /* name of celestial body */
    double mass;                    /* mass */
    struct Vec3D r;                 /* positions */
    struct Vec3D v;                 /* velocities */
    struct Vec3D a;                 /* acceleration */
    struct Vec3D half_velocity;     /* half velocity for velocity-Verlet computation */
};

#endif