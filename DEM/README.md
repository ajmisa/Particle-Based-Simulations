# DEM Simulation of Granular Flow in a Rotating Drum

This project models the flow of granular materials in a rotating drum using the Discrete Element Method (DEM). The work involves verification of particle collisions, implementation of rolling friction, and simulation of flow regimes observed in experimental studies. The simulations and analysis are conducted in the context of Assignment 4 for the Particle-Based Simulation course.

## Overview

Granular materials like spherical glass beads are simulated in a rotating drum setup. The primary goals include:
- Verifying particle-particle and particle-wall restitution coefficients
- Investigating rolling behavior with and without rolling friction
- Simulating flow regimes at different drum rotation speeds
- Quantifying the angle of repose and its dependence on particle properties and drum conditions

## Key Features

- DEM implementation based on soft-sphere interaction models
- Tangential and rolling friction support using torque models
- Cylindrical rotating drum with periodic boundary conditions
- Visualization-ready output (.xyz format) compatible with Ovito

## Tasks Performed

### 1. Collision Verification
- Modified parameters to match experimental restitution and friction values
- Simulated and measured restitution coefficients for particle-wall and particle-particle collisions

### 2. Rolling Verification
- Simulated single-particle rolling down a sloped surface
- Implemented rolling friction torque
- Compared observed accelerations with theoretical predictions

### 3. Rotating Drum Simulation
- Constructed a drum of 100 mm diameter and 16 mm depth
- Filled the drum to 35% and varied rotation speeds
- Observed transition between slumping, rolling, cascading, and cataracting regimes
- Measured the dynamic angle of repose at various Froude numbers
- Validated simulation results using experimental trends from Yang et al. (2008) and Pourandi et al. (2024)

## File Structure

- `src/` – Source code implementing DEM with wall functions
- `output/` – Simulation output files in `.xyz` format
- `A4_DEM.pdf` – Final report
- `PBS_4_group_10.pdf` – Assignment brief
- `Yang et al 2008.pdf`, `Pourandi et al 2024.pdf` – Reference articles

## References

- Yang, R.Y., et al. *Numerical simulation of particle dynamics in different flow regimes in a rotating drum*, Powder Technology, 188 (2008), 170–177.
- Pourandi, S., et al. *A mathematical model for the dynamic angle of repose of a granular material in the rotating drum*, Powder Technology, 446 (2024), 120176.

## Author

Adam Jordani Misa  
MSc Chemical Engineering – TU/e  
Email: aj.misa@outlook.com  
LinkedIn: https://www.linkedin.com/in/adam-misa
