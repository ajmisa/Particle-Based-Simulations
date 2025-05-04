# Molecular Dynamics of a Supercritical Methane–Ethane Mixture

This project simulates the structural and dynamic properties of a binary methane–ethane mixture using Molecular Dynamics (MD) techniques. It was completed as part of the Particle-Based Simulations (6EMA02) course at TU/e (Eindhoven University of Technology).

## Overview

The aim is to investigate a 50:50 mole fraction mixture of methane and ethane in the supercritical phase at 298 K and 200 kg/m³. The simulation builds upon a progressively developed C codebase and implements key MD features such as:

- Velocity-Verlet integration
- Thermostatting (Berendsen)
- Bonded and non-bonded interactions
- Radial distribution functions (RDF)
- Mean squared displacement (MSD) and diffusion analysis
- Velocity distribution and comparison to Maxwell–Boltzmann statistics

## Project Structure

- `src/` – Main C codebase with `.c` and `.h` files
- `data/` – Output files such as energy, RDF, MSD, and velocity logs
- `plots/` – MATLAB scripts for postprocessing and visualization
- `report.pdf` – Full technical report documenting methods, analysis, and results

## Key Features

### Particle Models
- **Methane (CH₄):** Modeled as a single united atom (Lennard-Jones particle)
- **Ethane (C₂H₆):** Modeled as two bonded CH₃ groups using Hookean springs

### Interactions
- **Non-Bonded:** Lennard-Jones potential with Lorentz-Berthelot mixing rules
- **Bonded:** Harmonic bond potential with equilibrium distance \( r_0 = 1.54 \) Å

### Thermostat
- **Berendsen Thermostat** used for NVT ensemble simulations

### Analysis
- **Energy conservation** in NVE/NVT ensembles
- **Temperature control** validation via energy and temperature evolution
- **Velocity distributions** compared to Maxwell-Boltzmann theory
- **RDF** to evaluate short-range order between particle types
- **MSD and diffusion coefficients** extracted from long-time dynamics

## Usage

1. Compile the code using your preferred C compiler (e.g., `gcc *.c -lm -o md_sim`)
2. Run the simulation with appropriate configuration parameters
3. Use MATLAB scripts in `plots/` to generate velocity histograms, MSD curves, and RDF plots

## Results

- Demonstrated stable temperature regulation via the thermostat
- Calculated diffusion coefficients consistent with literature
- Identified RDF patterns and velocity distributions characteristic of supercritical mixtures

## Author

Adam Misa  
TU/e – MSc in Chemical Engineering  
Email: aj.misa@outlook.com  
LinkedIn: https://www.linkedin.com/in/adam-misa

