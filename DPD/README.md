# Dissipative Particle Dynamics (DPD) Simulation

This project implements the Dissipative Particle Dynamics (DPD) method to study the mesoscopic behavior of binary mixtures and polymeric liquids. The method is based on the model introduced by Groot and Warren (J. Chem. Phys., 1997) and is developed by extending a Molecular Dynamics (MD) code framework.

## Overview

The simulation investigates equilibrium structure, dynamics, and phase separation in both monomeric and polymeric systems using DPD. The core features include soft conservative repulsion, dissipative friction, and stochastic random forces. These interactions are used to model realistic fluid behavior at a mesoscopic scale.

## Objectives

- Implement conservative, dissipative, and random forces as per the DPD method
- Verify correct implementation through:
  - Energy conservation
  - Ideal gas structure
  - Velocity distributions
- Validate the model through:
  - Radial distribution function (RDF)
  - Binary mixture phase separation
  - Flory–Huggins χ-parameter
  - Chain simulations for N = 2, 4, 8

## Implementation Highlights

- **Conservative Force**: Soft repulsion scaled by interaction parameter \( a_{ij} \)
- **Dissipative Force**: Velocity-dependent damping proportional to relative velocity
- **Random Force**: Thermostatting via noise balanced by the fluctuation-dissipation theorem
- **Bonded Interactions**: Harmonic springs for modeling chains
- **Integration**: Velocity-Verlet scheme adapted for DPD

## Verification

- **Energy Conservation**: Confirmed for monomers with only conservative forces
- **Ideal Gas Structure**: Observed RDF \( g(r) \approx 1 \) using only dissipative and random forces
- **Velocity Statistics**: Distribution matched the Maxwell–Boltzmann profile

## Validation

- **Radial Distribution Function**: Matches Groot & Warren’s results at ρ = 3
- **Phase Separation**: Binary mixtures with varying \( a_{AB} \) showed domain formation
- **χ-Parameter**:
  - Calculated using volume fractions
  - Results agree with theoretical relation \( \chi = 0.286 \cdot \Delta a \)
- **Chains**: Bond connectivity and phase behavior validated for N = 2, 4, and 8

## File Structure

- `src/` – Source code files (`.c`, `.h`) implementing the DPD simulation
- `output/` – Velocity histograms, RDF data, and χ-parameter values
- `PBS_assignment_3.pdf` – Final project report
- `A3-DPD.pdf` – Assignment handout
- `GROOT-RD-1997-J-CHEM-PHYS-V107-P11.pdf` – Reference paper

## References

Groot, R. D., & Warren, P. B. (1997). Dissipative particle dynamics: Bridging the gap between atomistic and mesoscopic simulation. *Journal of Chemical Physics*, 107(11), 4423–4435. https://doi.org/10.1063/1.474784

## Author

Adam Jordani Misa  
MSc Chemical Engineering – TU/e  
Email: aj.misa@outlook.com  
LinkedIn: https://www.linkedin.com/in/adam-misa
