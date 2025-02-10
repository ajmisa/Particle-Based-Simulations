# Planetary Motion

This repository contains the first assignment for the **Particle-Based Simulations** course in the **Master's program in Chemical Engineering** at **Eindhoven University of Technology (TU/e)**.

## 📌 Overview
This assignment focuses on simulating planetary motion using **Newton’s laws of motion** and **gravitational interactions**. The numerical integration methods used include:
- **Velocity-Verlet Algorithm** (primary method for integration)
- **Explicit Euler and Symplectic Euler methods** (for comparison)

## 📁 Contents
- **Report (PDF):** Contains detailed analysis, methodology, and results.
- **Source Code (C):** Implements planetary motion simulation.
- **Data Files:** Initial conditions for celestial bodies.
- **Results:** Energy conservation analysis and trajectory visualizations.

## 🛠 Implementation Details
- **Force Calculation:** Computes gravitational forces between unique body pairs.
- **Velocity-Verlet Integration:** Ensures stability and accuracy in long-term simulations.
- **Handling Zero-Mass Bodies:** Massless bodies follow gravitational forces but do not exert forces.
- **Energy Conservation Analysis:** Tracks kinetic and potential energy to verify numerical stability.

## 📊 Results and Analysis
- **Orbital Trajectories:** Comparison with NASA’s planetary data.
- **Energy Conservation:** Validates numerical stability.
- **Algorithm Comparison:** Evaluates trade-offs between different integration schemes.

📧 **Author:** Adam Misa and Group 18  
📅 **Date:** September 20, 2024

