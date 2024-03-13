# Classical-Monte-Carlo-Particle-Sim
Classical Monte Carlo Particle Dynamics Simulation.

This repository contains Python code for a Monte Carlo simulation of particle dynamics in a 2D box. The simulation implements the Metropolis criterion to explore the configuration space of the particles and calculate their energy distribution at a given temperature.

Features
    Lennard-Jones Potential: The potential energy between particles is modeled using the Lennard-Jones potential, which describes the interaction between neutral atoms or molecules.

    Periodic Boundary Conditions: The simulation applies periodic boundary conditions to ensure that particles do not interact with the box boundaries, effectively simulating an infinite system.

    Temperature Control: The temperature of the system can be adjusted to study its effect on particle behavior and energy distribution.

    Metropolis Criterion: The simulation uses the Metropolis criterion to determine whether a particle move is accepted based on the change in energy. Moves that lower the energy are always accepted, while moves that increase the energy have a probability of being accepted.

Output
  This program outputs final position data in a csv file to whichever directory the command prompt is currently in.
