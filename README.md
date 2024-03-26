# Classical-Monte-Carlo-Particle-Sim
Classical Monte Carlo Particle Dynamics Simulation.
This repository contains Python code for a Monte Carlo simulation of particle dynamics in a 2D box. The simulation implements the Metropolis criterion to explore the configuration space of the particles and calculate their energy distribution at a given temperature.

REQUIREMENTS---------------------------
This program utilizes numpy:
    To install numpy tun the following command
    pip --install numpy


Features
    Lennard-Jones Potential: The potential energy between particles is modeled using the Lennard-Jones potential, which describes the interaction between neutral atoms or molecules.

    Periodic Boundary Conditions: The simulation applies periodic boundary conditions to ensure that particles do not interact with the box boundaries, effectively simulating an infinite system.

    Temperature Control: The temperature of the system can be adjusted to study its effect on particle behavior and energy distribution.

    Metropolis Criterion: The simulation uses the Metropolis criterion to determine whether a particle move is accepted based on the change in energy. Moves that lower the energy are always accepted, while moves that increase the energy have a probability of being accepted.

Output
  This program outputs final position data in a CSV file named position_data.csv to whichever directory the command prompt is currently in.
  To visualize the data open the CSV file in excel and plot on a dot plot. Each row containts a particles X then Y coordinate.

Usage
    To run the simulation call the Monte_Carlo(number of steps, steps remaining increment)
    You can also change the systems various parameters such as sigma, length and width of box (L), epsilon, and number of particles (N).
