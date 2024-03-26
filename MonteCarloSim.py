import numpy as np
import time
import math
import os

'''
This program simulates particle interactions based on given constants utilizing a monte-carlo satistical averaging method.
Utilizes periodic boundary conditions for both the physical location and the energy calculation
'''

'''                
Defining constants
'''
sigma = 3.4                                #Essentially size of argon atom in units of angstrom
epsilon = 0.238                            #Value of potential energy well depth in kilocalories per mole of argon gas
kb = 1.3806503 * (10 ** -23)               #Value of Boltzman constant
L = 150                                    #Size of box in angstroms
N = 50                                     #Number of particles in box
T = 0.1                                    #Temperature of box in kelvin
positions = L * np.random.rand(N, 2)       #Creates a 2 dimensional array of positions for N number of particles

'''
Measures potenial energy of one particle given a distance to another particle using lennard jones potenial energies
'''
def measure_PE(distance):
    if distance == 0:
        return 0
    else:
        return 4 * epsilon * ((sigma / distance) ** 12 - (sigma / distance) ** 6)

'''
Inputs:
    positions: Passes the positions array to the function
    chosen_particle_index: passes the chosen particles index to the function

Measures the total energy of the system based on a chosen particle
Also implements periodic boundary conditions for the energy.
''' 
def measure_total_energy(positions, chosen_particle_index):
    energy = 0.0
    for i in range(N):
        
        X_distance = positions[i][0] - positions[chosen_particle_index][0] 
        Y_distance = positions[i][1] - positions[chosen_particle_index][1]

#Periodic boundary conditions for the energy below
        
        if X_distance > L / 2: 
            X_distance -= L

        if Y_distance > L / 2:
            Y_distance -= L

        distance = math.sqrt((X_distance ** 2) + (Y_distance ** 2))
            
        energy += measure_PE(distance)
    return energy

'''
Inputs:
    energy_before: energy before a particle move was made
    energy_after: energy after a particle move was made

calculates the probability of a move and decides whether to keep it based on if the energy of the system is lowered or increased. 
AKA Metropolis criterion. A move that lowers the energy of the system is automatically accepted. While a move that increases the energy 
has a low probability of being accepted
'''
def metropolis(energy_before, energy_after):
    delta_energy = energy_after - energy_before 
    if delta_energy < 0 or np.random.rand() < np.exp(-delta_energy / (kb * T)):
        return True
    else:
        return False 

'''
Generates a random move within a range of -0.5 to 0.5 Angstroms of a particles current position. Move is in 2 dimensions.
'''
def generate_move():
    move = np.random.uniform(-.5, .5, size=2)
    return move

'''
Inputs
    positions_new: postitions of particles 

This function enforces physical periodic boundary conditions
It checks both the X and Y values of a specific particle and sees if it is outside of the dimensions of the box.
If the particle is outside the allowable area it subtracts the length or height of box from the particles X or Y position 
so that the particle appears on the exact opposite side of the box.
'''

def physical_periodic_boundary(new_position):

    for i in range(2):
        if new_position[i] > L:
            new_position[i] -= L

        if new_position[i] < 0:
            new_position[i] += L

    return new_position


"""
Inputs:
    steps: The number of moves which you would like the simulation to make

    specified_print_step: The increment with which you would like the system to print how many steps are left.
      For example specified_print_step = 1000 would print steps remaining every 1000 steps.

The main function which runs the simulation. Calculates the energy of the system before a move, makes a move, calculates energy of the system after the move
then checks the metropolis criteria to see if the move will be accepted or not. If metropolis criteria is less than 0, move is accepted.
Also prints moves remaing every specified number of moves moves, and time elapsed once simulation is finished.
"""
def Monte_Carlo(steps, specified_print_step):

    start = time.time()
    step_counter = steps
    chosen_particle_index = 0
    delta_energy_cumulative = 0

    for i in range(steps):                   #Increments over specified number of steps. Sequentially measures energy of particles starting from 0 to N -1
        if chosen_particle_index > N - 1:    #Checks to see if the particle index number is greater than the allowable index.      
            chosen_particle_index = 0        #If it is it sets the index back to 0 so that the functions starts iterating from index 0.    
                                                  
        positions_new = positions.copy()

        energy_before = measure_total_energy(positions, chosen_particle_index)      
                                                                                    
        move = generate_move()                                                                          


        positions_new[chosen_particle_index] += move

        positions_new[chosen_particle_index] = physical_periodic_boundary(positions_new[chosen_particle_index])     
                                                                                    
        energy_after = measure_total_energy(positions_new, chosen_particle_index)   
        delta_energy = energy_after - energy_before    
                                                                                     
        if metropolis(energy_before, energy_after):                                 
            positions[chosen_particle_index] = positions_new[chosen_particle_index]
            delta_energy_cumulative += delta_energy
        chosen_particle_index += 1

        step_counter -= 1
        
        if step_counter % specified_print_step == 0:
             print(step_counter, "Steps left")
        

    end = time.time()
    print(f'Time elapsed: {end - start}')
    print("Cumulative delta energy", delta_energy_cumulative)

'''
Copies position array into a csv file to be visualized using excel. Generally saves it to user directory   
'''
def write_csv():
    script_dir = os.path.dirname(__file__)
    file_path = os.path.join(script_dir, 'position_data.csv')

    np.savetxt(file_path, positions, delimiter=',', header='X,Y')
    


Monte_Carlo(10, 10000)
write_csv()

'''  
TO DO:
Export positions as GIF
    DONE iterate through particles not random, 
keep array of N individual particle energies. one for before and one for after move
    DONE calculate delta E. E old - E new. Use with metropolis criteria to decide if move is kept or not
    DONE Add to running total of energy of system
Start from crystalin state. Also try running at T = very low to get a crystaline structure to use.
Hexagonal close packing
Also starting from crystaline state where particles are closer than sigma
DONE Implement phyisical and energy periodic boundary conditions         
Don't take total energy initially, let run for a bit, equilibrate for a bit, then start running. During equilibrating don't store any data
Implement cutoff distance
'''           


     





