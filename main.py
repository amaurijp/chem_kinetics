import os
import sys
root = r'/home/amaurijp/Dropbox/Minhas disciplinas/BCT Cinética Química/simulacao' #r'C:\Users\amauri.paula\Dropbox\Minhas disciplinas\BCT Cinética Química\simulacao'
sys.path.append(root + '/Modules')
from simulation import simulation


if __name__ == '__main__':
    #running the simulation
    print('running simulation')
    s = simulation(dt = 0.01, 
                   n_particles = 500, 
                   radius_input = [10, 10], 
                   max_velocity = 100, 
                   mass = 0.2,
                   temperature = 500,
                   box_width = 1500, 
                   box_height = 1500, 
                   n_dimension = 2, 
                   root_path = root, 
                   save_path = root + '/Data')
    s.run(n_iterations = 100)