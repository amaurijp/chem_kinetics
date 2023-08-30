import os
import sys
root = r'C:\Users\amauri.paula\Dropbox\Minhas disciplinas\BCT Cinética Química\reaction_sim' #r'/home/amaurijp/Dropbox/Minhas disciplinas/BCT Cinética Química/reaction_sim'
sys.path.append(root + '/Modules')
from simulation import simulation


if __name__ == '__main__':
    #running the simulation
    print('running simulation')
    s = simulation(dt = 0.01,
                   n_particles = 10,
                   box_width = 200, 
                   box_height = 200,
                   potential = None,
                   root_path = root, 
                   save_path = root + '/Data')
    s.run(n_iterations = 100)