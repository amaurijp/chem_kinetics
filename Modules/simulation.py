import os
import numpy as np
import matplotlib.pyplot as plt
import imageio
import math


from particle_props import particle_props
from functions import check_overlap, get_particles_distance, plot_particles, plot_v_histogram, correct_counter, get_force_between_particles
from initial_conds import initial_conds


class simulation():
    
    #instantiation
    def __init__(self, 
                 dt = 0.01, 
                 n_particles = 2, 
                 mass= 0.2, 
                 temperature = 500,
                 box_width = 200,
                 box_height = 200, 
                 n_dimension = 2, 
                 potential = None,
                 root_path = None, 
                 save_path = None):
        
        self.n_particles = n_particles
        self.n_dimension = n_dimension
        self.dt = dt
        self.mass = mass
        self.potential = potential
        self.temperature = temperature
        self.box_width = box_width
        self.box_height = box_height

        
        #diretorio
        self.save_path = save_path
        self.cwd = root_path


    #getting initial conditions
    def get_initial_conditions(self):
        
        #instance of p_props matrix
        p_props = particle_props(self.n_particles, self.n_dimension)

        #getting initial setup
        setup_simulation = initial_conds(p_props = p_props, 
                                        radius_input = [10, 20], 
                                        random_val_for_r = True,
                                        random_val_for_v = True,
                                        random_val_for_a = False,
                                        max_velocity = 150,
                                        max_accel = 10,
                                        n_particles = self.n_particles,
                                        n_dimension = self.n_dimension, 
                                        box_width = self.box_width, 
                                        box_height = self.box_height)
        
        self.p_props = setup_simulation.update_p_props()


    #calculate and sort all inter-particles distances
    def calc_particles_force_distances(self, sort_d = False):

        import itertools as itt

        #checking the distances between the particles
        dtype = [('i', int), ('j', int), ('d', float)]
        self.d_matrix = [ (i, j, 0) for i, j in list(itt.combinations(range(self.n_particles), 2)) ]
        self.d_matrix = np.array(self.d_matrix, dtype=dtype)

        for k in range(len(self.d_matrix)):
            
            #calculating the distance between the all pair of particles
            i, j = self.d_matrix[k][0] , self.d_matrix[k][1]
            d = get_particles_distance(self.p_props.r(i), self.p_props.r(j))
            self.d_matrix[k] = ( i, j, d )

        if sort_d is True:
            #organizing by crescent value of d
            self.d_matrix = np.sort(self.d_matrix, order='d')


    #updating accelerations
    def update_a(self):
        
        if self.potential == 'leonard_jones':

            for i, j, d in self.d_matrix:
                f_i_j = get_force_between_particles(d = d, r = self.p_props.r(i) - self.p_props.r(j), potential = 'leonard_jones')
                a_i_j = f_i_j / self.p_props.radius(i)
                a_j_i = - f_i_j / self.p_props.radius(j)
                self.p_props.a(i, update = True, new_value = self.p_props.a(i).copy() + a_i_j)
                self.p_props.a(j, update = True, new_value = self.p_props.a(j).copy() + a_j_i)


    #updating velocities
    def update_v(self):
            
            #checking the next ball positions
            r_next = self.p_props.r().copy() + (self.p_props.v().copy() * self.dt)

            #new velocities
            v_ = self.p_props.v().copy()

            #scanning the distance matrix
            for i, j, d in self.d_matrix:
                        
                #distance in the next dt
                d_next = get_particles_distance(r_next[i], r_next[j])
                m_i, m_j = self.p_props.radius(i), self.p_props.radius(j)

                if check_overlap(r_i = self.p_props.r(i), r_j = self.p_props.r(j), radius_i = m_i, radius_j = m_j) == True and ( d_next < d ):
                    
                    M = m_i + m_j
                    v_[i] = self.p_props.v(i) - ( ( ( 2 * m_j ) / M ) * ( np.dot( self.p_props.v(i) - self.p_props.v(j) , self.p_props.r(i) - self.p_props.r(j) ) * ( self.p_props.r(i) - self.p_props.r(j) ) ) / d**2 )
                    v_[j] = self.p_props.v(j) - ( ( ( 2 * m_i ) / M ) * ( np.dot( self.p_props.v(j) - self.p_props.v(i) , self.p_props.r(j) - self.p_props.r(i) ) * ( self.p_props.r(j) - self.p_props.r(i) ) ) / d**2 )
                
                #checking wall collision
                for i in range(self.n_particles):
                    
                    cond_wall_xmin = self.p_props.x(i) - self.p_props.radius(i) <= 0
                    cond_wall_v_xmin = v_[i, 0] < 0
                    cond_wall_xmax = self.p_props.x(i) + self.p_props.radius(i) >= self.box_width
                    cond_wall_v_xmax = v_[i, 0] > 0

                    if (False not in (cond_wall_xmax, cond_wall_v_xmax) ) or (False not in (cond_wall_xmin, cond_wall_v_xmin) ):
                        v_[i, 0] = ( v_[i, 0] * -1 )

                    cond_wall_ymin = self.p_props.y(i) - self.p_props.radius(i) <= 0
                    cond_wall_v_ymin = v_[i, 1] < 0
                    cond_wall_ymax = self.p_props.y(i) + self.p_props.radius(i) >= self.box_width
                    cond_wall_v_ymax = v_[i, 1] > 0

                    if (False not in (cond_wall_ymax, cond_wall_v_ymax) ) or (False not in (cond_wall_ymin, cond_wall_v_ymin) ):
                        v_[i, 1] = ( v_[i, 1] * -1 )

            #updating velocities as a function of acceleration
            v_ = v_ + (self.p_props.a().copy() * self.dt)
            self.p_props.v(update = True, new_value= v_)


    #updating positions
    def update_r(self):

        new_value = self.p_props.r().copy() + (self.p_props.v().copy() * self.dt)
        self.p_props.r(update= True, new_value = new_value)


    #maxwell-boltzmann distribution
    def maxwell_boltzmann(self, v, mass = None, T = None, factor = 5):
        
        k = 1 #here k is one to put the masses in a proper scale
        a = ( self.mass / ( 2 * k * self.temperature) )
        return factor * 4 * math.pi * ( a / math.pi)**(3/2) * (v**2) * math.exp( - a * (v**2) )


    #saving figures
    def save_fig(self, file_counter):
        
        vels = np.linspace(0, 600, 100)
        freq, bins = self.p_props.hist('v', bins = vels)
        freq = freq / self.n_particles
        
        fig, axes = plt.subplots(1, 2)
        plot_particles(axes[0], self.p_props, self.box_width, self.box_height, file_counter)
        plot_v_histogram(axes[1], freq, bins, hist_x_limit = 600, hist_y_limit = 0.1)
        
        #f = [ self.maxwell_boltzmann(v, mass = self.mass, T = self.temperature) for v in vels]
        #axes[1].plot(vels, f, color='red')

        fig.savefig(self.save_path + r'/p_{}.png'.format(correct_counter(file_counter)) )
        plt.close()



    #building gifs
    def build_gif(self):
        
        gif_name = 'movie'
        filenames = os.listdir(self.save_path)
        filenames = [self.save_path + f'/{filename}' for filename in filenames]
        filenames.sort()
        with imageio.get_writer(self.save_path + f'/{gif_name}.gif', mode='I') as writer:
            for filename in filenames:
                if filename[-3:] in ('png', 'jpg'):
                    image = imageio.v3.imread(filename)
                    writer.append_data(image)

        #Remove files
        #for filename in set(filenames):
        #    if filename[-3:] != 'gif':
        #        os.remove(filename)


    #run function
    def run(self, n_iterations = 100):
        
        if not os.path.exists(self.save_path):
            os.makedirs(self.save_path)

        #getting initial values
        self.get_initial_conditions()

        counter = 0
        while counter < n_iterations:
            self.calc_particles_force_distances()
            self.update_a()
            self.update_v()
            self.update_r()
            self.save_fig(counter)
            counter += 1
            print('iteration n ', counter)
        
        self.build_gif()
        #self.log.close()