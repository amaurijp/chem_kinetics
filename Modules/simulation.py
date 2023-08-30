import os
import itertools as itt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib import animation
from itertools import combinations
import imageio
import math


from particle_props import particle_props

class simulation(particle_props):
    
    #instantiation
    def __init__(self, dt = 0.1, n_particles = 2, radius_input = [1, 2], max_velocity = None, mass= None, temperature = None, 
                 box_width = None, box_height = None, n_dimension = 2, root_path = None, save_path = None):
        
        self.n_particles = n_particles
        self.n_dimension = n_dimension
        self.radius_input = radius_input
        self.dt = dt
        self.max_velocity = max_velocity
        self.mass = mass
        self.temperature = temperature
        self.box_width = box_width
        self.box_height = box_height
                
        #calculating a mean value for mass (if more than one is given) - mass will be equivalent to radius in this model
        self.mean_radius = np.mean(np.array(self.radius_input))
        
        #diretorio
        self.save_path = save_path
        self.cwd = root_path
        
        #log_file
        #self.log = open(self.cwd + r'/log.txt', 'w')

        #instance of p_props matrix
        particle_props.__init__(self, n_particles, n_dimension)


    #getting initial radius
    def get_balls_radius(self):
        #self.log.write('> getting ball_radius' + '\n')
        #when radii were given for each particle
        if self.n_particles == len(self.radius_input):
            for i in range(len(self.p_props)):
                self.radius(i, update = True, new_value = self.radius_input[i])

        #for random radii
        else:
            new_value = np.random.uniform(min(self.radius_input), max(self.radius_input), size=self.n_particles)
            self.radius(update = True, new_value = new_value)

        #self.log.write(str(self.radius()) + '\n')


    #checking walls proximity
    def check_wall_proximity(self, r = None, radius = None, factor = 0):
        cond1 = ( r[0] - radius ) > self.box_width * factor
        cond2 = ( r[0] + radius ) < self.box_width -  ( self.box_width * factor )
        cond3 = ( r[1] - radius ) > self.box_height * factor
        cond4 = ( r[1] + radius ) < self.box_height - ( self.box_height * factor )
        #self.log.write('> checking wall proximity:\n')
        #self.log.write(f'  cond1: x - radius = {r[0] - radius} >  box_width * factor = {self.box_width * factor} \n')
        #self.log.write(f'  cond2: x + radius = {r[0] + radius} <  box_width - box_width * factor = {self.box_width - (self.box_width * factor)} \n')
        #self.log.write(f'  cond1: y - radius = {r[1] - radius} >  box_height * factor = {self.box_height * factor} \n')
        #self.log.write(f'  cond2: y + radius = {r[1] + radius} <  box_height - box_height * factor = {self.box_height - (self.box_height * factor)} \n')
        #self.log.write('  (cond1, cond2, cond3, cond4): ' + str((cond1, cond2, cond3, cond4)) + '\n')
        return False not in (cond1, cond2, cond3, cond4)


    #checking balls distance
    def get_particles_distance(self, i, j, r_i, r_j):
        d = np.linalg.norm(r_i - r_j)
        #self.log.write(f'> calculating distance for particles {i}, {j}: r_i = {r_i} , r_j = {r_j}, d = {d} \n')
        return d


    #function for checking overlap
    def check_overlap(self, i = None, j = None, r_i = None, r_j = None, radius_i = None, radius_j = None, factor = 0):
        d = self.get_particles_distance(i, j, r_i, r_j)
        #self.log.write('> checking overlaping between particles: \n')
        #self.log.write(f'  cond: radius_i = {radius_i}, radius_j = {radius_j}, d = {d}, radius_i + radius_j = {radius_i + radius_j} \n')
        return d <= ( radius_i + radius_j ) + ( ( radius_i + radius_j ) * factor )


    #geting r_0
    def balls_initial_positions(self, random_val = True):
        
        #self.log.write(f'> getting balls initial positions (random values = {random_val}) \n')

        #positioning each particle
        attempts = 10000
        for i in range(self.n_particles):
            
            #random positions
            added_to_box = False
            if random_val is True:

                counter = 0
                while True:
                    
                    new_r = np.random.uniform(size=2)
                    new_r[0] = new_r[0] * self.box_width #pos x in the box
                    new_r[1] = new_r[1] * self.box_height #pos y in the box
                    
                    counter += 1
                    #self.log.write(f'> trying position: x, y = {new_r} , radius = {self.radius(i)} \n')
                    
                    if self.check_wall_proximity(r = new_r, radius = self.radius(i) , factor = 0.05 ) == False:
                        continue

                    #checking if ovelaps other particles
                    overlaped_particle = False
                    for j in range(len(self.p_props)):
                        #new random pos must be different from those already insertetd in the table and must be compared to those line in rvva table that were already fullfilled (no zero)
                        if ( new_r[0] != self.x(j) and self.x(j) != 0 )  and ( new_r[1] != self.y(j) and self.y(j) != 0 ):
                            if self.check_overlap(i = i, j = j, r_i = new_r, r_j = self.r(j), radius_i = self.radius(i), radius_j = self.radius(j), factor = 0.2) == True:
                                overlaped_particle = True
                                break
                        
                    #if not overlapped
                    if overlaped_particle is False:
                        self.r(i, update = True, new_value = new_r)                        
                        #self.log.write(f'> added particle {i} to the box at {new_r} \n')
                        added_to_box = True
                    
                    if added_to_box is True or counter > attempts:
                        break

            if self.r().min() > 0 :
                #self.log.write('> Succefful positioning of all particles.' + '\n')
                pass

            elif counter > attempts:
                print('> Error trying to add more particles. The box is alreay full. Change the number of particles or min-max particles radii.')
                #self.log.write('> Error trying to add more particles. The box is alreay full. Change the number of particles or min-max particles radii.' + '\n')
                break


    #geting v_0
    def balls_initial_velocities(self, random_val = False, val = (None, None, None)):
        #self.log.write('> getting balls initial velocities \n')
        if random_val is True:
            v = np.random.uniform( -self.max_velocity, self.max_velocity, size=(self.n_particles, self.n_dimension) )
            self.v(update = True, new_value = v)
        else:
            if self.n_dimension == 2:
                v = np.zeros((self.n_particles, self.n_dimension))
                v[:, 0] = np.repeat(val[0], self.n_particles)
                v[:, 1] = np.repeat(val[1], self.n_particles)
                self.v(update = True, new_value = v)
            elif self.n_dimension == 3:
                v = np.zeros((self.n_particles, self.n_dimension))
                v[:, 0] = np.repeat(val[0], self.n_particles)
                v[:, 1] = np.repeat(val[1], self.n_particles)
                v[:, 2] = np.repeat(val[2], self.n_particles)
                self.v(update = True, new_value = v)
        
        #self.log.write(str(self.v()) + '\n')


    #geting a (constant)
    def balls_initial_accelerations(self, random_val = False, val = (None, None, None) ):
        #self.log.write('> getting balls acceleration \n')
        if random_val is True:
            a = np.random.uniform(-self.max_velocity, self.max_velocity, size=(self.n_particles, self.n_dimension))
            self.a(update = True, new_value = a)
        else:
            if self.n_dimension == 2:
                a = np.zeros((self.n_particles, self.n_dimension))
                a[:, 0] = np.repeat(val[0], self.n_particles)
                a[:, 1] = np.repeat(val[1], self.n_particles)
                self.a(update = True, new_value = a)
            elif self.n_dimension == 3:
                a = np.zeros((self.n_particles, self.n_dimension))
                a[:, 0] = np.repeat(val[0], self.n_particles)
                a[:, 1] = np.repeat(val[1], self.n_particles)
                a[:, 2] = np.repeat(val[2], self.n_particles)
                self.a(update = True, new_value = a)
        
        #self.log.write(str(self.a()) + '\n')


    #updating positions
    def update_r(self, copy_matrix = False):
        #self.log.write('> uptading positions' + '\n')
        #self.log.write(f'  r before\n {self.r()} \n')
        
        new_value = self.r().copy() + (self.v().copy() * self.dt)
        self.r(update= True, new_value = new_value)
        
        #self.log.write(f'  r after\n {self.r()} \n')


    #updating velocities
    def update_v(self, collision_type = 'elastic'):

        #updating velocities as a function of acceleration
        #self.log.write('> Uptading velocities \n')
        #self.log.write(f'  v before {self.v()} \n') 
        v_ = self.v().copy() + (self.a().copy() * self.dt)
        #self.log.write('  v after acceleration\n{v_}\n')
        self.v(update = True, new_value= v_)

        #checking the distances between the particles
        dtype = [('i', int), ('j', int), ('d', float)]
        d_sorted = [ (i, j, 0) for i, j in list(itt.combinations(range(self.n_particles), 2)) ]
        d_sorted = np.array(d_sorted, dtype=dtype)

        #calculating the distance between the all pair of particles
        for k in range(len(d_sorted)):
            i, j = d_sorted[k][0] , d_sorted[k][1]
            d_sorted[k] = ( i, j, self.get_particles_distance(i, j, self.r(i), self.r(j)) )
        
        #organizing by crescent value of d
        d_sorted = np.sort(d_sorted, order='d')
        #self.log.write('d_sorted' + '\n')
        #self.log.write(str(d_sorted) + '\n')
        
        #checking the next ball positions
        
        #self.log.write(f'  r_next before\n {self.r()} \n')
        r_next = self.r().copy() + (self.v().copy() * self.dt)
        #self.log.write(f'  r_next after\n {r_next} \n')

        #new velocities
        v_ = self.v().copy()

        #loop while d <= max value of summed radius
        d = d_sorted[0][2]
        index = 0
        while d <= ( 2 * self.radius().max() ):

            i, j = d_sorted[index][0], d_sorted[index][1]
        
            #distance in the next dt
            d_next = self.get_particles_distance(i, j, r_next[i], r_next[j])
            
            #self.log.write(f'> Checking collisions between {i} , {j} \n')
            #self.log.write(f'  r_i = {self.r(i)}, r_j = {self.r(j)} \n')
            #self.log.write(f'  d = {d}, d_next = {d_next}, radius_i + radius_j = {self.radius(i) + self.radius(j)} \n')
            
            m_i, m_j = self.radius(i), self.radius(j)

            if self.check_overlap(i = i, j = j, r_i = self.r(i), r_j = self.r(j), radius_i = m_i, radius_j = m_j) == True and ( d_next < d ):
                
                #self.log.write(f'  collision detected between particles i = {i}, j = {j} \n')

                if collision_type == 'elastic':

                    M = m_i + m_j
                    v_[i] = self.v(i) - ( ( ( 2 * m_j ) / M ) * ( np.dot( self.v(i) - self.v(j) , self.r(i) - self.r(j) ) * ( self.r(i) - self.r(j) ) ) / d**2 )
                    v_[j] = self.v(j) - ( ( ( 2 * m_i ) / M ) * ( np.dot( self.v(j) - self.v(i) , self.r(j) - self.r(i) ) * ( self.r(j) - self.r(i) ) ) / d**2 )
            
            index += 1
            if index < len(d_sorted):
                d = d_sorted[index][2]
            else:
                break
        
        #checking wall collision
        for i in range(self.n_particles):
            
            cond_wall_xmin = self.x(i) - self.radius(i) <= 0
            cond_wall_v_xmin = v_[i, 0] < 0
            cond_wall_xmax = self.x(i) + self.radius(i) >= self.box_width
            cond_wall_v_xmax = v_[i, 0] > 0

            if (False not in (cond_wall_xmax, cond_wall_v_xmax) ) or (False not in (cond_wall_xmin, cond_wall_v_xmin) ):
                v_[i, 0] = ( v_[i, 0] * -1 )

            cond_wall_ymin = self.y(i) - self.radius(i) <= 0
            cond_wall_v_ymin = v_[i, 1] < 0
            cond_wall_ymax = self.y(i) + self.radius(i) >= self.box_width
            cond_wall_v_ymax = v_[i, 1] > 0

            if (False not in (cond_wall_ymax, cond_wall_v_ymax) ) or (False not in (cond_wall_ymin, cond_wall_v_ymin) ):
                v_[i, 1] = ( v_[i, 1] * -1 )

        
        #updating velocities as a function of the collisions
        #self.log.write(f'  v after collisions\n{v_}\n')
        self.v(update = True, new_value= v_)


    #maxwell-boltzmann distribution
    def maxwell_boltzmann(self, v, mass = None, T = None, factor = 5):
        k = 1 #here k is one to put the masses in a proper scale
        a = ( self.mass / ( 2 * k * self.temperature) )
        return factor * 4 * math.pi * ( a / math.pi)**(3/2) * (v**2) * math.exp( - a * (v**2) )


    #saving figures
    def save_fig(self, file_counter, insert_ball_number = False, hist_x_limit = 100):
        
        fig, axes = plt.subplots(1, 2)
        
        for s in ['top','bottom','left','right']:
            axes[0].spines[s].set_linewidth(2)

        axes[0].set_title('frame ' + self.correct_counter(file_counter) )
        axes[0].set_aspect('equal', 'box')
        axes[0].set_xlim(0, self.box_width)
        axes[0].set_ylim(0, self.box_width)
        axes[0].xaxis.set_ticks([])
        axes[0].yaxis.set_ticks([])

        index = 0
        styles = {'edgecolor': 'C0', 'facecolor':'royalblue' , 'linewidth': None, 'fill': True}        
        for i in range(len(self.p_props)):
            x = self.x(i)
            y = self.y(i)
            radius = self.radius(i)
            circle = Circle(xy=(x, y), radius=radius, **styles)
            axes[0].add_patch(circle)
            if insert_ball_number is True:
                axes[0].annotate(index, (x - 3, y - 3))
            index += 1
        
        vels = np.linspace(0, 600, 100)

        freq, bins = self.hist('v', bins = vels)
        freq = freq / self.n_particles
        axes[1].set_title('v; temperature' )
        axes[1].stairs(freq, bins, fill = True)
        axes[1].yaxis.set_ticks([])
        axes[1].set_xlim(0, hist_x_limit)
        axes[1].set_ylim(0, 0.1)
        
        f = [ self.maxwell_boltzmann(v, mass = self.mass, T = self.temperature) for v in vels]
        axes[1].plot(vels, f, color='red')

        fig.savefig(self.save_path + r'/p_{}.png'.format(self.correct_counter(file_counter)) )
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


    #correct numbering for filenames
    def correct_counter(self, number_int):

        if number_int < 10:
            return str('000000' + str(number_int))
        elif 10 <= number_int < 100:
            return str('00000' + str(number_int))
        elif 100 <= number_int < 1000:
            return str('0000' + str(number_int))
        elif 1000 <= number_int < 10000:
            return str('000' + str(number_int))
        elif 10000 <= number_int < 100000:
            return str('00' + str(number_int))
        elif 100000 <= number_int < 1000000:
            return str('0' + str(number_int))
        elif 1000000 <= number_int < 10000000:
            return str(str(number_int))

        

    #run function
    def run(self, n_iterations = 100):
        
        if not os.path.exists(self.save_path):
            os.makedirs(self.save_path)

        #getting initial values
        self.get_balls_radius()
        self.balls_initial_positions(random_val = True)
        self.balls_initial_velocities(random_val = True, val = (None, None, None))
        self.balls_initial_accelerations(random_val = False, val = (0, -100))

        counter = 0
        while counter < n_iterations:
            #self.log.write(f'\n\n\n> iteration {counter}\n')
            self.update_r()
            self.update_v()
            self.save_fig(counter, insert_ball_number = False, hist_x_limit= 600)
            counter += 1
            print('iteration n ', counter)
        
        self.build_gif()
        #self.log.close()