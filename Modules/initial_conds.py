import numpy as np
from functions import check_wall_proximity, check_overlap

class initial_conds(object):

    def __init__(self, 
                 p_props = None,
                 radius_input = None,                  
                 n_particles = None,
                 n_dimension = None,
                 box_width = None,
                 box_height = None,
                 random_val_for_r = True,
                 random_val_for_v = True,
                 max_velocity = 10,
                 random_val_for_a = True,
                 max_accel = 10
                 ):
        
        self.radius_input = radius_input
        self.p_props = p_props
        self.n_particles = n_particles
        self.n_dimension = n_dimension
        self.max_velocity = max_velocity
        self.max_accel = max_accel
        self.box_width = box_width
        self.box_height = box_height
        self.random_val_for_r = random_val_for_r
        self.random_val_for_v = random_val_for_v
        self.max_velocity = max_velocity
        self.random_val_for_a = random_val_for_a
        self.max_accel = max_accel
        



    #getting initial radius
    def get_balls_radius(self):
        #self.log.write('> getting ball_radius' + '\n')
        #when radii were given for each particle
        if self.n_particles == len(self.radius_input):
            for i in range(len(self.p_props)):
                self.p_props.radius(i, update = True, new_value = self.radius_input[i])

        #for random radii
        else:
            new_value = np.random.uniform(min(self.radius_input), max(self.radius_input), size=self.n_particles)
            self.p_props.radius(update = True, new_value = new_value)

        #self.log.write(str(self.radius()) + '\n')


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
                    
                    if check_wall_proximity(box_width = self.box_width, box_height = self.box_height, r = new_r, radius = self.p_props.radius(i) , factor = 0.05 ) == False:
                        continue

                    #checking if ovelaps other particles
                    overlaped_particle = False
                    for j in range(len(self.p_props.p_props)):
                        #new random pos must be different from those already insertetd in the table and must be compared to those line in rvva table that were already fullfilled (no zero)
                        if ( new_r[0] != self.p_props.x(j) and self.p_props.x(j) != 0 )  and ( new_r[1] != self.p_props.y(j) and self.p_props.y(j) != 0 ):
                            if check_overlap(r_i = new_r, r_j = self.p_props.r(j), radius_i = self.p_props.radius(i), radius_j = self.p_props.radius(j), factor = 0.2) == True:
                                overlaped_particle = True
                                break
                        
                    #if not overlapped
                    if overlaped_particle is False:
                        self.p_props.r(i, update = True, new_value = new_r)                        
                        #self.log.write(f'> added particle {i} to the box at {new_r} \n')
                        added_to_box = True
                    
                    if added_to_box is True or counter > attempts:
                        break

            if self.p_props.r().min() > 0 :
                #self.log.write('> Succefful positioning of all particles.' + '\n')
                pass

            elif counter > attempts:
                print('> Error trying to add more particles. The box is alreay full. Change the number of particles or min-max particles radii.')
                #self.log.write('> Error trying to add more particles. The box is alreay full. Change the number of particles or min-max particles radii.' + '\n')
                break


    #geting v_0
    def balls_initial_velocities(self, random_val = False, val = (None, None, None)):

        if random_val is True:
            v = np.random.uniform( -self.max_velocity, self.max_velocity, size=(self.n_particles, self.n_dimension) )
            self.p_props.v(update = True, new_value = v)
        else:
            if self.n_dimension == 2:
                v = np.zeros((self.n_particles, self.n_dimension))
                v[:, 0] = np.repeat(val[0], self.n_particles)
                v[:, 1] = np.repeat(val[1], self.n_particles)
                self.p_props.v(update = True, new_value = v)
            elif self.n_dimension == 3:
                v = np.zeros((self.n_particles, self.n_dimension))
                v[:, 0] = np.repeat(val[0], self.n_particles)
                v[:, 1] = np.repeat(val[1], self.n_particles)
                v[:, 2] = np.repeat(val[2], self.n_particles)
                self.p_props.v(update = True, new_value = v)


    #geting a (constant)
    def balls_initial_accelerations(self, random_val = False, val = (None, None, None) ):
        
        if random_val is True:
            a = np.random.uniform(-self.max_accel, self.max_accel, size=(self.n_particles, self.n_dimension))
            self.p_props.a(update = True, new_value = a)
        else:
            if self.n_dimension == 2:
                a = np.zeros((self.n_particles, self.n_dimension))
                a[:, 0] = np.repeat(val[0], self.n_particles)
                a[:, 1] = np.repeat(val[1], self.n_particles)
                self.p_props.a(update = True, new_value = a)
            elif self.n_dimension == 3:
                a = np.zeros((self.n_particles, self.n_dimension))
                a[:, 0] = np.repeat(val[0], self.n_particles)
                a[:, 1] = np.repeat(val[1], self.n_particles)
                a[:, 2] = np.repeat(val[2], self.n_particles)
                self.p_props.a(update = True, new_value = a)

    
    #returning the modified p_props
    def update_p_props(self):
        self.get_balls_radius()
        self.balls_initial_positions(random_val = self.random_val_for_r)
        self.balls_initial_velocities(random_val = self.random_val_for_v)
        self.balls_initial_accelerations(random_val = self.random_val_for_a, val = (0, 0, 0) )
        
        return self.p_props
