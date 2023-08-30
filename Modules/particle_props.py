import numpy as np

#getting a class for the matrix particles properties in the format [x, y, z, radius, vx, vy, vz, ax, ay, az]
class particle_props(object):


    def __init__(self, n_particles, n_dimension):
        
        self.p_props = np.zeros((n_particles, 3 + 1 + 3 + 3))
        self.n_dimension = n_dimension

    #functions for handling the rrva matrix properties
    #get positions (r)
    def r(self, *i:int, update:bool = False, new_value:np.array = None):
        if self.n_dimension == 2:
            if update == True:
                if len(i) == 1:
                    self.p_props[i[0], 0:2] = new_value.copy()
                else:
                    self.p_props[:, 0:2] = new_value.copy()
            
            else:
                if len(i) == 1:
                    return self.p_props[i[0], 0:2]
                else:
                    return self.p_props[:, 0:2]
        
        if self.n_dimension == 3:
            if update == True:
                if len(i) == 1:
                    self.p_props[i[0], 0:3] = new_value.copy()
                else:
                    self.p_props[:, 0:3] = new_value.copy()
                
            else:
                if len(i) == 1:
                    return self.p_props[i[0], 0:3]
                else:
                    return self.p_props[:, 0:3]

    #get x positions
    def x(self, *i:int):
        if self.n_dimension == 2:
            if len(i) == 1:
                return self.p_props[i[0], 0]
            else:
                return self.p_props[:, 0]
        if self.n_dimension == 3:    
            if len(i) == 1:
                return self.p_props[i[0], 0]
            else:
                return self.p_props[:, 0]

    #get y positions
    def y(self, *i:int):
        if self.n_dimension == 2:
            if len(i) == 1:
                return self.p_props[i[0], 1]
            else:
                return self.p_props[:, 1]
        if self.n_dimension == 3:    
            if len(i) == 1:
                return self.p_props[i[0], 1]
            else:
                return self.p_props[:, 1]

    #get z positions
    def z(self, *i:int):
        if self.n_dimension == 2:
            if len(i) == 1:
                return self.p_props[i[0], 2]
            else:
                return self.p_props[:, 2]
        if self.n_dimension == 3:    
            if len(i) == 1:
                return self.p_props[i[0], 2]
            else:
                return self.p_props[:, 2]
    
    #get radius
    def radius(self, *i:int, update:bool = False, new_value:np.array = None):
        if update == True:
            if len(i) == 1:
                self.p_props[i[0], 3] = new_value.copy()
            else:
                self.p_props[:, 3] = new_value.copy()
        
        else:
            if len(i) == 1:
                return self.p_props[i[0], 3]
            else:
                return self.p_props[:, 3]
    
    #get velocities
    def v(self, *i:int, update:bool = False, new_value:np.array = None):
        if self.n_dimension == 2:
            if update == True:
                if len(i) == 1:
                    self.p_props[i[0], 4:6] = new_value.copy()
                else:
                    self.p_props[:, 4:6] = new_value.copy()
            
            else:
                if len(i) == 1:
                    return self.p_props[i[0], 4:6]
                else:
                    return self.p_props[:, 4:6]
        
        if self.n_dimension == 3:
            if update == True:
                if len(i) == 1:
                    self.p_props[i[0], 4:7] = new_value.copy()
                else:
                    self.p_props[:, 4:7] = new_value.copy()
                
            else:
                if len(i) == 1:
                    return self.p_props[i[0], 4:7]
                else:
                    return self.p_props[:, 4:7]

    #get acelerations
    def a(self, *i:int, update:bool = False, new_value:np.array = None):
        if self.n_dimension == 2:
            if update == True:
                if len(i) == 1:
                    self.p_props[i[0], 7:9] = new_value.copy()
                else:
                    self.p_props[:, 7:9] = new_value.copy()
            
            else:
                if len(i) == 1:
                    return self.p_props[i[0], 7:9]
                else:
                    return self.p_props[:, 7:9]
        
        if self.n_dimension == 3:
            if update == True:
                if len(i) == 1:
                    self.p_props[i[0], 7:10] = new_value.copy()
                else:
                    self.p_props[:, 7:10] = new_value.copy()
                
            else:
                if len(i) == 1:
                    return self.p_props[i[0], 7:10]
                else:
                    return self.p_props[:, 7:10]


    #get vectors magnitude
    def get_magnitude(self, vectors: np.array):
        return np.apply_along_axis(np.linalg.norm, 1, vectors)


    #get histograms
    def hist(self, var, bins = 10, range = None):

        if var.lower() == 'r':
            mags = self.get_magnitude( self.r() )
            return np.histogram(mags, bins = bins, range = range)
        
        elif var.lower() == 'v':
            mags = self.get_magnitude( self.v() )
            return np.histogram(mags, bins = bins, range = range)

        elif var.lower() == 'a':
            mags = self.get_magnitude( self.a() )
            return np.histogram(mags, bins = bins, range = range)