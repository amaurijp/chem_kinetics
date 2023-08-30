import numpy as np
    

#function for checking overlap
def check_overlap(r_i = None, r_j = None, radius_i = None, radius_j = None, factor = 0):
    d = get_particles_distance(r_i, r_j)
    return d <= ( radius_i + radius_j ) + ( ( radius_i + radius_j ) * factor )


#checking walls proximity
def check_wall_proximity(box_width = None, box_height = None, r = None, radius = None, factor = 0):
    cond1 = ( r[0] - radius ) > box_width * factor
    cond2 = ( r[0] + radius ) < box_width -  ( box_width * factor )
    cond3 = ( r[1] - radius ) > box_height * factor
    cond4 = ( r[1] + radius ) < box_height - ( box_height * factor )

    return False not in (cond1, cond2, cond3, cond4)


#correct numbering for filenames
def correct_counter(number_int):

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


#calculating forces between particles
def get_force_between_particles(d:float = None, r:np.array = None, potential = 'leonard_jones'):
    
    #radius here is equivalent to mass, since we are using different particles sizes
    if potential == 'leonard_jones':
        k = 1000
        rm = 30
        return ( ( 12 * k / rm**2 ) * ( (rm / d)**14 - ( rm / d)**8 ) * r )


#checking balls distance
def get_particles_distance(r_i, r_j):
    d = np.linalg.norm(r_i - r_j)
    return d


#plot particles position
def plot_particles(ax, p_props, box_width, box_height, file_counter):

    from matplotlib.patches import Circle

    for s in ['top','bottom','left','right']:
        ax.spines[s].set_linewidth(1)

    ax.set_title('frame ' + correct_counter(file_counter) )
    ax.set_aspect('equal', 'box')
    ax.set_xlim(0, box_width)
    ax.set_ylim(0, box_height)
    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])

    styles = {'edgecolor': 'C0', 'facecolor':'royalblue' , 'linewidth': None, 'fill': True}        
    for i in range(len(p_props.p_props)):
        x = p_props.x(i)
        y = p_props.y(i)
        radius = p_props.radius(i)
        circle = Circle(xy=(x, y), radius=radius, **styles)
        ax.add_patch(circle)


#plot velocity histogram
def plot_v_histogram(ax, freq, bins, hist_x_limit = 1, hist_y_limit = 1):
    
    ax.set_title('v')
    ax.stairs(freq, bins, fill = True)
    ax.yaxis.set_ticks([])
    ax.set_xlim(0, hist_x_limit)
    ax.set_ylim(0, hist_y_limit)