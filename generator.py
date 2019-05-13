# Functions for random generator and for verifying if fibers are inside area.
# As the final aim is to generate fibers in 3D, the depth (z) is also considered

from numpy.random import random, randn
from math import pi
from numpy import cos as npcos
from numpy import sin as npsin
from numpy import array, ndarray


def generate_fiber_center(number_of_fibers, xlim, ylim, zlim):
    '''
    Returns a sequence os coordinates (x, y, z).
    '''
    xmin, xmax = xlim
    ymin, ymax = ylim
    zmin, zmax = zlim
    # xs
    qsi = random(number_of_fibers)
    xs = xmin * (1.0 - qsi) + xmax * qsi
    # ys
    qsi = random(number_of_fibers)
    ys = ymin * (1.0 - qsi) + ymax * qsi
    # zs
    qsi = random(number_of_fibers)
    zs = zmin * (1.0 - qsi) + zmax * qsi
    return xs, ys, zs

def generate_angles(number_of_fibers, beta_lim = (0.0,2.0*pi),
                    alpha_params = (30.0*pi/180.0, 15.0*pi/180.0)):
    '''
    Returns angles for generating fibers (in radians).
    
    Alpha is the angle around an axe in a plane
    perpendicular to the x-y plane;
    Beta is the angle around an axe in a plane 
    perpendicular to the x-z plane.

    The argument beta_lim indicates the minimum
    and maximum values of beta, in radians.
    
    The argument alpha_params indicates the mean value
    of alpha, in radians, and its stardard variation,
    assuming a normal distribution for such an angle.
    '''
    min_beta, max_beta = beta_lim
    mu, sigma = alpha_params
    alpha = sigma * randn(number_of_fibers) + mu
    qsi = random(number_of_fibers)
    beta = min_beta * (1.0 - qsi) + max_beta * qsi
    return alpha, beta

def compute_fiber_edges(xs, ys, zs, alphas, betas, Lf):
    '''
    Returns the coordinates of the begin and end of the fiber.
    ((xi, yi, zi),(xf, yf, zf)).
    - xs - x coordinates (middle of the fiber);
    - ys - y coordinates (middle of the fiber);
    - zs - z coordinates (middle of the fiber);
    - alphas and betas - angles (as describer in 'generate_angles');
    - Lf - Fiber lenght.
    '''
    dx = 0.5 * Lf * npcos(betas) * npsin(alphas)
    dy = 0.5 * Lf * npsin(betas) * npsin(alphas)
    dz = 0.5 * Lf * npcos(alphas)

    xi = xs + dx; yi = ys + dy; zi = zs + dz
    xf = xs - dx; yf = ys - dy; zf = zs - dz

    return ((xi, yi, zi), (xf, yf, zf))

def intercepts(segment, other_segment):
    '''
    Returns True if two segments intercept each other; returns False, otherwise).
    The segment is a tuple ((xi, yi), (yf, yf)); two different situations can be computed:
    1) both segment and other_segmente are tuples of tuple, with float values for x and y;
    2) one of the tuples contains numpy arrays instead of floats; in this case,
    the program returns an array with True or False, values, instead of a single
    boolean variable. This use recommend for several (thousands of) segments, for which
    the case 1) will be considerably less efficient.
    '''
    # Testing if input data is consistent (maybe this testing is not so smart... but...)
    try:
        ((xsi,ysi),(xsf,ysf)) = segment
        ((xosi,yosi),(xosf,yosf)) = other_segment
        if isinstance(xsi,float) and isinstance(xsf, float) and \
           isinstance(ysi,float) and isinstance(ysf, float):
            if isinstance(xosi, ndarray) and isinstance(xosf, ndarray) and \
               isinstance(yosi, ndarray) and isinstance(yosf, ndarray):
                pass
            elif isinstance(xosi, float) and isinstance(xosf, float) and \
                 isinstance(yosi, float) and isinstance(yosf, float):
                 pass
            else:
                raise Exception('Invalid values for input data!')
        elif isinstance(xsi, ndarray) and isinstance(xsf, ndarray) and \
             isinstance(ysi, ndarray) and isinstance(ysf, ndarray):
            if isinstance(xosi, float) and isinstance(xosf, float) and \
               isinstance(yosi, float) and isinstance(yosf, float):
                pass
            else:
                raise Exception('Invalid values for input data!')
    except:
        raise Exception('Invalid values for input data!')
    
    # if it gets here, data is OK!
    print("Data is consistent")



