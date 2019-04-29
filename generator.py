# Functions for random generator and for verifying if fibers are inside area.
# As the final aim is to generate fibers in 3D, the depth (z) is also considered

from numpy.random import random, randn
from math import pi


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
    dx = 0.5 * L

