# Functions for random generator and for verifying if fibers are inside area.
# As the final aim is to generate fibers in 3D, the depth (z) is also considered

from numpy.random import random, randn
from math import pi
from numpy import cos as npcos
from numpy import sin as npsin
from numpy import array, ndarray
from numpy import logical_and, logical_or
from numpy import abs as npabs
from matplotlib.path import Path


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

def intercepts(segment, other_segment, tol = 1.0e-12):
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
                case = "FA" # float-array
            elif isinstance(xosi, float) and isinstance(xosf, float) and \
                 isinstance(yosi, float) and isinstance(yosf, float):
                 case = "FF" # float-float
            else:
                raise Exception('Invalid values for input data!')
        elif isinstance(xsi, ndarray) and isinstance(xsf, ndarray) and \
             isinstance(ysi, ndarray) and isinstance(ysf, ndarray):
            if isinstance(xosi, float) and isinstance(xosf, float) and \
               isinstance(yosi, float) and isinstance(yosf, float):
                case = "FA" # float-array
            else:
                raise Exception('Invalid values for input data!')
    except:
        raise Exception('Invalid values for input data!')
    
    # if it gets here, data is OK!
    #print("Data is consistent")
    
    # Deriving the formulation

    # segment (xs, ys) for a given qsi (0<=qsi<=1)
    # xs(qsi) = xsi * (1 - qsi) + xsf * qsi
    # ys(qsi) = ysi * (1 - qsi) + ysf * qsi

    # other_segment (xos, yos) for a given eta (0<=eta<=1)
    # xos(eta) = xosi * (1 - eta) + xosf * eta
    # yos(eta) = yosi * (1 - eta) + yosf * eta
    
    # segment and other_segment cross each other when
    # xs(qsi) == xos(eta) [eq. (I)] 
    # and ys(qsi) == yos(eta) [eq. (II)]
    # if segment and other_segment are not parallel
    # and if 0<=qsi<=1 and 0<=eta<=1, they cross each other

    # Deriving eq. (I)
    # xsi * (1 - qsi) + xsf * qsi == xosi * (1 - eta) + xosf * eta
    # (xsf - xsi) * qsi - (xosf - xosi) * eta = (xosi - xsi)

    # Deriving eq. (II)
    # ysi * (1 - qsi) + ysf * qsi == yosi * (1 - eta) + yosf * eta
    # (ysf - ysi) * qsi - (yosf - yosi) * eta = (yosi - ysi)

    # The resulting system of equation is:

    # | (xsf - xsi)  -(xosf - xosi) |   | qsi |   | (xosi - xsi) |
    # |                             | * |     | = |              |
    # | (ysf - ysi)  -(yosf - yosi) |   | eta |   | (yosi - ysi) |

    # D = (xsf - xsi) * (yosi -yosf) - (ysf -ysi) * (xosi - xosf)
    # If the determinant D == 0, the segments are parallel segments
    # Changing term names:
    # A11 = (xsf - xsi); A12 = (xosi - xosf)
    # A21 = (ysf - ysi); A22 = (yosi - yosf)
    # B1 = (xosi - xsi); B2 = (yosi - ysi)
    # D = A11 * A22 - A12 * A21
    # qsi = (A22 * B1 - A12 * B2) / D (if D != 0)
    # eta = (A11 * B2 - A21 * B1) / D (if D != 0)

    A11 = (xsf - xsi); A12 = (xosi - xosf)
    A21 = (ysf - ysi); A22 = (yosi - yosf)
    B1 = (xosi - xsi); B2 = (yosi - ysi)
    D = A11 * A22 - A12 * A21

    if case == "FF": # Both segments are represented by float values
        if abs(D) < tol:
            # segment and segment are parallels, and can be collinear
            # However, this function is intended to be used after
            # testing if both edges are inside a given area. In this
            # situation, for the aimed purpose, if they are colinear
            # the status False describes better the situation to accept
            # the segment (one do not intercept the other)
            result = False
        else:
            qsi = (A22 * B1 - A12 * B2) / D # (if D != 0)
            eta = (A11 * B2 - A21 * B1) / D # (if D != 0)
        if (0.0 <= qsi <= 1.0) and (0.0 <= eta <= 1.0):
            result = True
        else:
            result = False
    elif case == "FA": # One of the segments are represented by numpy arrays
        # Identify which determinats are lower than the tolerance!
        prls = npabs(D) < tol # prls is an array of boolen indicating parallel lines
        # For the same reason indicated above, parallel lines return False status 
        # for the present function; one way to deal with this issue here is to modify
        # D such as the computations for this case result in a cross point far away
        # from the segment domain. In this case, a solution is to set D = tol
        # when abs(D) < tol
        D[prls] = tol    
        print('D->',D)
        # Computing the array of qsis and etas
        qsi = (A22 * B1 - A12 * B2) / D # (if D != 0)
        eta = (A11 * B2 - A21 * B1) / D # (if D != 0)
        # Testing interception
        qsis_ = logical_and(0.0<=qsi,qsi<=1.0)
        etas_ = logical_and(0.0<=eta,eta<=1.0)
        result = logical_and(qsis_, etas_)
        
    return result

    def points_inside_boundary(points, boundary_vertices):
        '''
        points - List of (x, y) coordinates of points;
        boundary_vertices - List  of (x, y) coordinates of the
        points defining the boundary of an area.
        Returns an numpy array of booleans; True if a given point
        is inside the defined area; False, otherwise.
        Example:
        points = ((0.0,0.0),(2.0,2.0)) # The first inside the area,
                                       # the second one outside.
        boundary_vertices = ((-1.0,-1.0),
                             ( 1.0,-1.0),
                             ( 1.0, 1.0),
                             (-1.0, 1.0),
                             (-1.0,-1.0))
        # For this case, the result is an array contains (True, False).
        '''
        try:
            n = len(points) - 1
            p = Path(boundary_vertices)
            codes = [p.MOVETO]
            codes += n * [p.LINETO]
            codes += [p.CLOSEPOLY]
            p.codes = codes
            result = p.contains_points(points)
        except:
            result = None
        return result








