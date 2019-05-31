# Functions to plot the tests of random_generator
from matplotlib.path import Path

class Plotter(object):
    '''
    Class to create a plot containing the area, the segments, and color them
    to identify which one of them are inside the area, which intercepts the
    boundary and other things like that.
    '''

    def __init__(self):
        self._boundary_path = None
        self._fibers = {}


    def define_boundary(vertices):
        try:
            n = len(vertices) - 2
            p = Path(vertices)
            codes = [p.MOVETO]
            codes += n * [p.LINETO]
            codes += [p.CLOSEPOLY]
            p.codes = codes
            self._boundary_path = {'path':p,'color':'black','lw':2.0}
        except:
            raise('Invalid data for defining boundary!')
            

