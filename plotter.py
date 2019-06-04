# Functions to plot the tests of random_generator
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import matplotlib.pyplot as plt


class Plotter(object):
    '''
    Class to create a plot containing the area, the segments, and color them
    to identify which one of them are inside the area, which intercepts the
    boundary and other things like that.
    '''

    def __init__(self):
        self._boundary_path = None
        self._fibers = {}


    def define_boundary(self,vertices):
        try:
            n = len(vertices) - 2
            p = Path(vertices)
            codes = [p.MOVETO]
            codes += n * [p.LINETO]
            codes += [p.CLOSEPOLY]
            p.codes = codes
            coords = p._vertices
            xmin = min(coords[:,0]);xmax = max(coords[:,0])
            ymin = min(coords[:,1]);ymax = max(coords[:,1])
            self._boundary_path = {'path':p,'color':'black','lw':2.0, 'limits':((xmin,ymin),(xmax,ymax))}
        except:
            raise('Invalid data for defining boundary!')
            
    def plot(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        # Adding boundary to the plot
        if self._boundary_path:
            p = self._boundary_path['path']
            color = self._boundary_path['color']
            lw = self._boundary_path['lw']
            ptc = PathPatch(p, facecolor='None', edgecolor=color, lw=lw)
            ax.add_patch(ptc)
            ((xmin,ymin),(xmax,ymax)) = self._boundary_path['limits']
            # Setting drawing limits
            lx = xmax - xmin; ly = ymax - ymin
            xm = 0.5 * (xmax + xmin); ym = 0.5 * (ymax + ymin)
            ax.set_xlim(xm-0.55*lx,xm+0.55*lx)
            ax.set_ylim(ym-0.55*ly,ym+0.55*ly)
            ax.set_aspect('equal')

        else:
            print('No boundary to plot! Exiting!')
            return None

        plt.show()


