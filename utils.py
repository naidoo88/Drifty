import numpy as np
import matplotlib.pylab as plt
#import boost_histogram as bh

###################################################
# Set up some global plot styling:                #
###################################################
plt.style.use('seaborn')
plt.rcParams['axes.grid'] = True
plt.rcParams['axes.formatter.limits'] = (-3,3)
plt.rcParams['axes.formatter.use_mathtext'] = True
plt.rcParams['font.size']= 16
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
# plt.rcParams['image.cmap'] = 'jet'
# plt.rcParams['figure.figsize'] = 8, 6
# plt.rcParams['legend.fontsize']= 12


###################################################
# Functions for plotting boost_histogram objects: #
###################################################

# boost_histogram is a histogramming library - NOT a plotting one!
# Therefore a couple of plotting functions will be useful.

def plot2D_bh(ax, h, cmap='viridis'):
    '''
    Function to plot 2D histograms.
    * Uses pylab.pcolormesh. Gets passed:
       -- the "quadrilateral corners" or bin-edges.
       -- the count-value of each bin.
       -- uses Viridis colormap unless specified when function is called.
    * A colorbar is drawn horizontally under the pcolormesh.
    '''
    pc = ax.pcolormesh(*h.axes.edges.T, h.view().T, cmap=cmap)
    plt.gcf().colorbar(pc, ax=ax, aspect=30, orientation="horizontal")
    return pc


def plot1D_bh(ax, h, color="", filled=True, alphafill=0.6, lparams={}, fparams={}):
    '''
    Function to plot 1D histograms.
    * Uses pylab.step. Gets passed:
       -- bin-edges, minus the last edge in order to match value dimensions.
       -- the count-value of each bin.
    * Color can be specified, or function simply uses the next item in colormap
    * To match a more typical histogram aesthetic, a vertical line is drawn at 
      either end of the step plot.
    * By default, the area under the step-plot is filled.  Can be set False.
    * Can pass any **kwarg dict via lparams/fparams for line/fill respectively.
    '''
    if not color:
        color = next(ax._get_lines.prop_cycler)['color']
    st = ax.step(h.axes[0].edges[:-1], h.view(), where='post', color=color, **lparams)
    ax.vlines(h.axes[0].edges[0], 0, h[0], color=color)
    ax.vlines(h.axes[0].edges[-1-1], 0, h[-1], color=color)
    if filled:
        ax.fill_between(h.axes[0].edges[:-1], 0, h,
            alpha=alphafill, 
            step='post', 
            color = color,
            **fparams
        )
    return st



###################################################
# Functions for fitting distributions with scipy: #
###################################################
# Define a few functions for the purpose of fitting

def gaussian(x, A, mu, sig):
    '''A Standard Gaussian'''
    return A * np.exp( - ((x - mu) / sig) ** 2)

def gaussian_offset(x, A, mu, sig, const):
    '''A Gaussian with a constant offset'''
    return A * np.exp( - ((x - mu) / sig) ** 2) + const
    
def linear(x, m, c):
    '''Simple linear function'''
    return m*x + c

def pol2(x, k2, k1, c):
    '''A second-order polynomial'''
    return k2*x**2 + k1*x + c

def pol3(x, k3, k2, k1, c):
    '''A third-order polynomial'''
    return k3*x**3 + k2*x**2 + k1*x + c

def pol4(x, k4, k3, k2, k1, c):
    '''A fourth-order polynomial'''
    return k4*x**4 + k3*x**3 + k2*x**2 + k1*x + c

def calc_chisq(f, fit_params, x, y_obs, sigma = 1):
    '''
    scipy.Optimise does not return Chi-Sq, but we want it in order to assess 
    goodness of fit.  Here we manually define the calculation, and return a 
    tuple containing the chi-sq, and the number of degrees of freedom in our 
    model.
    '''
    x = np.array(x)
    y_pred = np.array(f(x, *fit_params))
    y_obs = np.array(y_obs)
    sigma = np.array(sigma)
    X2 = np.sum(((y_obs-y_pred)/sigma)**2)
    dof = len(x) - len(fit_params)
    return (X2, dof)


