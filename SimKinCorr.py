import numpy as np
import matplotlib.pylab as plt
import boost_histogram as bh
from scipy.optimize import curve_fit
import functools as ft
import copy
import utils as utl

class CalcCorrections:
    '''
    Class is designed to be used for calculating simulation-driven kinematic 
    corrections which correct for offsets caused by reconstruction algorithms.     

    Contains methods required to create a 2D histogram of delta/diff vs recon, 
    split it into N-slices and create N projections, fit the means of each slice
    and determine the distribution of the mean correction required across the 
    kinematic range provided. 
    
    When an instance of the class is initialised, a 2D histogram is generated 
    and filled.  Slices are then created as instances of a Slice-class object, 
    which fits each slice and stores slice-specific data.  See Slice-class docs 
    for more info.  The fitted means for each slice are then plotted against the 
    center position of the slices, and this is fitted with a weighted least-
    squares to produce a relationship which can be used for the purposes of 
    applying whatever correction is necessary to correct the systematic shift 
    the reconstruction introduces.  All of the produced distributions are 
    plotted on a single figure such that the end-user can easilly determine the 
    quality of all of the fits and decide whether or not any adjustments need 
    to be made.

    The process is largely automated by means of the process method, which is 
    called in the __init__ method.  Although this is unconventional, it is 
    functional for the design purposes of the class.  That said, it is expected 
    that some tweaking will be required.
    
    Therefore the class contains several methods to "tweak" the fits once the
    algorithms have made their best estimation.  
    
    Individual slices can:
    -- Have their initial fit parameters hard coded
    -- Have their fit range changed
    -- Be rebinned if a particular slice has slow stats
    -- If you want to revert the changes you have made, the slice can be 
       restored.
    -- If a slice cannot be 'made to cooperate' it can be vetoed, meaning it 
       will not be included in the final fit of the means.  
       NB: slices can also be 'unveto'd'.

    Once changes have been made to an individual slice, the fitting of the means
    is automatically updated and the results displayed in the figure.

    Instances of the class must be supplied:
        -- data, a Pandas Dataframe containing the data which will be analysed
        -- setup, a Dictionary containing information regarding the variable 
        which is to be corrected.  
    User may also supply:
        -- cuts, a list of pandas df filters which are combined into a single 
        boolean series below.
        -- n_slices, how many slices/projections you would like to produce

    The setup dictionary should have the format:
        -- 'var': The name of the variable in the dataframe we wish to 
                    correct.
        -- 'v_bins_range': A tuple which contains (# of bins, low range lim, 
                            high range lim) for 'var'.
        -- 'v_label': A plot label for 'var' which will be used in the 
                        figure.
        -- 'delta': The name of the difference of the target variable we 
                    are correcting.
        -- 'd_bins_range': A tuple which contains (# of bins, low range lim, 
                            high range lim) for 'delta'.
        -- 'd_label': A plot label for 'delta' which will be used in the 
                        figure.
    All fields in setup must be provided.
    '''
    
    def __init__(self, data, setup, cuts=None, n_slices=10):   
        self.hist2D = None # 2D bh.histo
        self.n_slices = n_slices # No. of slices
        self.slices = [] # Slice-class objects (projections & properties)
        self.output_fig = None # arranged figure to show output
        self.output_ax = None # axes/subplots of output_fig
        self.fit_model = utl.linear # model for fit of means:- default linear
        self.fit_vals = None # results of means fit
        self.fit_errs = None # errors of means fit results
        self.fit_cov = None # covariance for means fit
        self.chisqdof = None # tuple containing chi-sq and dof for means model
        self.cached_slices = [] # keeps copy of original slice when rebinning

        self.process(data, setup, cuts)
        
    def generate_2Dhisto(self, data, setup, cuts):
        '''
        Using the provided dataframe, and parameters from the setup dictionary, 
        generates the 2D boost_histogram object and fills them.

        If a list of cuts is provided, these are combined such that they
        produce a single boolean result for each event/entry in the dataframe.
        '''
        h = bh.Histogram(
            bh.axis.Regular(*setup['v_bins_range']), # x-axis
            bh.axis.Regular(*setup['d_bins_range']) # y-axis
        )
        if cuts is not None:
            # NB: functools.reduce accepts a function as input, and repeats the 
            # function iteratively over the range of arguments.
            # ie: if we pass values [a,b,c,d], reduce evaluates:
            #     * a&b
            #     * (a&b)&c
            #     * (a&b&c)&d 
            # Below we pass every item in the cuts list.
            cut = ft.reduce(lambda x, y: x & y, cuts[:])
            h.fill((data[setup['var']][cut].values), 
                (data[setup['delta']][cut].values)
            )    
        else:
            h.fill((data[setup['var']].values), (data[setup['delta']].values))    
        
        self.hist2D = h


    def create_slices(self):
        '''
        Takes 2D histogram, and slices it N-times into slices in 'x', and 
        projects the slices in 'y'.  Each slice/projection is initialised as a
        "Slice" object, which contains all relevant data for the slice. 

        If #2D-histo-xbins/#slices is indivisible, the modulo are abitrarily 
        split amongst the first bins.
        '''
        total_xbins = len(self.hist2D.axes[0]) 
        #first max/min returns tuple containing bin edges, hence we repeat.
        total_xrange = (min(min(self.hist2D.axes[0])),
            max(max(self.hist2D.axes[0]))
        )
        bin_len = (total_xrange[1] - total_xrange[0])/total_xbins

        #calc quotient/mod to handle indivisable n_xbins/n_slices
        bin_q,bin_mod = divmod(total_xbins, self.n_slices)
        lim_low_buff = 0
        for n in range(self.n_slices):
            if bin_mod == 0:
                slice_bins = bin_q
            else: 
                # if indivisible, (abitrarily) add 1 to quotient and subtract 1 
                # from mod until remaining bins are divisible.
                slice_bins = bin_q+1
                bin_mod = bin_mod-1
            
            slice_len = slice_bins*bin_len
            
            lim_low = lim_low_buff
            lim_high = lim_low+slice_bins
            lim_low_buff = lim_high
            
            slice_lim = (lim_low,lim_high)
            slice_centre = total_xrange[0] + lim_low*bin_len + (slice_len/2)
            
            # sums between limits in bins on x-axis, for all bins in y-axis.
            # -- see boost_histogram docs
            slice_projection = self.hist2D[lim_low:lim_high:sum, :]
            
            self.slices.append(Slice(N=n, 
                                     proj=slice_projection, 
                                     lims=slice_lim, 
                                     centre=slice_centre)
            )

            
    def output_grid(self):
        '''
        Generate the output figure for displaying all results.
        This is achieved using matplotlib.gridspec.  The below code produces:
            * Left/Right half split:
                - Left: halved vertically for 2D histo, plotted fitted means.
                - Right: Nx2 grid to display slice/projection distributions
        '''
        self.output_fig = plt.figure(figsize=(12,8), constrained_layout=True)
        grd = self.output_fig.add_gridspec(1,2, wspace=0.0, hspace=0.0)
        grd_left = grd[0].subgridspec(2,1)
        self.output_fig.add_subplot(grd_left[0,0])
        self.output_fig.add_subplot(grd_left[1,0])
        
        if self.n_slices % 2 != 0:
            # if off number of slices being used, we create a grid of N+1.
            N = np.int((self.n_slices + 1)/2)
        else:
            N = np.int(self.n_slices/2)
        grd_right = grd[1].subgridspec(N,2)
        for x in range(N):
            for y in range(2):
                self.output_fig.add_subplot(grd_right[x,y])
                
        self.output_ax = self.output_fig.get_axes()
        
        
    def fit_means(self):
        '''
        Take fitted means, and fit them to ascertain result to correct. Fit is 
        weighted by the error in the fit result of mean for each slice.
        '''
        means = []
        mean_errs = []
        centres = []
        for slc in self.slices:
            if not slc.veto: 
                #Avoids veto'd values
                means.append(slc.fit_vals['mean'])
                mean_errs.append(slc.fit_errs['mean'])
                centres.append(slc.centre)
        

        self.fit_vals, self.fit_cov = curve_fit(
                                        self.fit_model, 
                                        centres, 
                                        means, 
                                        sigma=mean_errs)
        #represents the actual fit errors (see curv_fit docs)
        self.fit_errs = np.sqrt(np.diag(self.fit_cov)) 
        
        ######  chi-sq:
        self.chisqdof = utl.calc_chisq(
                            f = self.fit_model, 
                            fit_params = self.fit_vals, 
                            x = centres, 
                            y_obs = means, 
                            sigma = mean_errs)

        
    def plot_all(self, setup):
        '''
        Method is called at the end of the process method, and plots all results
        from the automatic first-pass.
        '''
        self.plot_2D(setup)
        self.plot_means()            
        for slc in self.slices:
            self.plot_slice(slc)

    
    def fit_result_string(self, results, errors):
        '''
        Create a string containing fit results on plot for predefined models
        results if a user-defined model is being used.

        pm() is defined such that a string can be formed with appropriate signs.
        ie. y = 3m - 4; not y = 3m + -4.
        '''
        def pm(number):
            return '+' if number >= 0 else '-'

        if self.fit_model is utl.linear:
            fit_str = r"y = {:.1E}($\pm${:.1E})x {} {:.1E}($\pm${:.1E})".format(
                       results[0], 
                       errors[0],
                       pm(results[1]), 
                       abs(results[1]), 
                       errors[1])

        elif self.fit_model is utl.pol2:
            fit_str = (r"$y$ = {:.1E}($\pm${:.1E})$x^{}$ "
                       r"{} {:.1E}($\pm${:.1E})$x$ "
                       r"{} {:.1E}($\pm${:.1E})".format(
                       results[0], 
                       errors[0],
                       2,
                       pm(results[1]), 
                       abs(results[1]), 
                       errors[1],
                       pm(results[2]), 
                       abs(results[2]), 
                       errors[2])
            ) 

        elif self.fit_model is utl.pol3:
            fit_str = ("$y$ = {:.1E}($\pm${:.1E})$x^{}$\n"
                       "{} {:.1E}($\pm${:.1E})$x^{}$\n"
                       "{} {:.1E}($\pm${:.1E})$x$\n"
                       "{} {:.1E}($\pm${:.1E})".format(
                       results[0], 
                       errors[0],
                       3,
                       pm(results[1]), 
                       abs(results[1]), 
                       errors[1],
                       2,
                       pm(results[2]), 
                       abs(results[2]), 
                       errors[2],
                       pm(results[3]), 
                       abs(results[3]), 
                       errors[3])
            ) 
        
        elif self.fit_model is utl.pol4:
            fit_str = ("$y$ = {:.1E}($\pm${:.1E})$x^{}$\n"
                       "{} {:.1E}($\pm${:.1E})$x^{}$\n"
                       "{} {:.1E}($\pm${:.1E})$x^{}$\n"
                       "{} {:.1E}($\pm${:.1E})$x$\n"
                       "{} {:.1E}($\pm${:.1E})".format(
                       results[0], 
                       errors[0],
                       4,
                       pm(results[1]), 
                       abs(results[1]), 
                       errors[1],
                       3,
                       pm(results[2]), 
                       abs(results[2]), 
                       errors[2],
                       2,
                       pm(results[3]), 
                       abs(results[3]), 
                       errors[3],
                       pm(results[4]),
                       abs(results[4]),
                       errors[4])
            ) 

        else:
            strs = []
            for i,(p,e) in enumerate(zip(results,errors)):
                strs.append("p{}: {:.1E}$\pm${:.1E}\n".format(i, p, e))
            fit_str = "".join(strs)
                        

        return fit_str



    def plot_2D(self, setup):
        '''
        Plot the 2D histogram in the top left position in the figure, and 
        decorate with labels and lines showing the position of the slice limits 
        which have been applied. Title is added to the top to give context for 
        which parameter is being corrected.
        '''
        utl.plot2D_bh(self.output_ax[0], self.hist2D)
        #self.output_ax[0].set_xlabel(setup['v_label'])
        #self.output_ax[0].set_ylabel(setup['d_label'])
        self.output_ax[0].set_title(
            "{} vs. {}".format(setup['d_label'], setup['v_label'])
        )

        # decorate with lines over slices, and number slices
        for slc in self.slices:
            # 'axes[0].bin...' returns tuple of x-values of bin edges at bin 
            # position, and we take lower edge.
            mark_pos = self.hist2D.axes[0].bin(slc.lims[1])[0]
            y_low, y_high = (min(min(self.hist2D.axes[1])), 
                             max(max(self.hist2D.axes[1]))
            )

            self.output_ax[0].vlines(mark_pos, y_low, y_high, 
                 color='red', 
                 linestyles='dashed', 
                 linewidths=0.5, 
                 alpha=0.8
            )
            
            self.output_ax[0].text(
                slc.centre, 
                y_low+abs(0.1*y_low), 
                str(slc.N), 
                color='red', 
                ha='center'
            )
            
            
    def plot_means(self):
        '''
        Plots the distribution and of the fitted mean for each slice, with error
        bars, as well as the fit result. This is decorated with the result of
        the fit, as well as the chi-sq for the model.
        '''
        self.output_ax[1].clear()
        mean_centres = []
        means = []
        mean_errs = []
        for slc in self.slices:
            if not slc.veto: #Avoids veto'd values
                mean_centres.append(slc.centre)
                means.append(slc.fit_vals['mean'])
                mean_errs.append(slc.fit_errs['mean'])
            
        self.output_ax[1].errorbar(
            mean_centres, 
            means, 
            yerr = mean_errs, 
            marker='o', 
            ms=5, 
            ls='None', 
            capsize=1)

        total_xrange = (min(min(self.hist2D.axes[0])), 
            max(max(self.hist2D.axes[0]))
        )
        x_fitplot_vals = np.linspace(total_xrange[0], total_xrange[1], 300) 
        self.output_ax[1].plot(
            x_fitplot_vals, 
            self.fit_model(x_fitplot_vals, *self.fit_vals),
            label='fit', 
            c='red'
        )
        
        # if self.fit_vals[1] > 0:
        #     fit_str = r"y = {:.3}($\pm${:.3})x + {:.3}($\pm${:.3})".format(
        #         self.fit_vals[0], 
        #         self.fit_errs[0], 
        #         self.fit_vals[1], 
        #         self.fit_errs[1])
        # else:
        #     fit_str = r"y = {:.3}($\pm${:.3})x - {:.3}($\pm${:.3})".format(
        #         self.fit_vals[0], 
        #         self.fit_errs[0], 
        #         abs(self.fit_vals[1]), 
        #         self.fit_errs[1])

        fit_str = self.fit_result_string(self.fit_vals, self.fit_errs)

        chisq_str = r"$\chi^{}/dof$ = {:.3}/{}".format(
            2,
            self.chisqdof[0], 
            self.chisqdof[1])

        # Attempt to dodge the plotted data with the labels:
        #if self.fit_model is utl.linear or self.fit_model is utl.pol2:
        first_yval = self.slices[0].fit_vals['mean']
        last_yval = self.slices[-1].fit_vals['mean']
        if first_yval < last_yval:
            # if self.fit_vals[0] > 0: 
            rel_fit_x_pos = 0.05
            rel_fit_y_pos = 0.9
            fit_ha = 'left'
            rel_chi_x_pos = 0.95
            rel_chi_y_pos = 0.1
            chi_ha = 'right'

        else: 
            rel_fit_x_pos = 0.95
            rel_fit_y_pos = 0.9
            fit_ha = 'right'
            rel_chi_x_pos = 0.05
            rel_chi_y_pos = 0.1
            chi_ha = 'left'

        # else:
        #     rel_fit_x_pos = 0.95
        #     rel_fit_y_pos = 0.9
        #     fit_ha = 'right'
        #     rel_chi_x_pos = 0.05
        #     rel_chi_y_pos = 0.1
        #     chi_ha = 'left'

        self.output_ax[1].text(rel_fit_x_pos, rel_fit_y_pos, fit_str, 
            color='k', 
            ha=fit_ha, 
            va='top', 
            transform=self.output_ax[1].transAxes)
        self.output_ax[1].text(rel_chi_x_pos, rel_chi_y_pos, chisq_str, 
            color='k', 
            ha=chi_ha, 
            va='top', 
            transform=self.output_ax[1].transAxes)
    
    def plot_slice(self, slc):
        '''
        Plot a given slice of the distribution, as well as the fit result for 
        that slice.  Decorates plot with the slice number, and the reduced 
        chi-sq for the model.
        '''
        self.output_ax[slc.N+2].clear()
        fitvals = [slc.fit_vals['amp'],
                    slc.fit_vals['mean'],
                    slc.fit_vals['sigma'],
                    slc.fit_vals['const']
        ]
        if slc.fit_range is not None:
            x_low, x_high = slc.fit_range
        else:
            x_low,x_high = (min(min(slc.hist.axes[0])), 
                max(max(slc.hist.axes[0])))

        utl.plot1D_bh(self.output_ax[slc.N+2], slc.hist)

        x_fitplot_vals = np.linspace(x_low, x_high, 300)
        self.output_ax[slc.N+2].plot(
            x_fitplot_vals, 
            utl.gaussian_offset(x_fitplot_vals, *fitvals), 
            label='fit', 
            c='red')
        self.output_ax[slc.N+2].text(0.1, 0.75, str(slc.N), 
            color='red', 
            ha='center', 
            transform=self.output_ax[slc.N+2].transAxes)

        redchisq = slc.chisqdof[0] / slc.chisqdof[1]
        redchisq_str = r"$\chi^{}_{}$ = {:.3}".format(2, '{red}', redchisq)
        self.output_ax[slc.N+2].text(0.95, 0.75, redchisq_str, 
            color='k', 
            ha='right', 
            size='small', 
            transform=self.output_ax[slc.N+2].transAxes)

    
    def process(self, data, setup, cuts):
        '''
        Called inside __init__.  Serves to 'automate' initial processing of a
        given parameter when the class is initialised. 
        '''
        self.generate_2Dhisto(data, setup, cuts)
        self.create_slices()
        self.output_grid()  
        self.fit_means()
        self.plot_all(setup)
        
        
    def set_slice_fitparams(self, slc_idx, amp, mean, sigma, const):
        '''
        Allows the initial parameter values used in fitting of a slice to be 
        hard-coded, if the fit has failed. 
        Overall result is then reavaluated and the figure is updated.
        '''
        fitparams = {
            'amp': amp,
            'mean': mean,
            'sigma': sigma,
            'const': const
        }
        self.slices[slc_idx].hardset_fitestimates(fitparams)
        self.plot_slice(self.slices[slc_idx])
        self.fit_means()
        self.plot_means()
        
    
    def set_slice_fitrange(self, slc_idx, low, high):
        '''
        Allows the range of values used in fitting of a slice to be changed, 
        if the fit has failed, or the extremities of the distribution are
        biasing the result.
        Overall result is then reavaluated and the figure is updated.
        '''
        self.slices[slc_idx].fit_limits((low, high))
        self.plot_slice(self.slices[slc_idx])
        self.fit_means()
        self.plot_means()
        
        
        
    def rebin_slice(self, slc_idx, factor): 
        '''
        Allows the projected histogram for a given slice to be rebinned if low
        statistics are affecting the quality of the fit.  As the original 
        distribution is lost when rebinning, a deep-copy (duplicate) of the 
        original is first made, such that it can be restored if need be.
        Overall result is then reavaluated and the figure is updated.
        '''
        cache = copy.deepcopy(self.slices[slc_idx])
        self.cached_slices.append(cache)
        self.slices[slc_idx].rebin(factor)
        self.plot_slice(self.slices[slc_idx])
        self.fit_means()
        self.plot_means()
        
    def restore_slice(self, slc_idx):
        '''
        Restores a given slice to the original result of the automatic 
        processing.  
        For a rebinned slice this means checking the cache for the
        slice .by checking the slice number.
        For a "tweaked" fit (hard-coded init params or adjusted range) this 
        means removing the range limit and reassessing the initial parameters 
        via the algorithms provided.
        In each case, the overall result is then reavaluated and the figure is 
        updated.

        '''
        if self.slices[slc_idx].rebinned is True:
            for i,slc in enumerate(self.cached_slices):
                if slc.N == slc_idx:
                    self.slices[slc_idx] = copy.deepcopy(slc)
                    self.cached_slices.pop(i)
                    self.plot_slice(self.slices[slc_idx])
                    self.fit_means()
                    self.plot_means()
                # else:
                #     print("That slice has not been rebinned yet.  Nothing to restore")
        elif self.slices[slc_idx].tweaked_fit is True:
            self.slices[slc_idx].fit_range = None
            self.slices[slc_idx].tweaked_fit = False
            self.slices[slc_idx].fit_estimates()
            self.slices[slc_idx].fit_slice()
            self.plot_slice(self.slices[slc_idx])
            self.fit_means()
            self.plot_means()
        else:
            print("There are no changes to be undone in slice {}".format(slc_idx))
        
    def veto_slice(self, slc_idx):
        '''
        Allows a slice to be discounted from the final result.  This is simply 
        a boolean flag which the relevant methods check.
        Overall result is then reavaluated and the figure is updated.
        '''
        self.slices[slc_idx].veto = True
        self.output_ax[slc_idx+2].text(0.5, 0.5, 'VETOED', color='black', ha='center', weight='bold', 
            transform=self.output_ax[slc_idx+2].transAxes, 
            bbox=dict(facecolor='red', edgecolor='black', boxstyle='round, pad=0.5', alpha=0.7)
        )
        self.fit_means()
        self.plot_means()
        
    
    def unveto_slice(self, slc_idx):
        '''
        Allows a slice to be un-vetoed to be included in the by resetting the
        veto-flag.
        Overall result is then reavaluated and the figure is updated.
        '''
        self.slices[slc_idx].veto = False
        self.plot_slice(self.slices[slc_idx])
        self.fit_means()
        self.plot_means()

    def change_fit_model(self, func_str):
        '''
        Allows the model used to fit the slice-means to be changed.  
        Function will accept a string which matches the name of one of the 
        fitting functions defined in utils.py
        '''
        valid_models = {
            'linear': utl.linear, 
            'pol2': utl.pol2, 
            'pol3': utl.pol3,
            'pol4': utl.pol4
        }
        try: 
            self.fit_model = valid_models[func_str]
        except KeyError:
            print("Accepted models are: linear, pol2, pol3, pol4.")
            raise KeyError
        self.fit_means()
        self.plot_means()

    def user_def_model(self, function):
        '''
        Allows user to choose an arbitrary model with which to fit the means
        '''
        self.fit_model = function
        self.fit_means()
        self.plot_means()


    
class Slice:
    '''
    This class is used by the SimKinCorr class to store slice/projections of a
    2D histogram, and the data/parameters pertinent to that slice. Similar to 
    the SimKinCorr class, it is largely automated by the process_slices method.
    When a slice is initialised its initial fit parameters are algorithmically 
    estimated, and it is fitted with these parameters.

    There are three utility methods which are used by SimKinCorr to make
    adjustments to a given slice when fitting has failed:
        * fit_limits
        * hardset_fitestimates
        * rebin
    '''
    def __init__(self, N, proj, lims, centre):
        self.N = N # index identifying which slice in the array
        self.hist = proj # 1D histogram - projected slice
        self.lims = lims # bin range of each slice for projections
        self.centre = centre # bin centres for plotting
        self.fit_range = None # range withthin which to fit
        self.fit_vals = None # fit results
        self.fit_errs = None # fit errors
        self.fit_cov = None # fit cov
        self.fit_guesses = None # initial guesses for fit 
        self.chisqdof = None # tuple containing (X2, dof) from fit result
        self.veto = False # flag to ignore in final result 
        self.rebinned = False # flag for "restoring" histo
        self.tweaked_fit = False # flag for "restoring" histo

        self.process_slice()


    def fit_estimates(self):
        '''
        Generates initial guesses on parameters for Gaussian fits
          * amp - value of the bin with max counts.
          * const - arbirtarily the height of the last bin in range.
          * mean - x-value at position of max-bin
          * sig - calculated from crude estimate of FWHM
             -- FWHM estimated by starting at bin with max counts 
                and incrementing bin till half-value is found.
        Method adjusts for if a fit range has been hard-coded by the end-user.
        '''
        # if range has been specified, slice histogram       
        if self.fit_range is not None:
            #bh.loc converts data-val to bin-number
            hist = self.hist[
                bh.loc(self.fit_range[0]):bh.loc(self.fit_range[1])
            ]
        # otherwise reference normal self.hist
        else: 
            hist = self.hist
        
        amp_guess = np.max(hist.view())           
        const_guess = hist.view()[-1]
        
        #step over bins to find bin with half-max value FWHM
        n_bins = len(hist.axes[0])
        lims = (min(min(hist.axes[0])),max(max(hist.axes[0])))
        bin_len = (lims[1]-lims[0])/n_bins
               
        max_cnts = amp_guess
        max_idx = np.argmax(hist.view())
        halfmax_idx = np.int(max_idx / 2) # incase no half is reached below
        for j in np.arange(max_idx, n_bins, 1):
            counts = hist.view()[j]
            if counts <= max_cnts/2:
                halfmax_idx = j
                break

        mean_guess = lims[0] + max_idx*bin_len
        sig_guess = (abs(2*(max_idx-halfmax_idx) * bin_len)
                     / (2*np.sqrt(2*np.log(2))))

        self.fit_guesses = {
            'amp': amp_guess,
            'mean': mean_guess,
            'sigma': sig_guess,
            'const': const_guess
        }


    def fit_slice(self):
        #TODO IS SLICING HISTOGRAM REDUNDANT HERE??
        '''
        Fits slice distribution.  It is assumed that the distribution of a 
        'clean' data sample will be roughly Gaussian.  As we are solely 
        concerned with the mean position of the distribution, and not with the
        shape of the background, a Gaussian offset with a constant is used.
        Chi-sq of the model is also calculated.
        '''
        # if range has been specified, slice histogram      
        if self.fit_range is not None:
            #bh.loc converts data-val to bin-number
            hist = self.hist[
                bh.loc(self.fit_range[0]):bh.loc(self.fit_range[1])
            ]
        else: # otherwise reference normal self.hist
            hist = self.hist
        
        bin_centres = hist.axes[0].centers
        bin_vals = hist.view()
        #Error assumed to be Poissonian (counting) ie. sqrt(counts)
        bin_err = np.sqrt(hist.view())
        # set zero values (empty bins) to one to avoid zero-division in chi-sq
        bin_err[bin_err==0] = 1 
        p0=[
            self.fit_guesses['amp'], 
            self.fit_guesses['mean'], 
            self.fit_guesses['sigma'], 
            self.fit_guesses['const'] 
            ]
        popt, pcov = curve_fit(
            utl.gaussian_offset, 
            bin_centres, 
            bin_vals, 
            p0=p0, 
            sigma = bin_err)
            
        self.fit_vals = {
            'amp': popt[0],
            'mean': popt[1],
            'sigma': popt[2],
            'const': popt[3]
        }
        #represents the actual fit errors (see curv_fit docs)
        perr = np.sqrt(np.diag(pcov)) 
        self.fit_errs = {
            'amp': perr[0],
            'mean': perr[1],
            'sigma': perr[2],
            'const': perr[3]
        }
        self.fit_cov = pcov
        
        self.chisqdof = utl.calc_chisq(
            f = utl.gaussian_offset, 
            fit_params = popt, 
            x = bin_centres, 
            y_obs = bin_vals, 
            sigma = bin_err)
        
    
    def fit_limits(self, fit_lim):
        '''
        Set fit limits and reavaluate slice.
        '''
        self.tweaked_fit = True
        self.fit_range = fit_lim
        self.fit_estimates()
        self.fit_slice()
        
    
    def hardset_fitestimates(self, fit_params):
        '''
        Set initial fit parameters and reavaluate slice.
        '''
        self.tweaked_fit = True
        self.fit_guesses = fit_params
        self.fit_slice()
    
    def rebin(self, n):
        '''
        Rebin distribution by a factor n and reavaluate slice.
        '''
        self.rebinned = True
        self.hist = self.hist[::bh.rebin(n)]
        self.fit_estimates()
        self.fit_slice()
        

    def process_slice(self):
        '''
        Called inside __init__.  Serves to 'automate' evaluation of a slice. 
        '''
        self.fit_estimates()
        self.fit_slice()
        
