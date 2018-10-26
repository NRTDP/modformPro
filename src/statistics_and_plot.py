import os
import numpy as np
import matplotlib.pyplot as plt

import globalfuncs as gbfunc

class StatisticsAndPlot(object):
    
    """This class provides methods to compute statistics and for plotting.
    """
    
    
    def plot_bars(self, D_lower_bounds, D_upper_bounds, Output_dir=None):
        """Instance method
        
        Plot bars showing lower and upper bounds of percentage concentration.
        
        Parameters
        ----------
        D_lower_bounds : dict
          Python dictionary containing modform :class:`.pulp` variables mapped to their lower bounds.
        
        D_upper_bounds : dict
          Python dictionary containing modform :class:`.pulp` variables mapped to their upper bounds.
                        
        Output_dir : str
          Output directory name. Default is ``Output``.
        """
        L_lower_bounds = []
        L_bounds_range = []
        
        for k in sorted(D_lower_bounds.keys(), key=lambda x: gbfunc.natural_keys(x)):
            #if v==0: continue

            L_lower_bounds.append(float(D_lower_bounds[k]))
            r = float(D_upper_bounds[k]) - float(D_lower_bounds[k])
            L_bounds_range.append(r)
        #end for k,v in D_lower_bounds

        N = len(D_lower_bounds)
        ind = [2*i+2 for i in np.arange(N)]
        
        width = 1.5
        
        plt.figure(4)
        p1 = plt.bar(ind, L_lower_bounds, width, color='white', edgecolor='white',align='center')
        p2 = plt.bar(ind, L_bounds_range, width, bottom=L_lower_bounds, facecolor='lightgray', edgecolor='gray',align='center', label="% conc. bounds")
        plt.ylim(-0.5,50)
        plt.xlim(0,2*N+4)
        plt.xlabel('Modform concentration variable')
        plt.ylabel('Bounds of percentage concentration')
        #plt.title('modform region calculation for simulated data')
        plt.xticks(ind, sorted(D_lower_bounds.keys(), key=lambda x:gbfunc.natural_keys(x)), rotation=90)
        plt.legend()
        #fout = "mod-form-region-bounds-" + options.output_file_suffix + ".pdf"

        # create output directory if it doesn't exist
        dirname = Output_dir
        if not os.path.exists(dirname):
            os.makedirs(dirname)

        fplot_string = 'modform-region-bounds.pdf'
        fout = "%s/%s"%(Output_dir,fplot_string)
        print "\n[Bar plot generated] %s -- Bars representing bounds on percentage abundance for modforms: "%fout
        plt.savefig(fout)

        return fplot_string
    
    def get_statistics_on_bounds(self, D_lower_bounds, D_upper_bounds):
        
        L_bounds_ranges = []
        for k,v in D_lower_bounds.iteritems():
            diff = float(D_upper_bounds[k]) - float(v)
            L_bounds_ranges.append(diff)
        #end for k,v in D_lower_bounds.iteritems()
        
        mean  = np.mean(L_bounds_ranges)
        stdev = np.std(L_bounds_ranges)

        return (mean,stdev)

    def get_sensitivity_statistics_for_one_sample(self, D_lower_bounds, D_upper_bounds, D_variables_simulated_values):
        
        True_hits= 0; False_hits=0; L_False_hits_variables = []
        for k,v in D_variables_simulated_values.iteritems():
            
            lower_bound = D_lower_bounds[k]
            upper_bound = D_upper_bounds[k]
            
            if float(lower_bound) <= float(v) and float(v) <= float(upper_bound):
                True_hits += 1
            else:
                False_hits += 1
                L_False_hits_variables.append(k)
            #end if
        #end for
        
        return (True_hits, False_hits, L_False_hits_variables)
        

    def get_sensitivity_performance_statistics_over_all_samples(self, L_D_lower_bounds, L_D_upper_bounds, D_vars_simulated_abundances, Output_dir=None):
        '''
        This function provides statistics on true hits for simulated data
        '''
        # List of sensitivity values over all samples
        L_sensitivity_values = []
        
        # Dictionary for variables mapped to number of false hits
        D_variables_False_hits = {}
        for key in L_D_lower_bounds[0].keys():
            D_variables_False_hits[key]=0
        
        
        for i in range(len(L_D_lower_bounds)):
            True_hits, False_hits, L_False_hits_variables = self.get_sensitivity_statistics_for_one_sample(L_D_lower_bounds[i], L_D_upper_bounds[i], D_vars_simulated_abundances)
            Total_hits = int(True_hits)+int(False_hits)
            sensitivity = int(True_hits)/(1.0*Total_hits)
            L_sensitivity_values.append(sensitivity)
            
            for v in L_False_hits_variables:
                if v in D_variables_False_hits.keys():
                    new_value = int(D_variables_False_hits[v]) + 1
                    D_variables_False_hits[v] = new_value
                else:
                    D_variables_False_hits[v] = 1

                #end if-else

            #end for v in L_False_hits_variables
            
        #end for i in range(len(L_D_lower_bounds))

        # statistics on overall Sensitivity
        mean_sensitivity  = "%.3f"%round(np.mean(L_sensitivity_values),3)
        stdev_sensitivity = "%.3f"%round(np.std(L_sensitivity_values),3)

        print "mean_sensitivity: ", mean_sensitivity
        print "stdev sensitivity: ", stdev_sensitivity
        
        sample_size = len(L_D_lower_bounds)

        L_variable_sensitivity_values = []
        L_variables = []
        for k,v in sorted(D_variables_False_hits.iteritems(), key=lambda x: gbfunc.natural_keys(x[0])):
            var_sensitivity = (sample_size - int(v))/(1.0*sample_size)
            L_variable_sensitivity_values.append(var_sensitivity)
            L_variables.append(k)
        #end for
        
        N = len(L_variables)
        ind = np.arange(N)
        width = 0.3
        
        plt.figure(2)
        p = plt.bar(ind, L_variable_sensitivity_values, width, facecolor='red',align='center')
        
        plt.xlabel('Modform concentration variable')
        plt.ylabel('Hit rate for a variable')
        plt.title('Variable wise sensitivity over all samples')
        plt.xticks(ind, L_variables, rotation=90)
        #plt.show()
        
        # create output directory if it doesn't exist
        dirname = Output_dir
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        fout="%s/variable-wise-sensitivity-over-all-samples-simulated-data.pdf"%Output_dir
        plt.savefig(fout)
        
        return (fout, mean_sensitivity, stdev_sensitivity)
    
    
    def get_dictionary_of_mean_and_standard_deviation_of_bounds(self, list_dictionary_variables_bounds):
        '''
        This function returns mean and standard deviation of lower bounds OR upper bounds when measurments are done multiple times. For simulated data, this amounts to sample multiple noise values.
        '''
        D_variables_bounds_mean  = {}
        D_variables_bounds_stdev = {}
        
        for key in list_dictionary_variables_bounds[0].keys():
            L_bounds_key = [float(d[key]) for d in list_dictionary_variables_bounds]
            mean_bounds  = np.mean(L_bounds_key)
            stdev_bounds = np.std(L_bounds_key)
            
            D_variables_bounds_mean[key]  = "%.3f"%round(mean_bounds,3)
            D_variables_bounds_stdev[key] = "%.3f"%round(stdev_bounds,3)
            
        #end for key in list_dictionary_variables_bounds[0].keys()
        
        return (D_variables_bounds_mean, D_variables_bounds_stdev)


    def plot_bars_for_simulated_data(self, D_mean_lower_bounds, D_mean_upper_bounds, D_stdev_lower_bounds, D_stdev_upper_bounds, D_variables_simulated_values=None, Output_dir=None):
                
        '''
        Ref: https://matplotlib.org/api/pyplot_api.html
        
        # a stacked bar plot with errorbars
        N = 5
        menMeans = (20, 35, 30, 35, 27)
        womenMeans = (25, 32, 34, 20, 25)
        menStd = (2, 3, 4, 1, 2)
        womenStd = (3, 5, 2, 3, 3)
        ind = np.arange(N) # the x locations for the groups
        width = 0.35       # the width of the bars: can also be len(x) sequence

        p1 = plt.bar(ind, menMeans, width, color='#d62728', yerr=menStd)
        p2 = plt.bar(ind, womenMeans, width, bottom=menMeans, yerr=womenStd)
        
        plt.ylabel('Scores')
        plt.title('Scores by group and gender')
        plt.xticks(ind, ('G1', 'G2', 'G3', 'G4', 'G5'))
        plt.yticks(np.arange(0, 81, 10))
        plt.legend((p1[0], p2[0]), ('Men', 'Women'))
        
        plt.show()
        '''
        
        L_mean_lower_bounds  = []
        L_bounds_range       = []
        L_stdev_lower_bounds = []
        L_stdev_upper_bounds = []
        L_simulated_values   = []
        L_simulated_top      = []
        for k in sorted(D_mean_lower_bounds.keys(), key=lambda x: gbfunc.natural_keys(x)):
            L_mean_lower_bounds.append(float(D_mean_lower_bounds[k]))
            L_stdev_lower_bounds.append(float(D_stdev_lower_bounds[k]))
            L_stdev_upper_bounds.append(float(D_stdev_upper_bounds[k]))
            r = float(D_mean_upper_bounds[k]) - float(D_mean_lower_bounds[k])
            L_bounds_range.append(r)
            L_simulated_values.append(float(D_variables_simulated_values[k]))
            L_simulated_top.append(0)
        #end for k in sorted(D_mean_lower_bounds.keys(), key=gbfunc.natural_keys)
        
        N = len(D_mean_lower_bounds)
        ind = [2*i+2 for i in np.arange(N)]
        
        width = 1.5
        
        plt.figure(1)
        
        p1 = plt.bar(ind, L_mean_lower_bounds, width, yerr=L_stdev_lower_bounds, color='white', edgecolor='white',align='center', error_kw={'ecolor':'black'})
        p2 = plt.bar(ind, L_bounds_range, width, bottom=L_mean_lower_bounds,yerr=L_stdev_upper_bounds, facecolor='lightgray', edgecolor='gray',align='center', error_kw={'ecolor':'black'}, label='% conc. bounds')

        #plot for original simulated values
        p3 = plt.bar(ind, L_simulated_top, width, linewidth=1.5, bottom=L_simulated_values, color='red', edgecolor='red', align='center', label='simulated % conc.')
        plt.ylim(-0.5,50)
        plt.xlim(0,2*N+4)
        plt.xlabel('Modform concentration variable')
        plt.ylabel('Bounds of percentage concentration')
        #plt.title('modform region calculation for simulated data')
        plt.xticks(ind, sorted(D_mean_lower_bounds.keys(), key=gbfunc.natural_keys), rotation=90)
        
        plt.legend()
        #plt.show()
        
        # create output directory if it doesn't exist
        dirname = Output_dir
        if not os.path.exists(dirname):
            os.makedirs(dirname)
            
        fplot_string = 'modform-region-bounds-simulated-data.pdf'
        fout="%s/%s"%(Output_dir,fplot_string)
        plt.savefig(fout)
    
        return fplot_string
    


    def plot_simulated_values(self, D_variables_simulated_values, Output_dir=None):
        
        '''
        Ref: https://matplotlib.org/api/pyplot_api.html
        
        # a stacked bar plot with errorbars
        N = 5
        menMeans = (20, 35, 30, 35, 27)
        womenMeans = (25, 32, 34, 20, 25)
        menStd = (2, 3, 4, 1, 2)
        womenStd = (3, 5, 2, 3, 3)
        ind = np.arange(N) # the x locations for the groups
        width = 0.35       # the width of the bars: can also be len(x) sequence

        p1 = plt.bar(ind, menMeans, width, color='#d62728', yerr=menStd)
        p2 = plt.bar(ind, womenMeans, width, bottom=menMeans, yerr=womenStd)
        
        plt.ylabel('Scores')
        plt.title('Scores by group and gender')
        plt.xticks(ind, ('G1', 'G2', 'G3', 'G4', 'G5'))
        plt.yticks(np.arange(0, 81, 10))
        plt.legend((p1[0], p2[0]), ('Men', 'Women'))
        
        plt.show()
        '''
        L_simulated_values   = []
        L_simulated_top      = []
        for k in sorted(D_variables_simulated_values.keys(), key=gbfunc.natural_keys):
            
            L_simulated_values.append(float(D_variables_simulated_values[k]))
            L_simulated_top.append(0)
        #end for k in sorted(D_mean_lower_bounds.keys(), key=gbfunc.natural_keys)
        
        N = len(D_variables_simulated_values)
        ind = [2*i+4 for i in np.arange(N)]
        
        width = 1.5
        
        plt.figure(3)
        
        p1 = plt.bar(ind, L_simulated_top, width, linewidth=2, bottom=L_simulated_values, facecolor='red', align='center')
        plt.ylim(-0.5,50)
        plt.xlim(0,2*N+4)
        
        plt.xlabel('Modform concentration variable')
        plt.ylabel('percentage concentration')
        #plt.title('modform region calculation for simulated data')
        plt.xticks(ind, sorted(D_variables_simulated_values.keys(), key=gbfunc.natural_keys), rotation=90)
        #plt.show()
        
        # create output directory if it doesn't exist
        dirname = Output_dir
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        
        fplot_string = 'modform-distribution-simulated-data.pdf'
        fout="%s/%s"%(Output_dir,fplot_string)
        plt.savefig(fout)
    
        return fplot_string
