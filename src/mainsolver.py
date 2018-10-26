import matplotlib.pyplot as plt
import os
import re
from collections import defaultdict

import ms_measurments_extractor as ms
import modform_region_estimator as estimator
import data_simulator as simulate
import globalfuncs as gbfunc
import linear_constraints as lc
import statistics_and_plot as stats_plot


class MainSolver(object):
    
    def __init__(self, specifications):
        
        self._specifications               = specifications
        
        self._ms_measurments_extractor = ms.MSMeasurementsExtractor()
        self._simulated_data           = simulate.DataSimulator(specifications)
        
        self._tex_output_file          = gbfunc.TexFileGenerator(specifications.tex_output_file, specifications.output_dir)
        
    
    def run_simulated_mode(self, protein_seq, D_sites_modtypes):
        
        """Instance method
        
        This function performs simulation, compute bounds for a given peak area measurment noise and number of samples and calls relevant print output functions.
        
        Parameters
        ----------
        protein_seq : str
          Protein amino-acid sequence string.
        
        D_sites_modtypes : defaultdict
          Python dictionary having aa residues mapped to list of modifications on respective residues.
        
        Notes
        -----
        
        """
        # use try-catch to ensure appropriate input is given and also linear constraint object has sought structure
        
        # Dictionary to store modform variables and their simulated abundances
        D_vars_simulated_abundances = {}
        
        # populating the dictionary via a data structure read from an input file
        self._simulated_data.D_sites_modtypes = D_sites_modtypes
                        
        self._simulated_data.generate_modform_distribution()
        
        L_D_lower_bounds = []; L_D_upper_bounds = []
        
        #linear_constraints = lc.LinearConstraints(D_sites_modtypes=D_sites_modtypes, specifications=self._specifications)
        
        sample_size = int(self._specifications.sample_size)
        for i in range(sample_size):
            
            if self._specifications.file_tdms_data:
                self._simulated_data.generate_TDMS_data()
            if self._specifications.file_bums_mrm_data:
                self._simulated_data.generate_BUMS_MRM_data()
            
            for var,mform in self._simulated_data.modform_space.D_variables_modforms.iteritems():
                D_vars_simulated_abundances[var] = mform.percentage_abundance
            #end for var,mform in D_vars_modforms.iteritems()
            
            if self._specifications.file_tdms_data:
                self._ms_measurments_extractor.get_list_of_isomeric_modforms_peak_areas(self._specifications.file_tdms_data)
                        
            if self._specifications.file_bums_mrm_data:
                self._ms_measurments_extractor.get_dictionary_of_peptides_and_isomeric_peak_areas(self._specifications.file_bums_mrm_data)
            
            linear_constraints = lc.LinearConstraints(D_sites_modtypes=D_sites_modtypes, specifications=self._specifications)
            
            linear_constraints.derive_linear_constraints_from_MS_measurments(self._ms_measurments_extractor.L_isomeric_modforms_peak_areas, self._ms_measurments_extractor.D_peptide_isomeric_peak_areas, D_variables_modforms_simulated=self._simulated_data.modform_space.D_variables_modforms, D_PuLP_variables_simulated=self._simulated_data.modform_space.D_PuLP_variables)
            # Verifying whether linear constraints object is populated
            if len(linear_constraints.L_linear_constraints)==0:
                raise Exception("Linear constraint list is empty")
            
            
            if self._specifications.variables_grouping and self._specifications.file_bums_mrm_data:
                L_grouped_lists = linear_constraints.get_grouped_variables()
                linear_constraints.get_linear_constraints_with_grouped_variables(L_grouped_lists)
            
            modform_region_estimator = estimator.ModformRegionEstimator(linear_constraints)
            
            (D_lower_bounds, D_upper_bounds) = modform_region_estimator.compute_bounds()
                        
            L_D_lower_bounds.append(D_lower_bounds)
            L_D_upper_bounds.append(D_upper_bounds)
            
        #end for i in range(sample_size)
        
        stats = stats_plot.StatisticsAndPlot()
        
        D_mean_lower_bounds, D_stdev_lower_bounds = stats.get_dictionary_of_mean_and_standard_deviation_of_bounds(L_D_lower_bounds) 
        D_mean_upper_bounds, D_stdev_upper_bounds = stats.get_dictionary_of_mean_and_standard_deviation_of_bounds(L_D_upper_bounds)
        
        print "\n####### Grouped variables ###########"
        for k,v in sorted(linear_constraints.D_grouped_variables.iteritems(), key=lambda x:gbfunc.natural_keys(x[0])):
            print "%s: %s"%(k,v)
        #end for k,v in sorted(linear_constraints.D_grouped_variables.iteritems(), key=lambda x:gbfunc.natural_keys(x[0]))
        
        print "\n####### Bounds on percentage abundance #########"
        for k,v in sorted(D_lower_bounds.iteritems(), key=lambda x:gbfunc.natural_keys(x[0])):
            print "%s <= %s <= %s"%(D_lower_bounds[k],k,D_upper_bounds[k])
        #end for
        
        (mean, stdev) = stats.get_statistics_on_bounds(D_lower_bounds, D_upper_bounds)
        print "mean and stdev of bounds: %.3f, %.3f"%(mean, stdev)

        for var,gp in linear_constraints.D_grouped_variables.iteritems():
            L_v = [v.strip() for v in re.split(r'[+-]', gp)]
            t_abundance = sum([float(D_vars_simulated_abundances[k]) for k in L_v])

            for vv in L_v:
                D_vars_simulated_abundances.pop(vv,None)
            #end for
            
            D_vars_simulated_abundances[var] = t_abundance
        
        #end for var,gp in linear_constraints.D_grouped_variables.iteritems()
        mform_distribution_plot = stats.plot_simulated_values(D_vars_simulated_abundances, Output_dir=self._specifications.output_dir)                
        mform_region_plot_pdf = stats.plot_bars_for_simulated_data(D_mean_lower_bounds, D_mean_upper_bounds, D_stdev_lower_bounds, D_stdev_upper_bounds, D_variables_simulated_values=D_vars_simulated_abundances, Output_dir=self._specifications.output_dir)

        stats.get_sensitivity_performance_statistics_over_all_samples(L_D_lower_bounds, L_D_upper_bounds, D_vars_simulated_abundances, Output_dir=self._specifications.output_dir)
        
        # printing output files in an Output folder
        # print file having modform variables(a_i) and modforms
        # print file with grouped variables and abundance bounds
        
        dirname = self._specifications.output_dir
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        #end if

        self._tex_output_file.print_content_in_tex_file(protein_seq, D_sites_modtypes, D_mean_lower_bounds, D_mean_upper_bounds, mform_region_plot_pdf, linear_constraints.modform_space.D_variables_modforms, linear_constraints.D_grouped_variables)
        
        
        fstring = 'Output/modforms.txt'
        fout_modforms = open(fstring,'w')
        print_string = "#Modform variables and their modifications map\n"
        for mvar,mform in sorted(linear_constraints.modform_space.D_variables_modforms.iteritems(), key=lambda x: gbfunc.natural_keys(x[0])):
            mform_string = '%s\n'%(mform)
            print_string += mform_string
                
        fout_modforms.write(print_string)
        print "\n[File created] %s -- Contains Modform variables and their respective modifications map"%fstring
        fout_modforms.close()
        
        fstring = 'Output/percentage_abundance_bounds.txt'
        fout_abundances = open(fstring,'w')
        print_string = '#Percentage abundance bounds\n\n'
        print_string += '###### Grouped variables ######\n'
        
        for gv,gp in sorted(linear_constraints.D_grouped_variables.iteritems(), key=lambda x:gbfunc.natural_keys(x[0])):
            print_string += '%s: %s\n'%(gv,gp)
        #end for gv,gp in sorted(linear_constraints.D_grouped_variables.iteritems(), key=lambda x:gbfunc.natural_keys(x[0]))
        
        print_string += "\n\n######### Percentange abundance bounds ###########\n"
        for k,v in sorted(D_lower_bounds.iteritems(), key=lambda x:gbfunc.natural_keys(x[0])):
            print_string += "%s <= %s <= %s\n"%(D_lower_bounds[k],k,D_upper_bounds[k])
        #end for k,v in sorted(D_lower_bounds.iteritems(), key=lambda x:gbfunc.natural_keys(x[0]))
        
        print_string += "\n\n###### Statistics on bounds##########\n"
        print_string += "Mean: %.3f\n"%round(mean,3)
        print_string += "Standard deviation: %.3f\n"%round(stdev,3)
        fout_abundances.write(print_string)
        print "\n[File created] %s -- Contains Bounds on percentage abundance for modforms"%fstring
        fout_abundances.close()
        
        
        
    
    def run_real_data_mode(self, protein_seq, D_sites_modtypes):
        
        """Instance method
        This function compute bounds for real data and calls relevant print output functions.
        
        Parameters
        ----------
        protein_seq : str
          Amino-acid sequence of the protein under investigation
        
        D_sites_modtypes : defaultdict(list)
          Python dictionary mapping aa residues to python list of modification types on the residues.
        """
                
        # get TDMS isomeric modforms peak areas and isomeric peptides peak areas from MS measurments
        self._ms_measurments_extractor.get_list_of_isomeric_modforms_peak_areas(self._specifications.file_tdms_data)
        self._ms_measurments_extractor.get_dictionary_of_peptides_and_isomeric_peak_areas(self._specifications.file_bums_mrm_data)

        linear_constraints = lc.LinearConstraints(D_sites_modtypes=D_sites_modtypes, specifications=self._specifications)
        linear_constraints.derive_linear_constraints_from_MS_measurments(self._ms_measurments_extractor.L_isomeric_modforms_peak_areas, self._ms_measurments_extractor.D_peptide_isomeric_peak_areas)

        # Verifying whether linear constraints object is populated
        if len(linear_constraints.L_linear_constraints)==0:
            raise Exception("Linear constraint list is empty")

        if self._specifications.variables_grouping and self._specifications.file_bums_mrm_data:
            L_grouped_lists = linear_constraints.get_grouped_variables()
            linear_constraints.get_linear_constraints_with_grouped_variables(L_grouped_lists)
        
        modform_region_estimator = estimator.ModformRegionEstimator(linear_constraints)
        
        (D_lower_bounds, D_upper_bounds) = modform_region_estimator.compute_bounds()
        
        print "\n####### Grouped variables ###########"
        for k,v in sorted(linear_constraints.D_grouped_variables.iteritems(), key=lambda x:gbfunc.natural_keys(x[0])):
            print "%s: %s"%(k,v)
        #end for k,v in sorted(linear_constraints.D_grouped_variables.iteritems(), key=lambda x:gbfunc.natural_keys(x[0]))
        
        print "\n####### Bounds on percentage abundance #########"
        for k,v in sorted(D_lower_bounds.iteritems(), key=lambda x:gbfunc.natural_keys(x[0])):
            print "%s <= %s <= %s"%(D_lower_bounds[k],k,D_upper_bounds[k])
        #end for

        stats = stats_plot.StatisticsAndPlot()
        
        (mean, stdev) = stats.get_statistics_on_bounds(D_lower_bounds, D_upper_bounds)
        print "mean and stdev of bounds: %.3f, %.3f"%(mean, stdev)

        mform_region_plot_pdf = stats.plot_bars(D_lower_bounds, D_upper_bounds, Output_dir=self._specifications.output_dir)
        
        dirname = self._specifications.output_dir
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        #end if

        self._tex_output_file.print_content_in_tex_file(protein_seq, D_sites_modtypes, D_lower_bounds, D_upper_bounds, mform_region_plot_pdf, linear_constraints.modform_space.D_variables_modforms, linear_constraints.D_grouped_variables)
        
        fstring = 'Output/modforms.txt'
        fout_modforms = open(fstring,'w')
        print_string = "#Modform variables and their modifications map\n"
        
        
        for mvar,mform in sorted(linear_constraints.modform_space.D_variables_modforms.iteritems(), key=lambda x: gbfunc.natural_keys(x[0])):
            mform_string = '%s\n'%(mform)
            print_string += mform_string
                        
        #end for mvar,mform in sorted(linear_constraints.D_variables_modforms.iteritems(), key=lambda x: gbfunc.natural_keys(x[0]))
        
        fout_modforms.write(print_string)
        print "\n[File created] %s -- Contains Modform variables and their respective modifications map"%fstring
        fout_modforms.close()
        
        fstring = 'Output/percentage_abundance_bounds.txt'
        fout_abundances = open(fstring,'w')
        print_string = '#Percentage abundance bounds\n\n'
        print_string += '###### Grouped variables ######\n'
        
        for gv,gp in sorted(linear_constraints.D_grouped_variables.iteritems(), key=lambda x:gbfunc.natural_keys(x[0])):
            print_string += '%s: %s\n'%(gv,gp)
        #end for gv,gp in sorted(linear_constraints.D_grouped_variables.iteritems(), key=lambda x:gbfunc.natural_keys(x[0]))
        
        print_string += "\n\n######### Percentange abundance bounds ###########\n"
        for k,v in sorted(D_lower_bounds.iteritems(), key=lambda x:gbfunc.natural_keys(x[0])):
            print_string += "%s <= %s <= %s\n"%(D_lower_bounds[k],k,D_upper_bounds[k])
        #end for k,v in sorted(D_lower_bounds.iteritems(), key=lambda x:gbfunc.natural_keys(x[0]))
        
        print_string += "\n\n###### Statistics on bounds##########\n"
        print_string += "Mean: %.3f\n"%round(mean,3)
        print_string += "Standard deviation: %.3f\n"%round(stdev,3)
        fout_abundances.write(print_string)
        print "\n[File created] %s -- Contains Bounds on percentage abundance for modforms"%fstring
        fout_abundances.close()
        
                
        
    def run(self):
        """Instance method
        
        This function sets the program running either in simulation mode or real data mode.
        """
        
        # reading amino acid sequence from the fasta file
        protein_seq = gbfunc.get_fasta_seq_from_the_file(self._specifications.file_fasta_sequence)
        
        # reading and storing the apriori information - sites and modifications
        D_sites_modtypes = gbfunc.read_sites_and_modification_types_from_file(self._specifications.file_modifications)
                
        if (self._specifications.data_mode=="") or (not self._specifications.data_mode=="simulation" and not self._specifications.data_mode=="real"):
            raise Exception("Select the data mode in the specification file -- 'real' OR 'simulation'")
                
        # populating linear constraints object
        if self._specifications.data_mode.lower()=="simulation":
            
            self.run_simulated_mode(protein_seq, D_sites_modtypes)
            
        elif self._specifications.data_mode.lower()=="real":
            
            self.run_real_data_mode(protein_seq, D_sites_modtypes)
                
                   
        #end if-elif self._specifications.data_mode..
                 
        
        
