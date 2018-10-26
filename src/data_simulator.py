from collections import defaultdict

import random
import globalfuncs   as gbfunc
import modform_space as mform_space

class DataSimulator(object):
    
    """
    This class simulates modform distribution and then simulates the MS process on the simulated modform distribution. 
    
    Based on fields set by the user in the specification file, peak areas for BUMS/MRM or TDMS or both are calculated which sets the stage for estimation of modform region from the simulated data.
    
    
    Parameters
    ----------
    specifications : :class:`.ModformProSpecifications`
      - ModformSpecifications object
      - Contains information about the options set by the user in the specification file.

    D_sites_modtypes : dict
      - Python dictionary containing aa residues mapped to the list of possible modifications.
      - In the current implementation, this is set after the creation of the instance object of this class.
    """
    def __init__(self, specifications=None, D_sites_modtypes=None):
        
        self.specifications    = specifications
        self._D_sites_modtypes = D_sites_modtypes
        
        if specifications:
            self._modform_space    = mform_space.ModformSpace(specifications, specifications.modes)
        
    
    @property
    def D_sites_modtypes(self):
        """Instance attribute
        
        * Python dictionary having aa residues mapped to the list of possible modifications
        * E.g., ``D = {'T32': ['P'], 'S55': ['Ac', 'P']}``
        """
        return self._D_sites_modtypes

    @D_sites_modtypes.setter
    def D_sites_modtypes(self, D_sites_modtypes):
        
        self._D_sites_modtypes = D_sites_modtypes
        
    
    @property
    def modform_space(self):
        """Instance attribute
        
        Instance of a :class:`.ModformSpace` object
        """
        return self._modform_space

    @modform_space.setter
    def modform_space(self, modform_space):
        self._modform_space = modform_space

    
    
    def generate_modform_distribution(self):
        """Instance method
         
           This method generates simulated modform distribution and populates the instance attribute ``modform_space``.
        """
        
        #asserting that nmodes(number of protein isomeric modes) must be less than nsites (number of sites)
        nsites = len(self._D_sites_modtypes)
        if int(self.specifications.number_modes) > int(nsites)+1:
            raise Exception("Number of protein isomeric modes simulated should be less than number of modification sites")
            
        
        if len(self.specifications.modes)!=int(self.specifications.number_modes):
            raise Exception("Number of modes must match modes provided")
        # generating modform space
        self._modform_space.generate_modform_space()
        
                        
        #specifying the percentage abundance of modforms
        total_molecules = sum([float(n) for n in self._modform_space.D_variables_nmolecules.values()])
        
        # setting the percentage abundance (two decimal places) for each modform
        for var, mform in self._modform_space.D_variables_modforms.iteritems():
            nmolecules = self._modform_space.D_variables_nmolecules[var]
            percentage_abundance = 100*float(nmolecules)/float(total_molecules)

            # storing float value upto three decimal places
            percentage_abundance_three_decimal_places = "%.3f"%round(percentage_abundance,3)
            mform.percentage_abundance = percentage_abundance_three_decimal_places
                    
        #end for var, mform in self._modform_space..

        if sorted(self._modform_space.D_variables_modforms.keys())!=sorted(self._modform_space.D_PuLP_variables.keys()):
                  raise Exception("Keys in the both the dictionaries must match")
                        
        print "\n#######simulated modform distribution percentage abundances#########"
        for k, v in sorted(self._modform_space.D_variables_modforms.iteritems(), key=lambda x: gbfunc.natural_keys(x[0])):
            print "%s: %s"%(v,v.percentage_abundance)
        
        
        
    def get_list_isomeric_modforms_peak_areas_lists(self, D_variables_modforms):
        """Instance method
        
        Its a subordinate method, i.e., called within the class only by a member method.
        

        Parameters
        ----------
         D_variables_modforms : dict
           Python dictionary of modform variables mapped to their modform objects.
        
        Returns
        -------
        list

        Examples
        --------
        >>> from data_simulator import DataSimulator
        >>> from modform import Modform
        >>>
        >>> D_variables_modforms = {}
        >>> D_variables_modforms['a1'] = Modform(var_name='a1', D_protein_sites_modifications={'0': 'UNMOD'})
        >>> D_variables_modforms['a2'] = Modform(var_name='a2', D_protein_sites_modifications={'S20': 'P'})
        >>> D_variables_modforms['a3'] = Modform(var_name='a3', D_protein_sites_modifications={'S15': 'P'})
        >>> D_variables_modforms['a4'] = Modform(var_name='a4', D_protein_sites_modifications={'S15': 'P', 'S20': 'P'})
        >>>
        >>> for k,v in D_variables_modforms.iteritems():
        ...     if k=='a1': v.percentage_abundance='0.000'
        ...     elif k=='a2': v.percentage_abundance='11.111'
        ...     elif k=='a3': v.percentage_abundance='44.444'
        ...     elif k=='a4': v.percentage_abundance='44.444'
        >>>
        >>> D_sites_modtypes = {'S20': ['P'], 'S15': ['P']}
        >>> sim_data = DataSimulator(D_sites_modtypes=D_sites_modtypes)
        >>> L_peptides = [(1, 5), (6, 10), (11, 20), (21, 32)]
        >>>
        >>> D = sim_data.get_dictionary_peptides_isomeric_peak_areas(D_variables_modforms, L_peptides)
        >>> print D
        defaultdict(<type 'list'>, {'11-20': [[{'UNMOD': '0'}, '0.000'], [{'P': 1}, '55.555'], [{'P': 2}, '44.444']]})
        """
        
        L_isomeric_modforms_peak_areas = []
        
        for k, v in D_variables_modforms.iteritems():

            # getting a modification frequency map for this modform
            D_map = v.get_modifications_frequency_map()

            
            percentage_peak_area_contribution = v.percentage_abundance

            if D_map in [t[0] for t in L_isomeric_modforms_peak_areas]:
                t_index =  [t[0] for t in L_isomeric_modforms_peak_areas].index(D_map)

                p_area = L_isomeric_modforms_peak_areas[t_index][1]

                percentage_peak_area = float(percentage_peak_area_contribution) + float(p_area)
                L_isomeric_modforms_peak_areas[t_index][1] = str(percentage_peak_area)

            else:
                L_isomeric_modforms_peak_areas.append([D_map, str(percentage_peak_area_contribution)])
        
        return L_isomeric_modforms_peak_areas

    
    def generate_TDMS_data(self):
        """Instance method
        
        * This method computes peak areas for Intact mass spectrometry for simulated modform distribution.
        * This is done only if the option ``FILE_TDMS_DATA`` is specified in the specification file. The content is then printed in the file specified.
        * Noise can also be introduced in the peak area values.
        
        
        Notes
        -----
        Prints a file having information of the isomeric modforms and their respective percentage peak areas.
        
        """
        D_variables_modforms = self._modform_space.D_variables_modforms
        L = self.get_list_isomeric_modforms_peak_areas_lists(D_variables_modforms)
        
        L_sorted = sorted(L, key=lambda x: sum([int(a) for a in x[0].values()]))        
        
        #adding measurment error in the simulated relative peak area values
        L_relative_peak_areas = [a[1] for a in L_sorted]
        for index, rarea in enumerate(L_relative_peak_areas):
            error = random.gauss(0,float(self.specifications.noise_stdev))
                        
            rarea = float(rarea) + float(error)
            if float(rarea) < 0:
                L_relative_peak_areas[index]= 0
            else:
                L_relative_peak_areas[index]= str(rarea)
            
            
        #end for index, rarea in enumerate(L_relative_peak_areas)
        
        sum_areas = sum([float(rarea) for rarea in L_relative_peak_areas])
        L_relative_peak_areas = [str(float(rarea)*100.0/sum_areas) for rarea in L_relative_peak_areas]
        
        #replacing relative peak area values (added error)
        for index,a in enumerate(L_sorted):
            L_sorted[index][1] = L_relative_peak_areas[index]
        #end for
        
        print_string = "#Isomeric modforms and their percentage peak areas\n\n"
        
        print_string += "START\n"
        for l in L_sorted:
            
            l_string = ""
            for k,v in l[0].iteritems():
                l_string += "%s=%s "%(k,v)
            #end for   
            
            print_string += "IsomericModformsPeakArea: %s%.3f\n"%(l_string, round(float(l[1]),3))
        #end for l in L_sorted
        
        print_string += "END\n"
        
        fout = open(self.specifications.file_tdms_data, 'w')
        fout.write(print_string)
        fout.close()
        
    
    def find(self, seq_string, character, proline_check=False):
        """Instance method
        
        * This is a subpordinate method, i.e., only called within the same class by a member method
        * Returs a list of positions of amino-acids (single alphabet representation) in a string
        * Note that position starts from 1
        * Note that the method is case-sensitive
        
        Parameters
        ----------
        seq_string : str
          amino-acid sequence.
        
        character : str
          Character whose position(s) in the sequence is to be found.
        
        proline_check : bool
          If this is set to be ``True`` then positions of the amino-acids preceding proline will be ignored. Default is ``False``.
        

        Examples
        --------
        >>> from data_simulator import DataSimulator
        >>>
        >>> sim_data = DataSimulator()
        >>> 
        >>> aa_sequence = 'MSRYSRT'
        >>> L_arginine_locations = list(sim_data.find(aa_sequence, "R", proline_check=False))
        >>> print L_arginine_locations
        [3, 6]
        """
        for i,ltr in enumerate(seq_string):
            if ltr==character:
                if proline_check and seq_string[i+1]=="P":
                    continue
                else:
                    yield (i+1)
    

    def get_peptides_from_GluC_digestion(self, aa_sequence):
        """Instance method
        
        * This is a subordinate method, i.e., only called within the same class
        * The function returns a list of tuples of start and end locations of peptides resulting from Gluc digestion
        * GluC digests at the carboxyl side of glutamic acid (E) and aspartic acid (D)
        * Note that peptides strictly larger than 4 amino-acid residues are reported
        * List of peptides cross-examined using ExPaSy PeptideCutter: http://web.expasy.org/peptide_cutter/

        Parameters
        ----------
        aa_sequence : str
          amino-acid sequence
        
        Returns
        -------
        list

        Examples
        --------
        >>> from data_simulator import DataSimulator
        >>>
        >>> sim_data = DataSimulator()
        >>>
        >>> aa_sequence = 'MFLTRSEYDRGVNTFSPEGRLF'
        >>> L_tuples = sim_data.get_peptides_from_GluC_digestion(aa_sequence)
        >>> print L_tuples
        [(1, 7), (10, 18)]
        """

        L_tuples = []
        first    = 1

        L_E = list(self.find(aa_sequence,"E", proline_check=False)) # Glutamic acid
        L_D = list(self.find(aa_sequence,"D", proline_check=False)) # Aspartic acid
        
        L_ED = L_E + L_D
        L_ED.sort()
        
        for p in L_ED:
            tup = (first,p)
            length_peptide = p-first+1
            if length_peptide >= 5:
                L_tuples.append(tup)

            first = p+1

        return L_tuples
    
    
    def get_peptides_from_Trypsin_digestion(self, aa_sequence):
        """Instance method
        
        * This is a subordinate method, i.e., only called within the same class
        * This method returns a list of tuples of start and end locations of peptides resulting from Trypsin digestion
        * With Proline check - Trypsin digests at the carboxyl side of lysine or arginine unless they are followed by proline
        * Note that peptides strictly larger than 4 amino acids are reported
        * List of peptides cross-examined using ExPaSy PeptideCutter: http://web.expasy.org/peptide_cutter/

        Parameters
        ----------
          aa_sequence : str
             amino-acid sequence

        Returns
        -------
        list

        Examples
        --------
        >>> from data_simulator import DataSimulator
        >>>
        >>> sim_data = DataSimulator()
        >>>
        >>> aa_sequence = 'MFLTRSEYDRGVNTFSPEGRLF'
        >>> L_tuples = sim_data.get_peptides_from_Trypsin_digestion(aa_sequence)
        >>> print L_tuples
        [(1, 5), (6, 10), (11, 20)]
        """

        L_tuples = []
        first    = 1

        L_K = list(self.find(aa_sequence,"K", proline_check=True)) # Lysine
        L_R = list(self.find(aa_sequence,"R", proline_check=True)) # Arginine
        L_KR = L_K + L_R
        L_KR.sort()

        for p in L_KR:
            tup = (first,p)
            length_peptide = p-first+1
            if length_peptide >= 5:
                L_tuples.append(tup)

            first = p+1

        return L_tuples
    
    def whether_peptide_contains_any_modification_site(self, N_loc, C_loc):
        """Instance method
        
        * This is a subordinate method, i.e., only called within the same class by a member method
        * The method tells whether the peptide (N_loc, C_loc) contains any modification site
        
        Parameters
        ----------
        N_loc : str or int
          N-terminal location of the peptide in the aa sequence
        
        C_loc : str or int
          C-terminal location of the peptide in the aa sequence
        
        Returns
        -------
        bool
        
        Notes
        -----
        A given peptide specified by its ``N_loc`` and ``C_loc`` residue positions is examined for containing any modification site. The information on modification sites is contained in the dictionary ``D_sites_modtypes``, which stores the data from an input file. 
        """
        L_sites_sorted = sorted(self._D_sites_modtypes.keys(), key=lambda x: int(x[1:]))
        
        modification_site_presence = False
        for site in L_sites_sorted:
            if int(site[1:]) >= int(N_loc) and int(site[1:]) <= int(C_loc):
                modification_site_presence = True

        return modification_site_presence
    
    def get_dictionary_modification_frequency_for_peptide_in_modform(self, N_loc, C_loc, modform_object):
        """Instance method
        
        * This is a subordinate method, i.e., only called within the same class by a member method
        * This function provides ``D_peptide_mod_frequency`` for a peptide (``N_loc``, ``C_loc``) in a given protein modform
        
        Parameters
        ----------
        N_loc : str or int
          N terminal position of the peptide in the protein sequence
        
        C_loc : str or int
          C terminal position of the peptide in the protein sequence
        
        modform_object : :class:`.Modform`
          Modform object
          
        Returns
        -------
        D_peptide_mod_frequency : dict
          Python dictionary containing modifications mapped to their frequencies for a peptide (``N_loc``, ``C_loc``).

        Examples
        --------
        >>> from data_simulator import DataSimulator
        >>> from modform import Modform
        >>>
        >>> sim_data = DataSimulator()
        >>> mform = Modform(var_name='a15', D_protein_sites_modifications={'S15': 'P', 'T16': 'P', 'S54': 'P'})
        >>>
        >>> N_loc = 1; C_loc = 24 # peptide - (1, 24)
        >>> D_peptide_mod_frequency = sim_data.get_dictionary_modification_frequency_for_peptide_in_modform(N_loc, C_loc, mform)
        >>> print D_peptide_mod_frequency
        {'P': '2'}
        >>>
        >>> N_loc = 50; C_loc = 58 # peptide - (50, 58)
        >>> D_peptide_mod_frequency = sim_data.get_dictionary_modification_frequency_for_peptide_in_modform(N_loc, C_loc, mform)
        >>> print D_peptide_mod_frequency
        {'P': '1'}
        """
        D_peptide_mod_frequency = {}
                
        # for an unmodified modform
        if modform_object.D_protein_sites_modifications=={'0':'UNMOD'}:
            return {'UNMOD': '0'}
                   
        else:
            
            # dictionary of mod and frequency for this modform
            L_sites_modifications_sorted = sorted(modform_object.D_protein_sites_modifications.items(), key=lambda x: int(x[0][1:]))
            
            for tup in L_sites_modifications_sorted:

                # checking if the modification site lies within the peptide (N_loc, C_loc)
                if int(tup[0][1:]) >= int(N_loc) and int(tup[0][1:]) <= int(C_loc):
                    
                    if tup[1] in D_peptide_mod_frequency.keys():
                        freq = D_peptide_mod_frequency[tup[1]]
                        D_peptide_mod_frequency[tup[1]] = int(freq)+1
                                
                    else:
                        D_peptide_mod_frequency[tup[1]] = 1
                                
                else: #if int(tup[0]) >= int(N_loc)...
                    continue

            #end for tup in L_sites_modifications_sorted

        #end if-else
        
        if not D_peptide_mod_frequency:
            return {'UNMOD': '0'}
        else:
            return D_peptide_mod_frequency
    
    
    def get_dictionary_peptides_isomeric_peak_areas(self, D_variables_modforms, L_peptides):      
        """Instance method
        
        * This is a subordinate method, i.e., only called within the same class by a member method.
        * This method computes the percentage peak areas for isomeric peptides from their abundance and returns a dictionary of peptides mapped to their isomeric peak areas.
                
        Parameters
        ----------
        D_variables_modforms : dict
          Python dictionary of modform variables mapped to their :class:`.Modform` objects.
        
        L_peptides : list
          Python list of tuples of peptides, e.g., ``L = [(1, 6), (7, 15), (16, 30)]``. Each tuple is starting and ending position of a peptide in the protein sequence.
        
        Returns
        -------
        defaultdict(list)
        
        Examples
        --------
        >>> from data_simulator import DataSimulator
        >>> from modform import Modform
        >>>
        >>> D_variables_modforms = {}
        >>> D_variables_modforms['a1'] = Modform(var_name='a1', D_protein_sites_modifications={'0': 'UNMOD'})
        >>> D_variables_modforms['a2'] = Modform(var_name='a2', D_protein_sites_modifications={'S20': 'P'})
        >>> D_variables_modforms['a3'] = Modform(var_name='a3', D_protein_sites_modifications={'S15': 'P'})
        >>> D_variables_modforms['a4'] = Modform(var_name='a4', D_protein_sites_modifications={'S15': 'P', 'S20': 'P'})
        >>>
        >>>for k,v in D_variables_modforms.iteritems():
        ...    if k=='a1': v.percentage_abundance='0.000'
        ...    elif k=='a2': v.percentage_abundance='11.111'
        ...    elif k=='a3': v.percentage_abundance='44.444'
        ...    elif k=='a4': v.percentage_abundance='44.444'
        >>>
        >>> D_sites_modtypes = {'S20': ['P'], 'S15': ['P']}
        >>> sim_data = DataSimulator(D_sites_modtypes=D_sites_modtypes)
        >>> L_peptides = [(1, 5), (6, 10), (11, 20), (21, 32)]
        >>>
        >>> D = sim_data.get_dictionary_peptides_isomeric_peak_areas(D_variables_modforms, L_peptides)
        >>> print D
        defaultdict(<type 'list'>, {'11-20': [[{'UNMOD': '0'}, '0.000'], [{'P': 1}, '55.555'], [{'P': 2}, '44.444']]})
        """
        D_peptide_isomeric_peak_areas = defaultdict(list)
        
        # finding percentage peak areas for isomeric peptides
        for peptide in L_peptides:
            
            N_loc = peptide[0]
            C_loc = peptide[1]
            
            L_isomeric_peptides_peak_areas = [] # list of lists [..[D_peptide_mod_frequency, peak_area]..]
            
            # skip peptides which do not contain any modification sites
            if not self.whether_peptide_contains_any_modification_site(N_loc, C_loc):
                continue
            
            else:
                
                # find isomers of this peptide in across all modforms
                # add peak area for same D_peptide_mod_frequency map

                for var, mform in D_variables_modforms.iteritems():
                    D_peptide_mod_frequency = self.get_dictionary_modification_frequency_for_peptide_in_modform(N_loc, C_loc, mform)

                    #updating the percentage peak area
                    
                    percentage_peak_area_contribution = mform.percentage_abundance
                    if D_peptide_mod_frequency in [x[0] for x in L_isomeric_peptides_peak_areas]:
                            x_index = [x[0] for x in L_isomeric_peptides_peak_areas].index(D_peptide_mod_frequency)
                            p_area = L_isomeric_peptides_peak_areas[x_index][1]
                            percentage_peak_area = float(percentage_peak_area_contribution)+ float(p_area)
                            L_isomeric_peptides_peak_areas[x_index][1] = str(percentage_peak_area)
                    else:
                        L_isomeric_peptides_peak_areas.append([D_peptide_mod_frequency ,str(percentage_peak_area_contribution)])
                        
                    #end if-else
                #end for var, mform in D_variables_modforms.iteritems()
                                                    
                peptide_string = "%s-%s"%(peptide[0], peptide[1])
                D_peptide_isomeric_peak_areas[peptide_string] = L_isomeric_peptides_peak_areas
                
        #end for peptide in L_peptides
        
        return D_peptide_isomeric_peak_areas

    
    def generate_BUMS_MRM_data(self):
        """Instance method
        
        * This method computes peak areas for BUMS-MRM method
        * This is done only if the option ``FILE_BUMS_MRM_DATA`` is specified in the specification file. The content is then printed in the file specified.
        * User can choose ``Trypsin`` or ``Gluc`` or both in the specification file
        * Noise can also be introduced in the peak area values
        """
        D_variables_modforms = self._modform_space.D_variables_modforms
        aa_sequence = gbfunc.get_fasta_seq_from_the_file(self.specifications.file_fasta_sequence)
        
        L_peptides = []
        
        if self.specifications.Trypsin:
            L_peptides_Trypsin = self.get_peptides_from_Trypsin_digestion(aa_sequence)
            L_peptides = L_peptides + L_peptides_Trypsin

        if self.specifications.GluC:
            L_peptides_GluC    = self.get_peptides_from_GluC_digestion(aa_sequence)
            L_peptides = L_peptides + L_peptides_GluC 
        
        D_peptide_isomeric_peak_areas = self.get_dictionary_peptides_isomeric_peak_areas(D_variables_modforms, L_peptides)
        
        print_string =  "# Peptides and percentage peak areas of isomeric peptides\n"
        print_string += "# Modification_type=frequency; e.g., P=2 Ac=1 MM=2\n"
        print_string += "# Last column reserved for percentage peak area value\n\n"
        
        for p,A in D_peptide_isomeric_peak_areas.iteritems():
            
            print_string += "\nPEPSTART\n"
            print_string += "peptide: %s\n"%p
            
            A_sorted = sorted(A, key=lambda x: sum([int(a) for a in x[0].values()]))

            #adding measurment error in the simulated relative peak area values
            L_relative_peak_areas = [a[1] for a in A_sorted]
            for index, rarea in enumerate(L_relative_peak_areas):
                error = random.gauss(0,float(self.specifications.noise_stdev))
                rarea = float(rarea) + float(error)
                if float(rarea)<0:
                    L_relative_peak_areas[index] = 0
                else:
                    L_relative_peak_areas[index]= str(rarea)
            #end for index, rarea in enumerate(L_relative_peak_areas)
        
            sum_areas = sum([float(rarea) for rarea in L_relative_peak_areas])
            L_relative_peak_areas = [str(float(rarea)*100.0/sum_areas) for rarea in L_relative_peak_areas]
            
            #modifying relative peak area values (added error)
            for index,a in enumerate(A_sorted):
                A_sorted[index][1] = L_relative_peak_areas[index]
            #end for
            
            for l in A_sorted:
                l_string = ""
                for k,v in l[0].iteritems():
                    l_string += "%s=%s "%(k,v)
                
                print_string += "IsomericPeptidesPeakArea: %s%.3f\n"%(l_string, round(float(l[1]),3))
            print_string += "PEPEND\n"
        
        #end for p, A in D_peptide_isomeric_peak_areas.iteritems()
        
        fout = open(self.specifications.file_bums_mrm_data, 'w')
        fout.write(print_string)

        fout.close()
 
