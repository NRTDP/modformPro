import pulp

from itertools import combinations
from itertools import product

import globalfuncs as gbfunc
import modform     as modform


class ModformSpace(object):
    
    """This class defines a modform space, i.e., set of modform objects.
    
    Parameters
    ----------
    specifications : :class:`.ModformProSpecifications`, optional
      - Instance of ModformProSpecifications object
      - Contains information about the option set by the user in the specification file
    
    L_isomeric_modforms_modes : list, optional
      - Python list of python dictionaries having modification types mapped to its frequency.
      - E.g.,  ``L = [{'UNMOD': '0'}, {'P': '1'}, {'Ac': '1', 'P': '2'}]``
        There are three isomeric modform modes in the list - unmodified, singly phosphorylated modforms and the last is referring to the modforms that are singly acetylated and doubly phosphorylated.
                
    """
    def __init__(self, specifications=None, L_isomeric_modforms_modes=None):
        
        self._specifications            = specifications
        self._L_isomeric_modforms_modes = L_isomeric_modforms_modes

                
        # Data structures to store modforms
        self._D_variables_modforms   = {}
        self._D_PuLP_variables       = {}
        
        self._D_variables_nmolecules = {}
        
        # reading and storing the apriori information - sites and modifications
        if specifications:
            self._D_sites_modtypes = gbfunc.read_sites_and_modification_types_from_file(specifications.file_modifications)
    
    @property
    def D_variables_modforms(self):
        """Instance attribute
        
        Python dictionary with modform variables (string format) being mapped to their respective :class:`.Modform` objects.
        """
        return self._D_variables_modforms

    @D_variables_modforms.setter
    def D_variables_modforms(self, D_variables_modforms):
        self._D_variables_modforms = D_variables_modforms

    @property
    def D_PuLP_variables(self):
        """Instance attribute
        
        * Python dictionary having modform variables (string format) being mapped to their respective :class:`.pulp` objects.
        * Note that modform variables in PuLP format are needed to perform linear programming computations with pulp library.
        """
        return self._D_PuLP_variables

    @D_PuLP_variables.setter
    def D_PuLP_variables(self, D_PuLP_variables):
        self._D_PuLP_variables = D_PuLP_variables

    @property
    def D_variables_nmolecules(self):
        """Instance attribute
        
        * Python dictionary is populated only during simulated data mode. Infact modform_space object is created while simulating modform distribution.
        * This has modform variables (in string format) being mapped to the number of molecules assigned to each modform variable. 
        """
        return self._D_variables_nmolecules

    @D_variables_nmolecules.setter
    def D_variables_nmolecules(self, D_variables_nmolecules):
        self._D_variables_nmolecules = D_variables_nmolecules

    
    
    def generate_modform_space(self):
        """Instance method
        
        * It defines the modform space (without specifying the percentage abundance of the modforms) and populate the instance object attributes ``D_variables_modforms``, ``D_PuLP_variables``.
        * If ``n_sites_full_modform_space`` option is set to be ``True`` then all the possible modforms are generated for n_sites and modification types.
        * This space will be restricted complying with the intact mass spectrometry data or user provides modes (for data simulation) if the option n_sites_full_modform_space is set to be False.
        * This method is called during simulation of modform distribution and also needed to write linear constraints from peak area calculations.

        Notes
        -----
        If the user option ``N_SITES_FULL_MODFORM_SPACE`` in the specification file is set to be ``True`` then full modform space is generated complying with the number of sites and their respective modification types. 
        If the above user option is set to be False, then it will generate a modform space based on the isomeric modform modes either derived from intact mass spectrometry data or in case of "simulation" data mode, the modes would be specified in the ``MODES`` field of the specification file.
        """
        D_variables_modforms = {}
        D_PuLP_variables     = {}
        
        if self._specifications.n_sites_full_modform_space:
            
            n_sites = len(self._D_sites_modtypes)
            start_index = 0
            for m in range(n_sites+1):
                
                L_D_protein_mod_frequency = self.get_list_of_modification_frequency_map(m)
                
                for D_protein_mod_frequency in L_D_protein_mod_frequency:
                    #determining the start_index value
                    if D_protein_mod_frequency.items()==[('UNMOD','0')] or m==0:
                        start_index=1
                    elif len(self._D_variables_modforms)==0:
                        start_index = 2
                    else:
                        # last index of modform in the dictionary of variables
                        last_index = sorted(self._D_variables_modforms.keys(), key=gbfunc.natural_keys)[-1].split('a')[1]
                        start_index = int(last_index) + 1
                    
                    #end if-elif-else D_protein_mod_frequency.items()==[('UNMOD','0')]..
                    
                    (D_vars_mforms, D_PuLP_vars)= gbfunc.get_modforms(start_index, D_sites_modtypes=self._D_sites_modtypes,D_protein_mod_frequency=D_protein_mod_frequency)
                    self._D_variables_modforms.update(D_vars_mforms)
                    self._D_PuLP_variables.update(D_PuLP_vars)

                    #Assigning number of molecules for modforms for simulation
                    if self._specifications.data_mode.lower()=='simulation' and self._L_isomeric_modforms_modes!=None:
                        
                        if D_protein_mod_frequency in self._L_isomeric_modforms_modes:
                            
                            D_vars_nmols = self.get_dictionary_of_variables_molecules(D_vars_mforms, state="principal")
                            self._D_variables_nmolecules.update(D_vars_nmols)
                            
                        else:
                            D_vars_nmols = self.get_dictionary_of_variables_molecules(D_vars_mforms, state="null")
                            self._D_variables_nmolecules.update(D_vars_nmols)
                        #end if-else
                    #end if self._specifications.data_mode..
                #end for D_protein_mod_frequency in L_D_protein_mod_frequency
            #end for m in range(n_sites+1)
            
        else: #modform space is restricted to the isomeric modform modes
            print "generating restricted modform space using: ", self._L_isomeric_modforms_modes
            start_index = 0
            for mode in self._L_isomeric_modforms_modes:

                total_freq = sum([int(a) for a in mode.values()])
                
                if mode.items()==[('UNMOD','0')] or total_freq==0:
                    start_index=1
                elif len(self._D_variables_modforms)==0:
                    start_index = 2
                else:
                    # last index of modform in the dictionary of variables
                    last_index = sorted(self._D_variables_modforms.keys(), key=gbfunc.natural_keys)[-1].split('a')[1]
                    start_index = int(last_index) + 1
                
                #end if-elif-else mode.items()==[('UNMOD','0')]..
                
                (D_vars_mforms, D_PuLP_vars)= gbfunc.get_modforms(start_index, D_sites_modtypes=self._D_sites_modtypes,D_protein_mod_frequency=mode)
                self._D_variables_modforms.update(D_vars_mforms)
                self._D_PuLP_variables.update(D_PuLP_vars)
                
                if self._specifications.data_mode=="simulation":
                    D_vars_nmols = self.get_dictionary_of_variables_molecules(D_vars_mforms, state="principal")
                    self._D_variables_nmolecules.update(D_vars_nmols)
                                        
                #end if self._specifications.data_mode=="simulation"       
            #end for mode in self._L_isomeric_modforms_modes
        #end if-else self._specifications.n_sites_full_modform_space

    
    def get_isomeric_modforms(self, D_protein_mod_frequency):
        """Instance method
        
        * It returns two python dictionaries having isomeric modforms variables (string format) mapped to modform objects and their PuLP counterparts.
        * These isomeric modforms are subset of the entire modform space.
        * The subset is chosen based on the isomeric modform mode which is defined by the dictionary D_protein_mod_frequency.

        Parameters
        ----------
        D_protein_mod_frequency: dict
          Python dictionary having modifications mapped to their frequency.
        
        Returns
        -------
        dict, dict
        """
        D_vars_modforms_isomeric = {}
        D_PuLP_vars_isomeric     = {}
        
        for var,mform in self._D_variables_modforms.iteritems():

            D_mod_freq = mform.get_modifications_frequency_map()
            if D_protein_mod_frequency==D_mod_freq:
                D_vars_modforms_isomeric[var]=mform
                D_PuLP_vars_isomeric[var] = self._D_PuLP_variables[var]
                
            else:
                continue

            #end if-else
        #end for
        
        return (D_vars_modforms_isomeric, D_PuLP_vars_isomeric)

    
    def get_dictionary_of_variables_molecules(self, D_variables_modforms, state="principal"):
        """Instance method
        
        * It assigns number of molecules to a given set of modform variables and returns a python dictionary of the same.
        * This method is only called in a member method of the class :class:`.DataSimulator`.
        
        Parameters
        ----------
        D_variables_modforms : dict
          Python dictionary of modforms variables (string format) mapped to the modform objects. 
        
        state : str, "principal" (default) or "subordinate" or "null"
          They determine number of molecules assigned to the modforms. 
          "Null" assigns 0 molecules to the modforms.
          Default is "principal".
        
        Returns
        -------
        D_variables_nmolecules : dict
          Python dictionary having modform variables names (str) mapped to the number of molecules assigned.
        
        Examples
        --------
        >>> from modform import Modform
        >>> from modform_space import ModformSpace
        >>>
        >>> D_variables_modforms = {}
        >>> D_variables_modforms['a1'] = Modform(var_name='a1', D_protein_sites_modifications={'0': 'UNMOD'})
        >>> D_variables_modforms['a2'] = Modform(var_name='a2', D_protein_sites_modifications={'S20': 'P'})
        >>> D_variables_modforms['a3'] = Modform(var_name='a3', D_protein_sites_modifications={'S15': 'P'})
        >>> D_variables_modforms['a4'] = Modform(var_name='a4', D_protein_sites_modifications={'S15': 'P', 'S20': 'P'})
        >>> 
        >>> mform_space = ModformSpace()
        >>>
        >>> D_nmolecules_principal = mform_space.get_dictionary_of_variables_molecules(D_variables_modforms, state="principal")
        >>> print D_nmolecules_principal
        {'a1': 1000, 'a3': 250, 'a2': 100, 'a4': 100}
        >>>
        >>> D_nmolecules_subordinate= mform_space.get_dictionary_of_variables_molecules(D_variables_modforms, state="subordinate")
        >>> print D_nmolecules_subordinate
        {'a1': 250, 'a3': 100, 'a2': 100, 'a4': 100}
        >>>
        >>> D_nmolecules_null = mform_space.get_dictionary_of_variables_molecules(D_variables_modforms, state="null")
        >>> print D_nmolecules_null
        {'a1': 0, 'a3': 0, 'a2': 0, 'a4': 0}
        """
        D_variables_nmolecules = {}

        # principal set of modforms
        if state=="principal":
            for i, var in enumerate(D_variables_modforms.keys()):
                #there is no ordering of items in the dictionary so they are assigned the number of molecules randomly
                nmol  = 0
                if   i==0: nmol=1000
                elif i==1: nmol=250
                elif i<6:  nmol=100
                else: nmol=50
                
                D_variables_nmolecules[var] = nmol
        
        # subordinate set of modforms
        elif state=="subordinate":
            for j, mvar in enumerate(D_variables_modforms.keys()):
                #there is no ordering of items in the dictionary so they are assigned the number of molecules randomly
                nmol  = 0
                if   j==0: nmol=250
                elif j<=6: nmol=100
                else: nmol=0
                
                D_variables_nmolecules[mvar] = nmol
            #end for j, mvar in enumerate(D_variables_modforms.keys())
        
        # null set of modforms
        elif state=="null":
            for j, mvar in enumerate(D_variables_modforms.keys()):
                #there is no ordering of items in the dictionary so they are assigned the number of molecules randomly
                nmol  = 0
                D_variables_nmolecules[mvar] = nmol
                
            #end for j, mvar in enumerate(D_variables_modforms.keys())
        #end if-elif state=="principal"
        
        return D_variables_nmolecules

    

    def get_list_of_modification_frequency_map(self, total_mod_frequency):
        """Instance method
        
        * This is a subordinate method, i.e. called inside this class method.
        * For a given total modification frequency, this function returns list of modifications to frequency map, L_D_protein_mod_frequency.
               
        Parameters
        ----------
        total_mod_frequeny:  int
          Total frequency of modifications or total number of sites being modified.
        
        Returns
        -------
        L_D_protein_mod_frequency : list
          Python list of python dictionary having modifications to their frequency.
        
        Examples
        --------
        >>> D_sites_modtypes = {'15': ['P'], '20': ['P'], '55': ['P', 'Ac']}
        >>>
        >>> mform_space = ModformSpace()
        >>>
        >>> total_mod_frequency = 0
        >>> L_D_protein_mod_frequency = mform_space.get_list_of_modification_frequency_map(total_mod_frequency)
        [{'UNMOD':'0'}]
        >>>
        >>> total_mod_frequency = 1
        [{'Ac': '1'}, {'P': '1'}]
        >>>
        >>> total_mod_frequency = 2
        [{'P': '1', 'Ac': '1'}, {'P': '2'}]
        >>>
        >>> total_mod_frequency = 3
        [{'P': '2', 'Ac': '1'}]
        """

        if total_mod_frequency==0:
            return [{'UNMOD':'0'}]
        else:
            L_D_protein_mod_frequency = []
        
            L = self._D_sites_modtypes.values()
        
            for l in combinations(L,total_mod_frequency):
    
                L_p = list(product(*l))
                
                for p in L_p:
                    D_m = {}
                    L_mod = list(set(p))
                    
                    for mod in L_mod:
                        c = p.count(mod)
                        D_m[mod] = str(c)
                    #end for mod in L_mod
                    if D_m in L_D_protein_mod_frequency: continue
                    else: L_D_protein_mod_frequency.append(D_m)
                #end for p in L_p
            
                
            #end for l in combinations(L, total_mod_frequency)
        
            return L_D_protein_mod_frequency
        
        #end if-else total_mod_frequency==0
        
