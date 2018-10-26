from itertools import product
from itertools import combinations

import re
import pulp

import modform       as modform
import globalfuncs   as gbfunc
import modform_space as mform_space

class LinearConstraints(object):
    
    """ This class provides framework to store python list of linear constraints in string format as well in PuLP format.
    
    Parameters
    ----------
    D_sites_modtypes : dict
      Python dictionary containing residues mapped to the list of possible modifications. E.g., ``D = {'T32': ['P'], 'S55': ['Ac', 'P']}``
    
    specifications : :class:`.ModformProSpecifications`
      Instance of ModformProSpecifications object.
    """
    
    
    def __init__(self, D_sites_modtypes, specifications=None):
        
        self._L_linear_constraints  = [] # da: will remove later if found futile
        
        self._L_PuLP_constraints    = []
                        
        self._D_sites_modtypes      = D_sites_modtypes
        self._specifications        = specifications
        
        # data structure for grouped variables
        
        self._D_grouped_variables                    = {}
        
        self._L_grouped_variables_linear_constraints = []
        self._L_grouped_variables_PuLP_constraints   = []

        # modform space
        self._modform_space = mform_space.ModformSpace(specifications)
      
    @property
    def L_linear_constraints(self):
        """ Instance attribute
        
        * Python list of linear constraints in string format
        * This is populated by a member method of this class
        * E.g., ``L = [('a1 + a3 + a5', '30'), ('a2 + a4 + a6', '70')]``
        """
        return(self._L_linear_constraints)

    @L_linear_constraints.setter
    def L_linear_constraints(self, list_constraints):
        self._L_linear_constraints = list_constraints

    @property
    def L_PuLP_constraints(self):
        """Instance attribute
        
        * Python list of linear constraints in PuLP format
        * PuLP format is needed to use pulp library for linear programming.
        * This is populated by a member method of this class
        
        >>> # Demo L_PuLP_constraints
        [LpConstraint(LpVariable(a1) + LpVariable(a3) + LpVariable(a5), 30.00), LpConstraint(LpVariable(a2) + LpVariable(a4) + LpVariable(a6), 70.00)]``
        """
        return(self._L_PuLP_constraints)

    @L_PuLP_constraints.setter
    def L_PuLP_constraints(self, list_PuLP_constraints):
        self._L_PuLP_constraints = list_PuLP_constraints
        

    @property
    def specifications(self):
        """Instance attribute
        
        Instance of a :class:`.ModformProSpecifications` object.
        """
        return self._specifications
        
    @property
    def D_sites_modtypes(self):
        """Instance attribute
        
        Python dictionary having aa residues mapped to the list ofpossible modifications on that residue.
        """
        return self._D_sites_modtypes
    
    @D_sites_modtypes.setter
    def D_sites_modtypes(self, D_sites_modtypes):
        self._D_sites_modtypes = D_sites_modtypes

    @property
    def D_grouped_variables(self):
        """Instance attribute
        
        * Python dictionary having a new coalesced variables mapped to the group of co-existing variables.
        * This is populated in a member method if ``VARIABLES_GROUPING`` is ``True`` in the specification file.
        * E.g., ``D = {'b1': 'a1 + a4', 'b2': 'a3 + a5 + a6'}``
        """
        return self._D_grouped_variables


    @D_grouped_variables.setter
    def D_grouped_variables(self, D_grouped_variables):
        self._D_grouped_variables = D_grouped_variables
    
    
    @property
    def L_grouped_variables_linear_constraints(self):
        """
        * Instance attribute
        * Python list of linear constraints with new coalesced variables
        * This is populated by a member method of this class.
        """
        return self._L_grouped_variables_linear_constraints

    @L_grouped_variables_linear_constraints.setter
    def L_grouped_variables_linear_constraints(self, L_grouped_variables_linear_constraints):
        self._L_grouped_variables_linear_constraints = L_grouped_variables_linear_constraints
    
    @property
    def L_grouped_variables_PuLP_constraints(self):
        """Instance attribute
        
        * Python list of linear constraints in PuLP format with new coalesced variables
        * This is populated by a member method of this class.
        """
        return self._L_grouped_variables_PuLP_constraints
    
    @L_grouped_variables_PuLP_constraints.setter
    def L_grouped_variables_PuLP_constraints(self, L_grouped_variables_PuLP_constraints):
        self._L_grouped_variables_PuLP_constraints   = L_grouped_variables_PuLP_constraints
    
    @property
    def modform_space(self):
        """Instance attribute
        
        * Instance of a :class:`.ModformSpace` object.
        * This needs to be defined for real data.
        * Predefined modform_space is used when running in simulation mode.
        """
        return self._modform_space
    
    @modform_space.setter
    def modform_space(self, modform_space):
        self._modform_space = modform_space

    def get_grouped_variables(self):
        """Instance method
        
        This function returns a python list of grouped variables.
        
        Returns
        -------
        list
        """
        
        L_list_of_lhs_variables = []
        for constraint in self._L_linear_constraints:
            lhs_string    = constraint[0]
            lhs_variables = [v.strip() for v in re.split(r'[+-]', lhs_string)]
            
            L_list_of_lhs_variables.append(lhs_variables)
        
        #end for constraint in self._L_linear_constraints
        
        # defining list of grouped sets of variables
        L_grouped_sets = []
        
        for variable in self._modform_space.D_variables_modforms.keys():
            variable = variable.strip()
            
            L_set_of_colleague_variables = []
            for lhs_vars in L_list_of_lhs_variables:
                
                if variable in lhs_vars: 
                    #print "variable, constraint_lhs: ",  variable, constraint_lhs
                    L_set_of_colleague_variables.append(set(lhs_vars))
                else: continue

            #end for lhs_vars in L_list_of_lhs_variables
            
            S = set.intersection(*L_set_of_colleague_variables)
            
            #grouping_all_constraints = True
            #for lhs_vars in L_list_of_lhs_variables:
            #    common =  S.intersection(set(lhs_vars))
            #    if len(common)>0 and len(common)<len(S):
            #        grouping_all_constraints=False
            
            #if len(S)>1 and S not in L_grouped_sets and grouping_all_constraints==True:
            #    L_grouped_sets.append(S)
            
            list_all_size_sets = list(gbfunc.powerset(S))
                        
            for sset in sorted(list_all_size_sets, key=lambda x:len(x), reverse=True):
                if len(sset)<2 or variable not in list(sset): continue
                else:
                    grouping_all_constraints = True;  count=0
                    for lhs_vars in L_list_of_lhs_variables:
                        common =  set(sset).intersection(set(lhs_vars))
                        if len(common)>0 and len(common)<len(sset):
                            grouping_all_constraints=False
                        elif len(common)>0:
                            count +=1
                        else:
                            continue
                        #end if-elif-else
                    #end for lhs_vars in L_list_of_lhs_variables
                    
                    if len(sset)>1 and count >1  and grouping_all_constraints==True:
                        L_binary = [set(sset)<=s for s in L_grouped_sets]
                        if True in L_binary: break
                        else:
                            L_grouped_sets.append(set(sset))
                            break
                            
        #end for variable in self._D_variables_modforms.keys()
        
        # sorting before returning
        L_grouped_lists = [sorted(v, key=gbfunc.natural_keys) for v in L_grouped_sets]
                
        return L_grouped_lists

    
    def get_sum_of_variables_string(self, L_variables):
        """Instance method

        * This is a subordinate method, i.e., only called in a member method of the same class.
        * The function returns a string of sum of variables given in the input list.

        Parameters
        ----------
        L_variables : list
          Python list of modform variables in string format.

        Returns
        -------
        str
        """
        if len(L_variables)==0:
            return ''
        
        L_variables = sorted(L_variables, key=gbfunc.natural_keys)
        
        sum_string ='%s'%L_variables[0]

        for i in range(1,len(L_variables)):

            sum_string += ' + %s'%L_variables[i]
            
        #end for i in range(1, len(L_variables))
        
        return sum_string
    
    
    def get_linear_constraints_with_grouped_variables(self, L_grouped_lists):
        """Instance method
        
        * This method creates List of linear constraints with grouped variables both in string as well as PuLP format thereby, populating the instance attributes, ``L_grouped_variables_linear_constraints``, ``L_grouped_variables_PuLP_constraints``
        * Note that additional variables are added to the existing dictionary ``D_PuLP_variables``, no additional dictionary is needed for grouped variables.
        
        Parameters
        ----------
        L_grouped_lists : list
          Python list of grouped variables. E.g., ``L = ['a1 + a3', 'a2 + a4 + a6']``
        """
        
        # initialising the list
        self._L_grouped_variables_linear_constraints = list(self._L_linear_constraints)
        self._L_grouped_variables_PuLP_constraints   = list(self._L_PuLP_constraints)
        
        count = 1
        for glist in L_grouped_lists:
            grouped_var_name = 'b%s'%count
            
            grouped_string = self.get_sum_of_variables_string(glist)

            self._D_grouped_variables[grouped_var_name] = grouped_string

            #delete keys from the self._D_PuLP_variables
            
            for k in [v.strip() for v in glist]:
                self._modform_space.D_PuLP_variables.pop(k, None)
            
            self._modform_space.D_PuLP_variables[grouped_var_name] = pulp.LpVariable(grouped_var_name,lowBound=0)
            
            count += 1
                                
            for ii, constraint in enumerate(self._L_linear_constraints):
                lhs_string = constraint[0]
                rhs_string = constraint[1]
                L_lhs_variables   = [v.strip() for v in re.split(r'[+-]', lhs_string)]
                set_lhs_variables = set(L_lhs_variables)
                
                if set_lhs_variables.intersection(set(glist))==set(glist):
                    set_diff = set_lhs_variables.difference(set(glist))
                    lhs_string_new = self.get_sum_of_variables_string(list(set_diff))
                    if lhs_string_new=='':
                        lhs_string_new = '%s'%grouped_var_name
                    else:
                        lhs_string_new += ' + %s'%grouped_var_name
                    
                    # replacing the existing linear constraint with an appropriate one
                    self._L_linear_constraints[ii] = (lhs_string_new, rhs_string)
                    
                    #replacing the existing linear constraint with an appropriate one
                    # assuming index are same for identical linear constraints
                    
                    LHS = pulp.lpSum(self._modform_space.D_PuLP_variables[var.strip()] for var in re.split(r'[+-]', lhs_string_new))
                    self._L_PuLP_constraints[ii] = pulp.LpConstraint(LHS,rhs=float(rhs_string))
                else:
                    continue

                #end if-else set_lhs_variables.intersection(gset)==gset
            #end for ii, constraint in enumerate(L_grouped_variables_linear_constraints)
        #end for gset in L_grouped_sets
                
    
    
    def define_modform_space_for_linear_constraints(self, D_variables_modforms_simulated=None, D_PuLP_variables_simulated=None):
        """Instance method
        
        * Its a subordinate method, only called within the same class.
        * This function defines modforms space for linear constraints. It updates the instance attribute, ``modform_space``
        * For simulation mode, predefined ``modform_space`` object is copied.
        * For real mode, a fresh ``modform_space`` object is created.
        
        Parameters
        ----------
        D_variables_modforms_simulated : dict
          Python dictionary containing modform variables mapped to their :class:`.Modform` objects. This is non-empty for simulation mode.
        
        D_PuLP_variables_simulated : dict
          Python dictionary containinng modform variables mapped to their :class:`.pulp` objects. This is non-empty for simulation mode.
        """

        if self._specifications.data_mode=='simulation':
            self._modform_space.D_variables_modforms = dict(D_variables_modforms_simulated)
            self._modform_space.D_PuLP_variables     = dict(D_PuLP_variables_simulated)

        elif self._specifications.data_mode=='real':
            self._modform_space.generate_modform_space()
        
    
    def derive_linear_constraints_from_MS_measurments(self, L_isomeric_modforms_peak_areas=None, D_peptide_isomeric_peak_areas=None, D_variables_modforms_simulated=None, D_PuLP_variables_simulated=None):
        """Instance method
        
        * This method derives the linear constraints from available peak area measurments either from BUMS-MRM or TDMS or both.
        * It populates the instance attributes, ``L_linear_constraints``, ``L_PuLP_constraints``
        
        Parameters
        ----------
        L_isomeric_modforms_peak_areas : list
          Python list containing peak areas of isomeric modforms. This data typically comes from TDMS measurements.
        
        D_peptide_isomeric_peak_areas : dict
          Python dictionary containing peak areas of peptides in their different isomeric states. This data typically comes from BUMS-MRM measurments.
        
        D_variables_modforms_simulated : dict
          Python dictionary containing modform variables mapped to their :class:`.Modform` objects for simulated data. This ensures usage of predefined modform variables for the linear constraints.
        
        D_PuLP_variables_simulated : dict
          Python dictionary containing modform variables mapped to their :class:`.pulp` objects for simulated data. This ensures usage of predefined modform variables for the linear constraints.
        """
        #defining modform space
        self.define_modform_space_for_linear_constraints(D_variables_modforms_simulated=D_variables_modforms_simulated, D_PuLP_variables_simulated=D_PuLP_variables_simulated)
        
        
        # dictionary of variables associated to modform objects
        # assuming lists are sorted -- UNMOD, Ac, M, MM, MMM, P, U, UU, UUU
        if self._specifications.file_tdms_data:
            print "#####TDMS LCs####"
            for modform_isomer in L_isomeric_modforms_peak_areas:
                D_vars_modforms_isomeric, D_PuLP_vars_isomeric = self._modform_space.get_isomeric_modforms(modform_isomer[0])
            
                # constructing linear equations for isomeric modforms
                L_sorted_vars = sorted(D_vars_modforms_isomeric.iterkeys(), key=gbfunc.natural_keys)
                                
                lhs_string_isomeric = L_sorted_vars[0]
                
                for k in range(1, len(L_sorted_vars)):
                    
                    lhs_string_isomeric = lhs_string_isomeric + ' + %s'%L_sorted_vars[k]
                #end for k in range(1, len(L_sorted_vars))
                
                rhs_string_isomeric = modform_isomer[1]
                print "TDMS Eq: %s = %s"%(lhs_string_isomeric, rhs_string_isomeric)
                self._L_linear_constraints.append((lhs_string_isomeric,rhs_string_isomeric))
                LHS = pulp.lpSum(v for v in D_PuLP_vars_isomeric.values())
                self._L_PuLP_constraints.append(pulp.LpConstraint(LHS, rhs=float(rhs_string_isomeric)))
            
            #end for modform_isomer in L_isomeric_modforms_peak_areas
        
        '''
        D_peptide_isomeric_peak_areas = {'15-25':[({'UNMOD':'0'},'25'),({'P':'1'},'50'),({'P':'2'},'25')],
        '50-60':[(['UNMOD':'0'},'50'),({'P':'1','Ac':'1'},'30'),({'P':'2'},'20')]}
        '''
        #print "D_peptide_isomeric_peak_areas: ", D_peptide_isomeric_peak_areas
        
        if self._specifications.file_bums_mrm_data:
            print "#####BUMS LCs####"
            for peptide, l_isomeric_peak_areas in D_peptide_isomeric_peak_areas.iteritems():
                
                N_loc, C_loc = peptide.split('-')
            
                for pep_isomer in l_isomeric_peak_areas:
                
                    # list of modform variables yielding the following isomer
                    L_modform_var = []
                    
                    for var, mform in  self._modform_space.D_variables_modforms.iteritems():
                        if mform.whether_modified_peptide_present_in_modform(N_loc, C_loc, D_peptide_mod_frequency=pep_isomer[0]):
                            L_modform_var.append(var)
                        
                        else:
                            continue
                    
                        #end if-else   
                    #end for var, mform
                    
                    # constructing a linear equation
                    L_modform_var = sorted(L_modform_var, key=gbfunc.natural_keys)
                    lhs_string_isomeric = L_modform_var[0]
                    
                    for i in range(1,len(L_modform_var)):
                        lhs_string_isomeric = lhs_string_isomeric + ' + %s'%L_modform_var[i]
                    #end for i in range..
                    
                    rhs_string_isomeric = pep_isomer[1]
                    print "BUMS Eq: %s = %s"%(lhs_string_isomeric, rhs_string_isomeric)
                    self._L_linear_constraints.append((lhs_string_isomeric, rhs_string_isomeric))
                    
                    LHS = pulp.lpSum(self._modform_space.D_PuLP_variables[v]  for v in L_modform_var )
                    self.L_PuLP_constraints.append(pulp.LpConstraint(LHS, rhs=float(rhs_string_isomeric)))

                #end for pep_isomer in l_isomeric_peak_areas
            #end for peptide, l_isomeric_peak_areas..
            
    
