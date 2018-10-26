import re
import pulp

import modform as modform
import linear_constraints as lc
import linear_programming_solver as lps

def natural_keys(text):
    return [int(c) if c.isdigit() else c for c in re.split('(\d+)', text)]


class FeasibleLinearConstraints(object):

    """This class provides attributes and methods to obtain feasible linear constraints by removing the inconsistency introduced possibly due to measurment errors.
    
    Parameters
    ----------
    linear_constraints : :class:`.LinearConstraints`
      Set of linear constraints.
    """
    
    
    def __init__(self, linear_constraints):
        
        self._linear_constraints = linear_constraints

        # Initializing feasible linear constraints object
        self._feasible_linear_constraints = lc.LinearConstraints(linear_constraints.D_sites_modtypes, linear_constraints.specifications)

        # Using the same modform space; only variable names are of importance.
        self._feasible_linear_constraints.modform_space = linear_constraints.modform_space

        
    @property
    def feasible_linear_constraints(self):
        """
        * Instance attribute
        * New instance of :class:`.LinearConstraints` object.
        * It is populated by member method
        """
        return self._feasible_linear_constraints

    @feasible_linear_constraints.setter
    def feasible_linear_constraints(self, feasible_lc):
        self._feasible_linear_constraints = feasible_lc

    
    def get_feasible_linear_constraints(self):
        """Instance method
        
        * This method obtains set of feasible linear constraints and populated the instance attribute, ``feasible_linear_constraints``.
        * Linear programming algorithm is employed with slack variables
        """
        E_constraints = self.get_elastic_constraints()

        # trying something new
        ###############################################
        # additional constraints
        
        L_evariables_sorted = sorted([k for k,v in E_constraints.modform_space.D_PuLP_variables.iteritems() if k[0]=='e'], key=natural_keys)
        L_PuLP_evariables_sorted = [E_constraints.modform_space.D_PuLP_variables[k] for k in L_evariables_sorted]
        
        L_extra_PuLP_constraints = [] # e.g., -t <= e_1-e_2 <= t
        t_var = pulp.LpVariable('t', lowBound=0)
        
        for ii in xrange(0, len(L_PuLP_evariables_sorted), 2):
    
            e1 = L_evariables_sorted[ii]
            e2 = L_evariables_sorted[ii+1]
            
            #e1_e2 - t <= 0; sense=-1 sets LE-- less than equal to inequality
            LHS_1 = E_constraints.modform_space.D_PuLP_variables[e1] - E_constraints.modform_space.D_PuLP_variables[e2] - t_var  
            L_extra_PuLP_constraints.append(pulp.LpConstraint(LHS_1,rhs=0, sense=-1))
            
            # e1-e2 +t >= 0; sense=1 sets GE-- greater than equal to inequality
            LHS_2 = E_constraints.modform_space.D_PuLP_variables[e1] - E_constraints.modform_space.D_PuLP_variables[e2] + t_var  
            L_extra_PuLP_constraints.append(pulp.LpConstraint(LHS_2,rhs=0, sense=1))
        #end for ii in xrange(0, len(L_PuLP_evariables_sorted), 2)
            
        #L_PuLP_constraints_extended = list(E_constraints.L_PuLP_constraints)
        E_constraints.L_PuLP_constraints.extend(L_extra_PuLP_constraints)
        
        #objective_function_PuLP_extended   = t_var + pulp.lpSum(L_PuLP_evariables_sorted)
        
        
        #################################################
        
        
        lp_solver = lps.LinearProgrammingSolver(E_constraints)
        
        #L_variables_sorted = sorted([k for k,v in E_constraints.D_PuLP_variables.iteritems() if k[0]=='e'], key=natural_keys)
        #L_PuLP_evariables_sorted = [E_constraints.D_PuLP_variables[k] for k in L_variables_sorted]
        
        objective_function_PuLP = pulp.lpSum(L_PuLP_evariables_sorted) + t_var
        
        lp_minimize_solution = lp_solver.linear_programming_solver(objective_function_PuLP, pulp.LpMinimize)
        
        #print(lp_minimize_solution)
        #print("status: ", LpStatus[lp_minimize_solution.status])
        
        #for v in lp_minimize_solution.variables():
        #    print(v, pulp.value(v))
        
        for i in range(1,len(self._linear_constraints.L_linear_constraints)+1):

            elastic_variable1 = "e_%s"%(2*i-1)
            elastic_variable2 =  "e_%s"%(2*i)

            value_PuLP_evariable1 = pulp.value(E_constraints.modform_space.D_PuLP_variables[elastic_variable1]) 
            
            value_PuLP_evariable2 = pulp.value(E_constraints.modform_space.D_PuLP_variables[elastic_variable2])
            diff         = value_PuLP_evariable1 - value_PuLP_evariable2
            adjusted_rhs = float(self._linear_constraints.L_linear_constraints[i-1][1]) - diff
            
            self._feasible_linear_constraints.L_linear_constraints.append((self._linear_constraints.L_linear_constraints[i-1][0], adjusted_rhs))
            
            lhs_PuLP = sum([v[0]*v[1] for v in self._linear_constraints.L_PuLP_constraints[i-1].items()])
            self._feasible_linear_constraints.L_PuLP_constraints.append(pulp.LpConstraint(lhs_PuLP, rhs=float(adjusted_rhs))) 

    
    def get_elastic_constraints(self):
        """Instance method
        
        This method creates an instance of the class LinearConstraints where elastic/slack variables are added.
        
        Returns
        -------
        E_linear_constraints : :class:`.LinearConstraints`
          Linear constraints with added pair of elastic (slack) variables to each constraint.
        """
        E_linear_constraints = lc.LinearConstraints(self._linear_constraints.D_sites_modtypes, self._linear_constraints.specifications)
        
        L_econstraints                 = [] # initializing list of elastic constraints
        L_PuLP_econstraints            = []
        D_variables_modforms_ex        = {} 
        D_PuLP_variables_ex            = {} # initializing extended dictionary

        # initializing extended dictionary for Modform objects
        D_variables_modforms_ex.update(self._linear_constraints.modform_space.D_variables_modforms)  
        
        # initializing extended dictionary for LpVariable objects
        D_PuLP_variables_ex.update(self._linear_constraints.modform_space.D_PuLP_variables)
        
        for i in range(1,len(self._linear_constraints.L_linear_constraints)+1):
            
            lhs = self._linear_constraints.L_linear_constraints[i-1][0]
                        
            elastic_variable1 = "e_%s"%(2*i-1)
            elastic_variable2 =  "e_%s"%(2*i)
            evariables_diff = " + " + elastic_variable1 + " - " + elastic_variable2

            PuLP_evariable1 = pulp.LpVariable(elastic_variable1, lowBound=0)
            PuLP_evariable2 = pulp.LpVariable(elastic_variable2, lowBound=0)
            
            D_variables_modforms_ex[elastic_variable1] = modform.Modform(elastic_variable1)
            D_variables_modforms_ex[elastic_variable2] = modform.Modform(elastic_variable2)
                        
            D_PuLP_variables_ex[elastic_variable1] = PuLP_evariable1
            D_PuLP_variables_ex[elastic_variable2] = PuLP_evariable2
                        
            lhs_new = lhs + evariables_diff
            
            L_econstraints.append((lhs_new, self._linear_constraints.L_linear_constraints[i-1][1]))

            # 
            lhs_PuLP = sum([v[0]*v[1] for v in self._linear_constraints.L_PuLP_constraints[i-1].items()])
            LHS = lhs_PuLP + PuLP_evariable1 - PuLP_evariable2
            
            L_PuLP_econstraints.append(pulp.LpConstraint(LHS, rhs=float(self._linear_constraints.L_linear_constraints[i-1][1])))           
        
        E_linear_constraints.L_linear_constraints = L_econstraints
        E_linear_constraints.L_PuLP_constraints   = L_PuLP_econstraints
        
        E_linear_constraints.modform_space.D_variables_modforms = D_variables_modforms_ex
        E_linear_constraints.modform_space.D_PuLP_variables     = D_PuLP_variables_ex
        
        
        return E_linear_constraints
