import os
import feasible_linear_constraints as flc
import linear_programming_solver as lps
import globalfuncs as gbfunc
import pulp

import re

class ModformRegionEstimator(object):
    
    """Computes modform region bounds given a set of 
    linear constraints.
    
    Parameters
    ----------
    linear_constraints : :class:`.LinearConstraints`
      Set of linear constraints.
    """
    
    def __init__(self, linear_constraints):
        
        self._linear_constraints = linear_constraints
        
    
    @property
    def linear_constraints(self):
        return(self._linear_constraints)
    
    
    def compute_bounds(self):
        """Instance method
        
        Compute bounds for each modform variable by maximizing and minimizing the variables in the linear program.
        
        Notes
        -----
        * Before getting feasible linear constraints, force all the modform variable in a liner constraint to be zero if rhs=0
        
        E.g., ``a13 + a16 + a18 = 0`` => ``a13=a16=a18=0``; ``set upBound=0``
        After doing that these variables wont play any role in fixing the infeasibility and finding the bounds subsequently
        
        * Create :class:`.LinearProgrammingSolver` object for feasible linear constraints
        * For each variable in the canonical ordered list run linear programming to minimize and maximize each modform variable
        """
        # Note: linear_constraints object has been been populated at this stage
        L_zero_var = []
        
        for constraint in self._linear_constraints.L_linear_constraints:
            lhs_string = constraint[0]
            rhs_string = constraint[1]
            if float(rhs_string)==0:
                #print "rhs=0: forcing the variables to zero"
                L_vars = re.split(r'[+-]',lhs_string)
                
                for var in L_vars:
                    modform_var = var.strip()
                    
                    # forcing all the variables in this constraint to be zero
                    self._linear_constraints.modform_space.D_PuLP_variables[modform_var] = pulp.LpVariable(modform_var, lowBound=0, upBound=0)
                    #print "var forced to zero: ", modform_var
                    L_zero_var.append(modform_var)
            else: #if float(rhs)==0
                continue
            
        if len(L_zero_var)>0:
            print "\n####### Variables forced to zero (rhs = 0) ##########"
            print "variables forced to zero: ", set(L_zero_var)
        
        feasible_lc = flc.FeasibleLinearConstraints(self._linear_constraints)
        
        feasible_lc.get_feasible_linear_constraints()
        
        feasible_linear_constraints = feasible_lc.feasible_linear_constraints
        
        lp_solver = lps.LinearProgrammingSolver(feasible_linear_constraints)
        
        D_lower_bounds = {}; D_upper_bounds = {}
        
        for v in [self._linear_constraints.modform_space.D_PuLP_variables[k] for k in sorted(self._linear_constraints.modform_space.D_PuLP_variables.keys(), key=gbfunc.natural_keys)]:
            
            if str(v) in L_zero_var:
                D_lower_bounds[str(v)] = '0'
                D_upper_bounds[str(v)] = '0'
                continue
            #end if str(v) in L_zero_var
            
            objective_function_PuLP = v
            
            list_values_minimize = lp_solver.linear_programming_solver(objective_function_PuLP, pulp.LpMinimize)
            D_lower_bounds[str(v)] = "%.3f"%round(pulp.value(v),3)
            
            list_values_maximize = lp_solver.linear_programming_solver(objective_function_PuLP, pulp.LpMaximize)
            D_upper_bounds[str(v)] = "%.3f"%round(pulp.value(v),3)

        #end for v in ..

        return((D_lower_bounds, D_upper_bounds))
