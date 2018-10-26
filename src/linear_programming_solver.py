import pulp

class LinearProgrammingSolver(object):
    
    """This class contains attributes and methods to solve linear programming problems.
    
    Parameters
    ----------
    linear_constraints : :class:`.LinearConstraints`
      Set of linear constraints.
    """
    
    def __init__(self, linear_constraints):

        # instance of the class LinearConstraints
        self._linear_constraints = linear_constraints
    
    @property
    def linear_constraints(self):
        """Instance attribute
        
        This is an instance of :class:`.LinearConstraints` object.
        """
        return(self._linear_constraints)
    
    
    # using puLP
    def linear_programming_solver(self, objective_function_PuLP, cmd):
        """Instance method
        
        * This method solves a linear program given list of LP constraints, objective function and command (LpMaximize or LpMinimize). 
        * It returns LpProblem object containing variables, constraints and solution to the linear problem.
        * Note that LP constraints and objective function must be according to :class:`.puLP` format.
          - LP constraint string: ``lhs = rhs``
          - PuLP constraint:  ``pulp.LpConstraint(lhs, rhs)``
          - Objective function should in terms of LP variables in puLP format, e.g., transforming a variable 'v' -- ``pulp.LpVariable(v, lowBound=0)`` 

        Parameters
        ----------
        objective_function_PuLP : Affine expression of :class:`.LpVariable` objects
            Objective function for a linear programming problem
          
        cmd : LpMaximize (or -1); LpMinimize (or 1)
          Commands maximization or minimization of the objective funtion in a linear programming problem
        """
        
        #L_PuLP_constraints   = self._linear_constraints.get_list_of_PuLP_constraints()
                
        LP_problem = pulp.LpProblem("problem",cmd)

        # objective function
        LP_problem += objective_function_PuLP
        
        for constraint in self._linear_constraints.L_PuLP_constraints:
            LP_problem += constraint

        #print "LP_problem: ", LP_problem        
        # solve the problem
        status = LP_problem.solve(pulp.solvers.PULP_CBC_CMD(msg=0))
        

        
        #print "status: ", LpStatus[status]
        #for v in prob.variables():
        #    print v, pulp.value(v)
        #print "objective value=", pulp.value(prob.objective)
        #print "\n"
        
        return LP_problem
