import pulp


class Modform(object):
    
    """This class defines modform, its variable name, representation and other methods to operate on the data-structures of this class.
    
    Demo representations:
    
    ``a1:  [('0', 'UNMOD')]``
    
    ``a11: [('T11', 'P'), ('S30', 'P'), ('S55', 'Ac')]``
    

    Parameters
    ----------
    var_name : str
      Modform variable name
    
    D_protein_sites_modifications : dict
      Python dictionary storing aa residues mapped to modifications on the residues.

    percentage_abundance : str or float
      Percentage abundance of this modform
    """
    
    
    def __init__(self, var_name, D_protein_sites_modifications=None, percentage_abundance=None):
        
        self._var_name              = var_name

        # modform variable in PuLP format, pulp.LpVariable instance
        self._PuLP_var              = pulp.LpVariable(var_name, 0, 100)
        
        # Dictionary of sites and associated modifications for this modform instance
        self._D_protein_sites_modifications = D_protein_sites_modifications
        
        self._percentage_abundance  = percentage_abundance

    
    def __repr__(self):
        
        """Instance method
        
        * printable representation of the object
        * E.g.,
        
        ``a1:  [('0', 'UNMOD')]``
        
        ``a11: [('T11', 'P'), ('S30', 'P'), ('S55', 'MM')]``
        """
        if self._D_protein_sites_modifications=={'0':'UNMOD'}:
            return (("%s: %s")%(self._var_name, self._D_protein_sites_modifications.items()))
        else:
            L_sorted = sorted(self._D_protein_sites_modifications.iteritems(), key=lambda x: int(x[0][1:]))
            
            return (("%s: %s")%(self._var_name, L_sorted))
        
    @property
    def var_name(self):
        """Instance attribute
        
        Variable name for the modform
        """
        return(self._var_name)
    
    @var_name.setter
    def var_name(self, var_name):
        self._var_name = var_name

    @property
    def PuLP_var(self):
        """Instance attribute
        
        * modform variable in PuLP format
        * E.g., ``LpVariable('a1')``
        """
        return  self._PuLP_var

    @PuLP_var.setter
    def PuLP_var(self, PuLP_var):
        self._PuLP_var = PuLP_var
    
    @property
    def D_protein_sites_modifications(self):
        """Instance attribute
        
        Python dictionary with aa residues mapped to a list of possible modifications on the residues
        """
        #L_sorted = sorted(self._D_protein_sites_modifications, key=lambda x: int(x[0]))
        return(self._D_protein_sites_modifications)

    @D_protein_sites_modifications.setter
    def D_protein_sites_modifications(self, D_protein_sites_modifications):
        self._D_protein_sites_modifications = D_protein_sites_modifications
    
    @property
    def percentage_abundance(self):
        """Instance attribute
        
        Percentage abundance of a modform
        """
        return(self._percentage_abundance)
    
    @percentage_abundance.setter
    def percentage_abundance(self, percentage_abundance):
        self._percentage_abundance = percentage_abundance
    
    
    def get_modifications_frequency_map(self):
        
        """Instance method
        
        * It returns a dictionary of modifications and their frequency in the modform OR 'UNMOD' in case of unmodified modform
        * E.g., 
          
          For ``a1: ['0':'UNMOD']``, the method returns ``{'UNMOD':'0'}``. 
          For ``a11: [('S11','P'),('S30','P'),('S55','Ac')]``, the method returns, ``{'P': '2', 'Ac': '1'}``
        """

        if self._D_protein_sites_modifications.items()==[('0','UNMOD')]:
            return {'UNMOD':'0'}
        
        else:
                        
            D_protein_mod_frequency = {}
            
            for site, mod in self._D_protein_sites_modifications.iteritems():
                if mod in D_protein_mod_frequency.keys():
                    freq =  int(D_protein_mod_frequency[mod])
                    D_protein_mod_frequency[mod] = str(freq+1)
                    
                else:
                    D_protein_mod_frequency[mod] = '1'
                    
            return D_protein_mod_frequency

    
    def whether_modified_peptide_present_in_modform(self, N_loc, C_loc, D_peptide_mod_frequency=None):
        """Instance method
        
        * Boolean method return True or False
        * peptide range is defined by N_loc, C_loc
        * modification state of peptide: D_peptide_mod_frequency
        * UNMOD modform can only generate UNMOD peptides
        * UNMOD peptides could also come from other modforms too
        
        Parameters
        ----------
        N_loc : str or int
          Position of the first residue of a peptide in the protein sequence
        
        C_loc : str or int
          Position of the last residue of a peptide in the protein sequence
        """
        
        if self._D_protein_sites_modifications.items()==[('0','UNMOD')]:
            if D_peptide_mod_frequency.items()==[('UNMOD','0')]: return True
            else: return False
            #end if-else
                
        else:
            L_mod_sorted = sorted(self._D_protein_sites_modifications.iteritems(), key=lambda x: int(x[0][1:]))
            
            #list of modifications in the peptide in this modform
            peptide_modifications = []
            
            for pair in L_mod_sorted:
                
                if int(pair[0][1:]) < int(N_loc) or int(pair[0][1:]) > int(C_loc):
                    continue
                
                else:
                    peptide_modifications.append(pair[1])
                #end if-else
            #end for pair in L_mod_sorted
            
            total_freq = sum([int(a) for a in D_peptide_mod_frequency.values()])            
            if D_peptide_mod_frequency.items()==[('UNMOD','0')] or total_freq==0:
                #print "N_loc, C_loc, D_pep_mod_freq, mform: ", N_loc, C_loc, D_peptide_mod_frequency, self._D_protein_sites_modifications
                if len(peptide_modifications)==0:
                    return True
                else: # len(peptide_modifications)>0
                    return False
                #end if-else len(peptide_modifications)==0
            else:
                # list of True,False upon cross-verification
                L_bin = [peptide_modifications.count(k)==int(v) for k,v in D_peptide_mod_frequency.items()]
                        
                if all(L_bin):
                    return True
                else:
                    return False
                #end if-else all(L_bin)
            #end if-else D_peptide_mod_frequency.items()==[('UNMOD','0')] or total_freq==0
        #end if-else self._D_protein_sites_modifications.items()==[('0','UNMOD')]
