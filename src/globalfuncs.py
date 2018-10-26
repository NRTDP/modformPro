from collections import defaultdict
from itertools import product, combinations, chain

import pulp
import re
import os

import modform as modform


class TexFileGenerator(object):
   """This class prints the output content in a tex file.
   
   Parameters
   ----------
   ftex_out : str
     File name with .tex extension
   
   output_dir : str
     Name of the output directory
   
   """
   def __init__(self, ftex_out, output_dir):
      
      # file name with .tex extension
      # e.g., foo.tex
      self._ftex_out = ftex_out
      
      self._output_dir = output_dir
      if not os.path.exists(output_dir):
            os.makedirs(output_dir)
      
      # add header to the file
      self.get_a_blank_tex_file_with_header()

   
   def print_content_in_tex_file(self, protein_seq, D_sites_modtypes, D_lower_bounds, D_upper_bounds, mform_region_plot_pdf, D_variables_modforms, D_grouped_variables=None):
      """Instance method
      
      Prints output content in a tex file.
      
      Parameters
      ----------
      protein_seq : str
        Amino-acid sequence of the protein under investigation.
      
      D_sites_modtypes : defaultdict(list)
        Python dictionary containing residues mapped to the list of possible modifications. E.g., ``D = {'T32': ['P'], 'S55': ['Ac', 'P']}``

      D_lower_bounds : dict
        Python dictionary mapping modform variables to their lower bouds of percentage concentration.

      D_upper_bounds : dict
        Python dictionary mapping modform variables to their upper bouds of percentage concentration.

      mform_region_plot_pdf : str
        Name of the file with a pdf extension containing modform region plot.
      
      D_variables_modforms : dict
        Python dictionary mapping modform variables (str) to their :class:`.Modform` objects.
      
      D_grouped_variables : dict, optional
        Python dictionary mapping new coalesced variable to sum of modform concentration variables.
      """
      # printing the formatted fasta seq in a tex file
      self.add_formatted_fasta_seq_in_the_file(protein_seq, D_sites_modtypes)
      newline_string = '\n\\bigskip\n'

      self.add_sites_modifications_table_in_tex_file(D_sites_modtypes)
      self.add_modform_region_plot_in_tex_file(mform_region_plot_pdf)
      if len(D_grouped_variables)!=0:
         self.add_grouped_variables_in_tex_file(D_grouped_variables)
      
      self.add_percentage_concentration_bounds_in_tex_file(D_lower_bounds, D_upper_bounds)
      
      self.add_modforms_in_tex_file(D_variables_modforms)
      
      #adding the end document directive in the tex file
      self.add_end_document_statement()
   
   def add_sites_modifications_table_in_tex_file(self, D_sites_modtypes):
      """Instance method
      
      Adds sites and modifications table in tex file.
     
      Parameters
      ----------
      D_sites_modtypes : defaultdict(list)
        Python dictionary mapping aa residues to a list of modification types to the respective residues.
      """
      table_string = '\\medskip \n'
      table_string += '\n\n\\begin{table}[h]\n\\centering\n\\begin{tabular}{|c|c|}\n \\hline \n'
      table_string += '\\textbf{Residue number} & \\textbf{possible modifications}\\\\ \\hline \n'
      for site,mtypes in sorted(D_sites_modtypes.iteritems(), key=lambda x: int(x[0][1:])):
         type_string = '%s'%mtypes[0]
         for i in range(1, len(mtypes)):
            type_string += ', %s'%mtypes[i]
         row_string = '%s & %s\\\\ \n'%(site, type_string)
         row_string += '\\hline \n'
         table_string += row_string
        
      #end for site,mtypes in D_sites_modtypes.iteritems()
        
      table_string += '\\end{tabular}\n\\caption{Residue number and possible modifications on them}\\end{table}\n\\medskip\n'
      self.add_latex_code_in_the_file(table_string)


   def add_modform_region_plot_in_tex_file(self, mform_region_plot_pdf):
      """Instance method
      
      Add modform region plot in the tex file.
      
      Parameters
      ----------
      mform_region_plot_pdf : str
        Name of the file with a pdf extension containing modform region plot.
      """
      # adding modform region plot in a tex file
      #plot_string  = '\n{\large Figure: Modform region plot}\n'
      plot_string  = '\n \\begin{figure}[h]\n\\centering\n'
      plot_string += '\\includegraphics[height=7.5cm]{%s}\n'%mform_region_plot_pdf
      plot_string +='\\caption{Modform region plot}\n\\end{figure}\n'
      self.add_latex_code_in_the_file(plot_string)

   
   def add_grouped_variables_in_tex_file(self, D_grouped_variables):
      """Instance method
      
      Add coalesced variables in the tex file.
      
      Parameters
      ----------
      D_grouped_variables : dict
        Python dictionary mapping new coalesced variable to a sum of variables.
      """
      grouped_variables_string = '\n \\newpage \n \\noindent \n \\textbf{\\large Coalesced modform concentration variables:}\\newline\n'
      for gv,gp in sorted(D_grouped_variables.iteritems(), key=lambda x:natural_keys(x[0])):
         grouped_variables_string += '%s: %s \\newline \n'%(gv,gp)
      #end for gv,gp in sorted(linear_constraints.D_grouped_variables.iteritems(), key=lambda x:gbfunc.natural_keys(x[0]))
      self.add_latex_code_in_the_file(grouped_variables_string)

   
   def add_percentage_concentration_bounds_in_tex_file(self, D_lower_bounds, D_upper_bounds):
      """Instance method
      
      Adds percentage concentration bounds of modforms in a tex file.
     
      Parameters
      ----------
      D_lower_bounds : dict
        Python dictionary mapping modform variables to their lower bouds of percentage concentration.
      
      D_upper_bounds : dict
        Python dictionary mapping modform variables to their upper bouds of percentage concentration.
      """
      percentage_bounds_string = '\n \\newpage \n \\noindent \n \\textbf{\\large Percentage bounds for modforms:}\\newline \n'
      
      for k,v in sorted(D_lower_bounds.iteritems(), key=lambda x:natural_keys(x[0])):
            percentage_bounds_string += "$%s \leq %s \leq %s$\\newline\n"%(D_lower_bounds[k],k,D_upper_bounds[k])

      self.add_latex_code_in_the_file(percentage_bounds_string)
      
   def add_modforms_in_tex_file(self, D_variables_modforms):
      """Instance method

      Adds modform variables and their modification maps in the tex file.
      
      Parameters
      ----------
      D_variables_modforms : dict
        Python dictionary mapping modform variables (str) to their :class:`.Modform` objects.
      """
      mform_tex_string = '\n \\newpage \n \\noindent \n \\textbf{\\large Modform variables and their modifications map} \n\\newline\n'
      #self._tex_output_file.add_latex_code_in_the_file(mform_header_string)
        
      for mvar,mform in sorted(D_variables_modforms.iteritems(), key=lambda x: natural_keys(x[0])):
         mform_string = '%s\n'%(mform)
         mform_tex_string += '%s\\newline\n'%mform_string
            
      #end for mvar,mform in sorted(linear_constraints.D_variables_modforms.iteritems(), key=lambda x: natural_keys(x[0]))
      self.add_latex_code_in_the_file(mform_tex_string)

   
   def get_a_blank_tex_file_with_header(self):
      ftex_out = '%s/%s'%(self._output_dir, self._ftex_out)
      fout = open(ftex_out,'w')
      
      print_string = '\\documentclass{article}\n\\usepackage[margin=0.5in]{geometry}\\usepackage{graphicx}\n\\usepackage[usenames, dvipsnames]{color}\n\\usepackage{float}\n'
      print_string += '\n\\begin{document}\n'
      fout.write(print_string)
      fout.close()
      

   def add_latex_code_in_the_file(self, tex_string):
      """Instance method
      
      Adds latex code in the tex file.
      
      Parameters
      ----------
      tex_string : str
        Latex code string to be included.
      """
      ftex_out = '%s/%s'%(self._output_dir, self._ftex_out)
      fout = open(ftex_out,'a')
      fout.write(tex_string)
      fout.close()
   
   def add_formatted_fasta_seq_in_the_file(self, protein_seq, D_sites_modtypes):
      """Instance method
      
      This function adds the fasta seq string in a particular format and also highlights modification sites by a bold blue and underlined.
      """
    
      # list of possible modification sites
      L_mod_sites = D_sites_modtypes.keys()
      
      L_amino_acids = list(protein_seq)

      #making every modification site bold, green, underlined
      for msite in L_mod_sites:
         index = int(msite[1:])-1
         L_amino_acids[index] = "\\underline{\\textbf{\\textcolor{OliveGreen}{\\Large %s}}}"%L_amino_acids[index]

      #end for msite in L_mod_sites

      #code for printing the amino acid sequence in blocks of 10 amino acids
    
      #making lists of 10 amino acid sequence
      L_lists_of_10_aa = [L_amino_acids[i:i+10] for i in range(0,len(L_amino_acids),10)]

      #joining each sublist of 10 amino-acids amino acids as one string after making relevant sites bold, green and underlined
      L_seq_blocks = []
      for sublist in L_lists_of_10_aa:
         L_seq_blocks.append(''.join(sublist))

   
      # printing 5 blocks per line, each block has 10 amino-acids
      formatted_fasta_seq = ''

      for bindex,block in enumerate(L_seq_blocks):
      
         block_count = bindex + 1
      
         #At the beginning of every line
         if block_count % 5 == 1:
            aa_num = (block_count-1)*10 + 1
            formatted_fasta_seq += '%s & %s'%(aa_num,block)
               
         #At the end of every line
         elif block_count % 5 == 0:
            aa_num = block_count*10
            formatted_fasta_seq += "& %s &  %s \\\\ \n"%(block,aa_num)
               
         #For intermediate blocks (Otherwise)
         else:
            formatted_fasta_seq += "& %s"%block
      
         #end if-else
      # end for bindex,block in enumerate(L_seq_blocks)

      # adding the formatted_fasta_seq in the tex file
      ftex_out = '%s/%s'%(self._output_dir, self._ftex_out)
      fout = open(ftex_out, 'a')
      print_string = '\\noindent\n'
      print_string += '{\\large Amino acid sequence highlighting modification residue}\n'
      print_string += '\\begin{table}[h]'
      #print_string += '\\hspace{-2cm}\n'
      print_string += '\\begin{tabular}{lllllll}\n'
      fout.write(print_string)
      fout.write(formatted_fasta_seq)
      fout.write('\n\\end{tabular}\n')
      fout.write('\\end{table}\n')
      fout.close()

      #return formatted_fasta_seq 
      
   
   def add_end_document_statement(self):
      ftex_out = '%s/%s'%(self._output_dir, self._ftex_out)
      fout = open(ftex_out,'a')
      print_string = '\n\\end{document}'
      fout.write(print_string)
      fout.close()

   

      


def read_sites_and_modification_types_from_file(fname_modifications):
   """Global method
   
   * Function populated instance variable, a dictionary, D_sites_modtypes with sites associated to their modfication types
   * Read from an input file
   * This information must be known apriori
   
   Parameters
   ----------
   fname_modifications : str
     Name of an input file containing amino-acid residues and their respective modification types.
   
   Returns
   -------
   D_sites_modtypes : defaultdict(list)
     Python dictionary mapping aa residues to a list of modification types to the respective residues. 
   E.g., ``{'13':['P']; '16':['P']; '19':['P']; '40':['P', 'Ac']}``
   """

   # initializing variable for number of sites
   D_sites_modtypes = defaultdict(list)
    
   nsites = 0

   #opening a file in read mode
   fin_handle = open(fname_modifications)

   # Binary variables to detect tags in the file are initialized to be False
   read_nsites = False; read_modification_sites = False
   
   for line in fin_handle:
      
      line = line.strip()
      
      # skipping the blank line
      if not line.strip():
         continue
      # skipping the comment line
      if line[0]=="#":
         continue
      
      if   line=="NSITESBEGIN": read_nsites = True
      elif line=="NSITESEND"  : read_nsites = False
      
      if read_nsites and line !="NSITESBEGIN":
         nsites = line

      if   line=="MODSITESBEGIN": read_modification_sites = True
      elif line=="MODSITESEND"  : read_modification_sites = False

      if read_modification_sites and line!="MODSITESBEGIN":
         site_modtype = line.split(':')

         site    = site_modtype[0].strip()
         modtype = site_modtype[1].strip()
         
         D_sites_modtypes[site].append(modtype)

   #end for line in fin_handle
   
   if  int(nsites)!=len(D_sites_modtypes):
      raise Exception("Number of sites not identical to the size of dictionary of sites associated with their modification types")

   
   return D_sites_modtypes


def get_fasta_seq_from_the_file(file_fasta_sequence):
   '''
   The function reads fasta sequence from the file
   '''

   try:
      fin = open(file_fasta_sequence)

   except:
      raise IOError("provide fasta sequence file name in the specification file")
   
   
   seq = ''
   for line in fin:
      line = line.strip()
      
      #skip the blank line and commented lines
      if not line.strip() or line[0]=="#" or line[0]==">": continue
      
      #skip the commented line
      
      seq += line

   return seq



def get_modforms_from_simulated_data(D_variables_modforms_simulated=None, D_PuLP_variables_simulated=None, D_protein_mod_frequency=None):
   """Global method
   
   * The function is alternative to the function get_modforms.
   * The difference however is that the modform variables and objects are used from simulated data and not generated from the dictionary D_sites_modtypes
   * This helps in comparing the results of computed modform distribution to that of the simulated modform distribution.
   
   Parameters
   ----------
   D_variables_modforms_simulated : dict
     Python dictionary containing modform variables mapped to their :class:`.Modform` objects for simulated data.
   
   D_PuLP_variables_simulated : dict
     Python dictionary containing modform variables mapped to their :class:`.pulp` objects for simulated data.
   
   D_protein_mod_frequency : dict
     Python dictionary containing modifications mapped to their frequencues.
   
   Returns
   -------
   D_vars_modforms_isobaric : dict
     Python dictionary contaitning isomeric modform variables mapped to their :class:`.Modform` objects.  
   
   D_PuLP_variables_isobaric : dict
     Python dictionary containing isomeric modform variables mapped to their :class:`.pulp` objects.

   Notes
   -----
   """
   D_vars_modforms_isobaric = {}; D_PuLP_variables_isobaric = {}
   

   for var,mform in D_variables_modforms_simulated.iteritems():
      D_mform_mod_frequency = mform.get_modifications_frequency_map()

      if D_mform_mod_frequency == D_protein_mod_frequency:
         D_vars_modforms_isobaric[var]   = mform
         D_PuLP_variables_isobaric[var] = D_PuLP_variables_simulated[var] 
      else:
         continue

      #end if-else D_mform_mod_frequency == D_protein_mod_frequency

   #end for var,mform in D_variables_modforms_simulated.iteritems()

   return (D_vars_modforms_isobaric, D_PuLP_variables_isobaric)



    
def get_modforms(start_index, D_sites_modtypes=None, D_protein_mod_frequency=None):
   """Global method
   
   Provides a set of intact isomeric modforms having same modification frequency and types.
   
   Parameters
   ----------
   start_index : int
     Index of the modform variable name.
   
   D_sites_modtypes : defaultdict(list), optional
     Python dictionary containing residues mapped to the list of possible modifications. E.g., ``D = {'T32': ['P'], 'S55': ['Ac', 'P']}``
   D_protein_mod_frequency : dict
     Python dictionary containing modifications mapped to their frequencies for isomeric modforms. E.g., ``D = {'P': '2', 'Ac': '1'}``
   
   Returns
   -------
   D_vars_modforms_isobaric : dict
     Python dictionary contaitning isomeric modform variables mapped to their :class:`.Modform` objects.  
   
   D_PuLP_variables_isobaric : dict
     Python dictionary containing isomeric modform variables mapped to their :class:`.pulp` objects.
   
   Examples
   --------
   >>> 
   >>> D_sites_modtypes = {'15':['P']; '16':['P']; '52':['Ac', 'P']; '54':['P']; '59':['P']}
   >>> 
   >>> get_modforms(start_index=1, D_sites_modtypes)
   {'0': 'UNMOD'}
   >>>
   >>> D_protein_mod_frequency = {'P': '2'}
   >>> get_modforms(start_index=1, D_sites_modtypes, D_protein_mod_frequency)
   """
   
   D_vars_modforms_isobaric    = {}
   D_PuLP_variables_isobaric   = {}
   
   # for unmodified modform
   total_freq = sum([int(a) for a in D_protein_mod_frequency.values()])
   if D_protein_mod_frequency==None or D_protein_mod_frequency.items()==[('UNMOD','0')] or total_freq==0:
      
      modform_unmod    = modform.Modform("a1", D_protein_sites_modifications={'0':'UNMOD'})
      PuLP_modform_var = pulp.LpVariable('a1', lowBound=0)
      
      D_vars_modforms_isobaric['a1']   = modform_unmod
      D_PuLP_variables_isobaric['a1']  = PuLP_modform_var
      
      return (D_vars_modforms_isobaric, D_PuLP_variables_isobaric)
   
   else:
      # creating list of all possible dictionaries of isobaric occupancy
      L_dictionaries = [{j: D_sites_modtypes[j] for j in i} for i in combinations(D_sites_modtypes, sum([int(v) for v in D_protein_mod_frequency.values()]))]
      
      
      count = start_index #modform index starting at 2 (by default) for modified modforms
      
      for k in L_dictionaries:
         
         # list of modforms a sub-dictionary
         L_modforms = [dict(zip(k,v)) for v in product(*k.values())]
         
         for mform in L_modforms:
            
            # set of modifications
            set_modification_types = set(mform.values())
            
            # The type of modifications must be identical
            set_modification_types_sought = set(D_protein_mod_frequency.keys())
            if len(set_modification_types.symmetric_difference(set_modification_types_sought))>0:
               continue
            
            else:
               #list of modifications in the mform
               l_mod = mform.values()
               
               # list of True,False upon cross-verification
               L_bin = [l_mod.count(k)==int(v) for k,v in D_protein_mod_frequency.items()]
               if False in L_bin:
                  continue
               
               else:
                  
                  var_name = "a%s"%count #variable name
                  
                  # creating a modform object
                  modform_object = modform.Modform(var_name, D_protein_sites_modifications=mform)
                  PuLP_modform_var = pulp.LpVariable(var_name, lowBound=0)
                  
                  # populating the dictionary of isobaric modforms
                  D_vars_modforms_isobaric[var_name]    = modform_object
                  D_PuLP_variables_isobaric[var_name] = PuLP_modform_var
                  
                  count += 1

               #end if-else False in L_bin

            #end if-else len(set_modification_types.symmetric_difference(set_modification_types_sought))>0

         #end for mform in L_modforms

      #end for k in L_dictionaries
         
   #end if-else D_protein_mod_frequency==None or D_protein_mod_frequency.items()==[('UNMOD','0')]

   return (D_vars_modforms_isobaric, D_PuLP_variables_isobaric)

def natural_keys(text):
    """Global method
    
    Examples
    --------
    >>> L = ['a11', 'a9', 'a2']
    >>> print sorted(L)
    ['a11', 'a2', 'a9']
    >>>
    >>> print sorted(L, key=lambda x: natural_keys(x))
    ['a2', 'a9', 'a11']
    """
    return [int(c) if c.isdigit() else c for c in re.split('(\d+)', text)]
 
def powerset(iterable):
   """Global method
   
   Returns a chain of elements of powerset of the iterable.
   
   Examples
   --------
   >>> print list(powerset([1,2,3]))
   [() (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)]
   """
   s = list(iterable)
   return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))
