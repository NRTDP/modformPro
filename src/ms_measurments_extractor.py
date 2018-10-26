from collections import defaultdict



class MSMeasurementsExtractor(object):
    """This class extracts measurement values (peak areas) from TDMS, BUMS and MRM data."""
    
    def __init__(self):
        
        self._L_isomeric_modforms_peak_areas = []
        self._D_peptide_isomeric_peak_areas = defaultdict(list)
        
    
    @property
    def L_isomeric_modforms_peak_areas(self):
        """Instance attribute
        
        * Python list of isomeric modforms and their percentage peak areas
        * This is populated by a member method
        
        >>> # Demonstrative example of L_isomeric_modforms_peak_areas
        [({'UNMOD': '0'}, '20'), ({'P': '1'}, '15'), ({'Ac': '1', 'P': '2'}, '25'), ({'P': '3'}, '10'), ({'Ac': '1', 'P': '4'}, '30')]
        """
        return self._L_isomeric_modforms_peak_areas

    @L_isomeric_modforms_peak_areas.setter
    def L_isomeric_modforms_peak_areas(self, L_isomeric_modforms_peak_areas):
        self._L_isomeric_modforms_peak_areas = L_isomeric_modforms_peak_areas

    @property
    def D_peptide_isomeric_peak_areas(self):
        """Instance attribute
        
        * Python dictionary having peptide mapped to its isomeris states and their respective percentage peak areas
        * This is populated by a member method
        
        >>> # Demonstrative example of D_peptide_isomeric_peak_areas
        D = {'15-25': [({'UNMOD': '0'}, '25'), ({'P': '1'}, '50'), ({'P': '2'}, '25')], '50-60': [({'UNMOD': '0'}, '50'), ({'Ac': '1', 'P': '1'}, '30'), ({'P': '2'}, '20')]}
        """
        return self._D_peptide_isomeric_peak_areas
    
    @D_peptide_isomeric_peak_areas.setter
    def D_peptide_isomeric_peak_areas(self, D_peptide_isomeric_peak_areas):
        self._D_peptide_isomeric_peak_areas = D_peptide_isomeric_peak_areas
    
    def get_list_of_isomeric_modforms_peak_areas(self, fout_isomeric_modforms_peak_areas):
        """Instance method
        
        * It populates the instance attribute, ``L_isomeric_modforms_peak_areas``
        * Note that the final list is sorted based on total frequency for isobaric modforms 
        
        >>> # Demonstrative example of  L_isomeric_modforms_peak_areas
        [({'UNMOD': '0'}, '20'), ({'P': '1'}, '15'), ({'Ac': '1', 'P': '2'}, '25'), ({'P': '3'}, '10'), ({'Ac': '1', 'P': '4'}, '30')]``

        Parameters
        ----------
        fout_isomeric_modforms_peak_areas : str
          File containing percentage peak areas of isomeric modforms. This data typically comes from TDMS measurments.
        
        Returns
        -------
        L_isomeric_modforms_peak_areas : list
          python list of isomeric modforms and their percentage areas. This data is typically known from intact MS measurements.
        """
        
        L_isomeric_modforms_peak_areas = []
        
        try:
            fin_handle = open(fout_isomeric_modforms_peak_areas)

        except IOError:
            raise("Provide a file containing percentage peak areas of isomeric modforms")

        block_start = False
        
        for line in fin_handle:
            
            line = line.strip()
        
            # skipping the blank line
            if not line.strip():
                continue
             
            # skipping the comment line
            if line[0]=="#":
                continue
            
            
            if   line=="START": block_start=True
            elif line=="END"  : block_start=False
            
            if block_start and line != "START":
                
                right_side = line.split(":")[1]
                right_side = right_side.strip()
                
                L_modtypes_freq_peak_area = [m.strip() for m in right_side.split(" ")]
                percentage_peak_area      = L_modtypes_freq_peak_area[-1] # last column
                
                D_modtype_freq = {}
                
                # running the loop so as to skip the last element
                for i, m in enumerate(L_modtypes_freq_peak_area[:-1]):
                    mtype = (m.split('=')[0]).strip()
                    freq  = (m.split('=')[1]).strip()
                    D_modtype_freq[mtype] = freq

                #end for

                L_isomeric_modforms_peak_areas.append((D_modtype_freq, percentage_peak_area))

            else: #if block_start and line != "START" and line[0]!="#"
                continue
        
        #end for line in fin_handle             
        
        # sorting the local list and storing it in a member attribute
        # sorting is based on sum of frequencies for isobaric modforms
        
        self._L_isomeric_modforms_peak_areas = sorted(L_isomeric_modforms_peak_areas, key=lambda x: sum([int(f) for f in x[0].values()]))
        

    
    def get_dictionary_of_peptides_and_isomeric_peak_areas(self, fout_peptides_isomeric_peak_areas):
        """Instance method
        
        * It populates the instance attribute, D_peptide_isomeric_peak_areas
        
        >>> # Demonstrative example of D_peptide_isomeric_peak_areas
        D = {'15-25': [({'UNMOD': '0'}, '25'), ({'P': '1'}, '50'), ({'P': '2'}, '25')], '50-60': [({'UNMOD': '0'}, '50'), ({'Ac': '1', 'P': '1'}, '30'), ({'P': '2'}, '20')]}

        Parameters
        ----------
        fout_peptides_isomeric_peak_areas : str
          File containing peptides, their isomeric states and respective percentage peak areas. This data typically comes from BUMS/MRM measurments.
        """
                
        try:
            fin_handle = open(fout_peptides_isomeric_peak_areas)

        except IOError:
            raise("Provide a file containing percentage peak_area of isomeric peptides")

        # local list; this appends all lines within a block and is emptied at the end of the block
        L_peptide_isomeric_peak_area = []

        block_start = False; pep_NC = ""
                
        for line in fin_handle:
            
            line = line.strip()
        
            # skipping the blank line
            if not line.strip():
                continue
            
            # skipping the comment line
            if line[0]=="#":
                continue
            
            
            if   line=="PEPSTART": block_start=True

            elif line=="PEPEND"  :
                block_start=False
                                
            #end elif
            
            if block_start and line!="PEPSTART":
                L = line.split(":")
                if L[0].strip() == "peptide":
                    pep_NC = L[1].strip() #e.g, '15-25'

                elif L[0].strip()=="IsomericPeptidesPeakArea":
                    right_side = L[1].strip()

                    L_modtypes_freq_peak_area = [m.strip() for m in right_side.split(" ")]
                    percentage_peak_area      = L_modtypes_freq_peak_area[-1] # last column
                    D_modtype_freq            = {}

                    # running the loop so as to skip the last element
                    for i, m in enumerate(L_modtypes_freq_peak_area[:-1]):
                        mtype = (m.split('=')[0]).strip()
                        freq  = (m.split('=')[1]).strip()
                        D_modtype_freq[mtype] = freq

                    #end for

                    L_peptide_isomeric_peak_area.append((D_modtype_freq, percentage_peak_area))
                    
            # end if block_start and line!="PEPSTART"    

            # pushing into the dictionary after end of each block

            if line=="PEPEND":

                # sorting the list based on total frequency of isomeric peptides
                L_sorted = sorted(L_peptide_isomeric_peak_area, key=lambda x: sum([int(f) for f in x[0].values()]))
                
                self._D_peptide_isomeric_peak_areas[pep_NC] = L_sorted
                
                #emptying the list for next block
                L_peptide_isomeric_peak_area = []
                                
                # emptying the peptide N_loc, C_loc string at the end of the block
                pep_NC = ""
                
    
