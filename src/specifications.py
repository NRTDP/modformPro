class ModformProSpecifications(object):
    """This class contains specifications for modformPro package.
    
    """
    def __init__(self):
        self.file_fasta_sequence = ""
        self.file_modifications  = ""
        
        self.data_mode           = ""
                
        self.file_tdms_data      = ""
        self.file_bums_mrm_data  = ""

        self.variables_grouping  = True
        
        self.number_modes        = 0
        self.modes               = []
        self.noise_stdev         = 2
        self.sample_size         = 30
        
        self.Trypsin             = True
        self.GluC                = False
        
        self.n_sites_full_modform_space = True
        
        self.output_dir          = 'Output'
        self.tex_output_file     = ""
        
    def read_specifications_from_file(self, fname_spec):
        """Instance method
        
        Parameters
        ----------
        fname_spec : str
          file name containing specification fields.
        
        Notes
        -----
        Parses the file line by line and stores the value of the fields in the respective member attributes.
        """
        fhandle = open(fname_spec)
        
        for line in fhandle:
            line = line.strip()
            
            # skipping the blank line
            if not line.strip():
                continue
            
            # skipping the first line
            if line[0]=="#":
                continue
            
            info = line.split(':')
            var  = info[0].strip()

            if var=='FILE_MODIFICATIONS':
                self.file_modifications = info[1].strip()
                        
            elif var=='FILE_FASTA':
                self.file_fasta_sequence = info[1].strip()
            
            elif var=='DATA_MODE':
                self.data_mode = info[1].strip().lower()
                
            
            elif var=='FILE_TDMS_DATA':
                self.file_tdms_data = info[1].strip()
                                
            elif var=='FILE_BUMS_MRM_DATA':
                self.file_bums_mrm_data= info[1].strip()
                            
            elif var=='VARIABLES_GROUPING':
                user_cmd = info[1].strip()
                if user_cmd.lower()=='true': self.variables_grouping=True
                else: self.variables_grouping = False
                
                
            elif var=='NUMBER_MODES':
                nmodes = info[1].strip()
                self.number_modes = int(nmodes)
            
            elif var=='MODES':
                modes_string = info[1].strip()
                
                # setting instance variables if modes_string is not empty
                if modes_string:
                    
                    L_modes   = [v.strip() for v in modes_string.split(';')]
                    for mode in L_modes:
                        L_modifications_freq = mode.split(',')
                        D_mod_freq = {}
                        for mod_freq in L_modifications_freq:
                            mod,freq = mod_freq.split('=')
                            mod = mod.strip()
                            freq = freq.strip()
                            D_mod_freq[mod]=freq
                            
                        #end for mod_freq in L_modifications_freq

                        #setting the instance variable self._modes
                        self.modes.append(D_mod_freq)

                    #end for mode in L_modes
                
                #end if not modes_string
                
            elif var=='NOISE_STDEV':
                stdev = info[1].strip()
                self.noise_stdev = stdev

            elif var=='SAMPLE_SIZE':
                sample_size = info[1].strip()
                self.sample_size = sample_size
            
            elif var=='Trypsin':
                user_cmd = info[1].strip()
                if user_cmd.lower()=='true': self.Trypsin=True
                else: self.Trypsin = False
            
            elif var =='GluC':
                user_cmd = info[1].strip()
                if user_cmd.lower()=='true': self.GluC=True
                else: self.GluC = False
            
            elif var=='N_SITES_FULL_MODFORM_SPACE':
                user_cmd = info[1].strip()
                if user_cmd.lower()=='true': self.n_sites_full_modform_space=True
                else: self.n_sites_full_modform_space = False
                            
            elif var=='OUTPUT_DIRECTORY':
                self.output_dir = info[1].strip()
            
            elif var=='TEX_OUTPUT_FILE':
                self.tex_output_file = info[1].strip()
            
            else:
                continue
            
            #end if-elif-else var=='NUMBER_MODES'

