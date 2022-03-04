from headers.StandardModules import *
from headers.DefaultFilePaths import *

class LoadData:
    
    def __init__(self):
        
        # Self Capacitance Experimental Data
        self.MedhurstData = np.genfromtxt(path_to_data+"Medhurst_Self_Capacitance_Experimental_Data.dat",
                                          skip_header=True,
                                          delimiter=',',
                                          unpack='False')                                                              # pF/cm, cm/cm 

        # 4NEC2 Data (Resistance/Reactance) from 1.5 - 80 MHz
        self.InputResistanceNEC = np.genfromtxt(path_to_data+"InputResistanceNEC.dat",
                                          skip_header=True,
                                          delimiter='	',
                                          unpack='False')                                                              # MHz, ohms         

        self.InputReactanceNEC = np.genfromtxt(path_to_data+"InputReactanceNEC.dat",
                                          skip_header=True,
                                          delimiter='	',
                                          unpack='False')                                                              # MHz, ohms         

        # 4NEC2 Data (Resistance/Reactance) at Self Resonant Frequency
        self.InputResistanceNEC_AtResonance = np.genfromtxt(path_to_data+"InputResistanceNEC_AtResonance.dat",
                                                            skip_header=True,
                                                            delimiter='	',
                                                            unpack='False')                                            # MHz, ohms  

        self.InputReactanceNEC_AtResonance = np.genfromtxt(path_to_data+"InputReactanceNEC_AtResonance.dat",
                                                           skip_header=True,
                                                           delimiter='	',
                                                           unpack='False')                                             # MHz, ohms                  
        
        # 4NEC2 Resonanting Capacitance
        self.ResonatingCapacitanceNEC = np.genfromtxt(path_to_data+"ResonatingCapacitance_NEC2.dat",
                                                      skip_header=True,
                                                      delimiter='	',
                                                      unpack='False')                                                  # MHz, pF 
        
        # 4NEC2 Phase
        self.PhaseNEC = np.genfromtxt(path_to_data+"Phase_NEC2.dat",
                                      skip_header=True,
                                      delimiter='	',
                                      unpack='False')                                                                  # MHz, Degrees 
        
        # EZNEC 2++ Pro Version 7.0 Data (Resistance/Reactance) from 1.5 - 80 MHz
        self.InputResistanceReactance_EZNEC = np.genfromtxt(path_to_data+"InputResistanceReactance_EZNEC.dat",
                                                            skip_header=True,
                                                            delimiter=',',
                                                            unpack='False')                                            # MHz, ohms, ohms

        # NanoVNA Measurements from 1 - 90 MHz
        self.NanoVNA = np.genfromtxt(path_to_data+"1MHz_to_90MHz_Combined.s2p",
                                     skip_header=True,
                                     delimiter=' ',
                                     unpack='False')