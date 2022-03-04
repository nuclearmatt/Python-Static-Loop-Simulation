from headers.StandardModules import *

class LocalConstants:
    
    def __init__(self):
        
        self.PermeabilityVacuum = 4*np.pi*10**-7                                                    # H/m (Vacuum)
        self.PermeativityVacuum = 1/(self.PermeabilityVacuum*sp.constants.speed_of_light**2)        # F/m (Vacuum)
        self.PermeabilityCopper = 1.256629*10**-6                                                   # H/m (Copper)
        self.CopperConductivity = 5.80*10**7                                                        # S/m (Copper)
        self.PermeabilityAir = 1.25663753*10**-6                                                    # H/m (Air)