from headers.StandardModules import *
from headers.DefaultFilePaths import *
from methods.ConvertUnits import ConvertUnits

class LoopAntenna:
    
    def __init__(self, 
                 ConductorDiameter,
                 AntennaCircumferance,
                 WidthCorrection,
                 NumberTurns,
                 UpperOperatingFrequency,
                 LowerOperatingFrequency,
                 IntendedOperatingFrequency,
                 Model,
                 LocalConstants,
                 Data,
                 **kwargs):
        
        self.LocalConstants = LocalConstants
        self.Data = Data
        self.ConductorDiameter = ConvertUnits(Attribute=ConductorDiameter,
                                              ConvertFrom='Inches',
                                              ConvertTo='Meters')                                    # m
        self.ConductorRadius = 0.5*self.ConductorDiameter                                            # m
        self.NumberTurns = NumberTurns 
        self.AntennaCircumferance = ConvertUnits(Attribute=AntennaCircumferance,   
                                                 ConvertFrom='Feet',
                                                 ConvertTo='Meters')                                 # m        
        self.Model = Model
        self.IntendedOperatingFrequency = IntendedOperatingFrequency
        self.DerivedRadius = self.AntennaCircumferance/(2*np.pi)

        # This calculation assumes equal sides
        self.WidthCorrection = ConvertUnits(Attribute=WidthCorrection,
                                            ConvertFrom='Feet',
                                            ConvertTo='Meters')                                     # m
        self.NumberSides = 8
        self.AdjustedAntennaCircumferance = self.AntennaCircumferance - self.WidthCorrection        # m
        self.Side = self.AdjustedAntennaCircumferance/self.NumberSides                              # m
        self.X = self.Side/np.sqrt(2)                                                               # m
        self.Diameter = self.Side + 2*self.X                                                        # m
        self.CircumscribedRadius = self.Side/(2*np.sin(45*np.pi/180/2))                             # m
        
        # This is the corrected loop area, factorin in the additional length
        self.Area = 2*(1+np.sqrt(2))*self.Side**2 + 0.5*self.WidthCorrection*(self.Side+2*self.X)   # m^2
        
        self.OperatingFrequencies = self.Data.InputReactanceNEC[0,:]                                # MHz
        self.AngularFrequencies = 2*np.pi*self.OperatingFrequencies*10**6                           # Rad/s
        self.Wavelengths = sp.constants.speed_of_light/(self.OperatingFrequencies*10**6)            # m
        self.LengthDiameterRatios = np.arange(start=.1,
                                              stop=60,
                                              step=0.001)
        
    def SelfResonantFrequency(self):
            
        return self.Data.InputReactanceNEC_AtResonance[0,np.argmin(np.abs(self.Data.InputReactanceNEC_AtResonance[1,:]))]
  
    def Convert_S1P(self):
            
        VNA = self.Z_in
            
        z0 = 50.0
        ReS11 = self.Data.NanoVNA[1,:]
        ImS11 = self.Data.NanoVNA[2,:]
            
        # Where do these expressions come from?
        
        ReZ = z0*(1-ReS11**2-ImS11**2)/((1-ReS11)**2+ImS11**2)
        
        ImZ = z0*2*ImS11/((1-ReS11)**2+ImS11**2)
        
        VNA.Reactance = ImZ
        VNA.Resistance = ReZ
        VNA.Frequency = ConvertUnits(Attribute=self.Data.NanoVNA[0,:],
                                     ConvertFrom='Hz',
                                     ConvertTo='MHz')
        
        return VNA
        
        
    class Z_in:
        def __init__(self):
            self.Reactance
            self.Resistance
            self.Phase
            self.Frequency
    
    class ExternalInductance_Equations:
        def __init__(self):        
            self.ARRL
            self.Grover
            self.Knight
            self.RayleighNiven
    
    class SelfCapacitance_Equations:
        def __init__(self):
            self.Medhurst
            self.CL_DAE
            self.Macro
            self.NEC2
            self.SheathHelix
            self.ButlerUpper
            self.ButlerLower
            self.Smythe
            self.McDonald
    
    def WavelengthFraction(self):
        
        return self.AntennaCircumferance/self.Wavelengths    
               
    def RadiationResistance(self):
        
        return 20.0*np.pi**2*(self.AntennaCircumferance/self.Wavelengths)**4*self.NumberTurns**2 
        
    def HighFrequencyResistance(self):
        
        return self.AntennaCircumferance/(2*np.pi*self.ConductorRadius)*self.SurfaceResistance()
        
    def SurfaceResistance(self):
        
        return (self.AngularFrequencies*self.LocalConstants.PermeabilityVacuum/(2*self.LocalConstants.CopperConductivity))**0.5
           
    def ExternalInductance(self):
             
        Method = self.ExternalInductance_Equations
        
        B2R = self.ConductorRadius/self.DerivedRadius
        SNL = self.Side*self.NumberTurns/((self.NumberTurns+1)*self.ConductorDiameter)
        RRw = self.DerivedRadius/self.ConductorRadius
        L2D = self.ConductorDiameter/(self.AntennaCircumferance/np.pi)
        
        Method.ARRL = 0.016*self.NumberTurns**2*self.Side*100*(np.log(2.613*SNL)+0.75143+(0.07153*1/SNL))                                                        # Micro Henries
        
        Method.Grover = 0.016*self.NumberTurns**2*self.Side*100*(np.log(1/B2R)+0.75143+0.18693*B2R+0.11969*B2R**2-0.08234*B2R**4)                                # Micro Henries
        
        Method.Knight = self.LocalConstants.PermeabilityAir*self.DerivedRadius*self.NumberTurns**2*(np.log(1+np.pi/(2*L2D)) + 1/(2.3004+3.2219*L2D+1.7793*L2D*2))*10**6    # Micro Henries
        
        Method.RayleighNiven = self.LocalConstants.PermeabilityAir*self.DerivedRadius*self.NumberTurns**2*(np.log(8*RRw)-2+1/8*RRw**-2*(np.log(8*RRw)+1/3))*10**6          # Micro Henries
        
        return Method
    
    def InternalInductance(self):
        
        return 10**6*self.AntennaCircumferance/(self.AngularFrequencies*2*np.pi*self.ConductorRadius)*(self.AngularFrequencies*self.LocalConstants.PermeabilityCopper/(2*self.LocalConstants.CopperConductivity))**0.5   # Micro Henries            
            
    def InputImpedance(self):
        
        ImpedanceIn = self.Z_in
        
        ImpedanceIn.Reactance = (self.AngularFrequencies*(self.InternalInductance() + self.ExternalInductance().Knight)*10**-6 - 1/(self.SelfCapacitance().McDonald*self.AngularFrequencies)*10**-6)    
        
        ImpedanceIn.Resistance = self.RadiationResistance() + self.HighFrequencyResistance()
        
        ImpedanceIn.Phase = np.arctan(ImpedanceIn.Reactance/ImpedanceIn.Resistance*np.pi/180)*180/np.pi
        
        return ImpedanceIn
    
    def ResonatingCapacitance(self):
        
        # Antenna Theory Analysis and Design 3rd Edition, pg. 257, eq. 5-35
        return self.AngularFrequencies**-1*self.InputImpedance().Reactance/(self.InputImpedance().Resistance**2+self.InputImpedance().Reactance**2)*10**12       # Pico Farads

    def SelfCapacitance(self):
        
        Method = self.SelfCapacitance_Equations
        
        L2D = self.LengthDiameterRatios

        Kc = 0.717439/L2D + 0.933048/L2D**1.5 + 0.106/L2D**2
        Km = 0.8249/L2D + 2.329/L2D**1.5
        
        eps0 = self.LocalConstants.PermeativityVacuum

        eps_poly = 2.56          # Polystyrene, as reported by Medhurst
        eps_air = 1.00058986     # Air at STP, 900 KHz
        eps_vac = 1.0            # Vacuum
        
        eps_rel = eps_poly
        
        psi = np.arctan(L2D/(self.NumberTurns*np.pi)/180)
        
        Method.Medhurst = 4/np.pi*eps0*L2D*(1+Km)*10**12                                                             # pF/m
        
        Method.CL_DAE = 4/np.pi*eps0*eps_air*L2D*(1+Kc*(1+eps_poly/eps_air)/2)*10**12*1/(np.cos(psi))**2             # pF/m
        
        Method.Macro = 4/np.pi*eps0*eps_air*L2D*(1+Kc*(1+eps_poly/eps_air)/2)*10**12                                 # pF/m
        
        Method.NEC2 = 1/self.ExternalInductance().Knight*10**6*(2*np.pi*self.SelfResonantFrequency()*10**6)**-2*10**12 # pF
        
        kD = self.NumberTurns*np.pi*0.5*1/(L2D)*eps_rel
        
        lnkD = np.log(kD)
        
        den = lnkD - 0.8091
        
        num1 = lnkD - 1.309
        
        num2 = lnkD**2 - 2.668*lnkD + 1.904
        
        Method.SheathHelix = 2*np.pi*eps0*(1/den + 0.1250*num1/den**2*kD**2 + 0.009766*num2/den**3*kD**4)*10**12 
        
        Method.ButlerLower = 2*np.pi**2*eps0/np.log(16/L2D*eps_rel)*10**12
        
        Method.ButlerUpper = -2*np.pi*1/L2D*eps0/(1+np.log(L2D*1/(eps_rel*4)))*10**12
        
        Method.Smythe = eps0*(4+3.475*(L2D/eps_rel)**0.76)*10**12
        
        Method.McDonald = eps0*self.DerivedRadius/np.log(self.DerivedRadius/self.ConductorRadius)
        
        return Method
    
    def Plot(self):
         
        # Setting Plot Formats
        plt.rcParams["font.family"] = "Palatino Linotype"
        
        # Creating Plot of Calculated Resonating Capacitance vs Operating Frequency
        fig, ax1 = plt.subplots(figsize=(12,6))
        
        ax1.set_yscale('log')
        #ax1.set_xscale('log')
        
        ax1.set_xlabel(r'Operating Frequency [MHz]',
                       fontname="Palatino Linotype",
                       fontsize=12)
    
        ax1.set_ylabel(r'Resonating Capacitance [pF]',
                       fontname="Palatino Linotype",
                       fontsize=12)
    
        ax1.tick_params(right=True, 
                        labelright=True)
        
        ax1.plot(self.OperatingFrequencies,
                 self.ResonatingCapacitance(),
                 "-",
                 markersize=0.1,
                 label="Resonating Capacitance (Calculated)",
                 alpha=0.7)
    
        ax1.plot(self.Data.ResonatingCapacitanceNEC[0,:],
                 self.Data.ResonatingCapacitanceNEC[1,:],
                 ".",
                 markersize=1,
                 label="NEC2 Paralell Equivalent Capacitance",
                 alpha=0.9)

        ax1.axvline(x=self.IntendedOperatingFrequency,
                    ymin=0.0, 
                    ymax=1.0, 
                    color='black',
                    label='Operating Frequency = %s MHz' %self.IntendedOperatingFrequency,
                    linewidth=0.1,
                    linestyle='--')
    
        ax1.legend(loc=0)
    
        # Saving Plot Output
        plt.tight_layout()
        plt.savefig(path_to_plots+"Resonating_Capacitance_Required.pdf")

        # Creating Plot of Calculated Radiation Resistance vs Operating Frequency
        fig, ax1 = plt.subplots(figsize=(12,6))
            
        ax1.set_yscale('log')
        
        ax1.set_xlabel(r'Operating Frequency [MHz]',
                       fontname="Palatino Linotype",
                       fontsize=12)
    
        ax1.set_ylabel(r'Resistance [Ohms]',
                       fontname="Palatino Linotype",
                       fontsize=12)
    
        ax1.tick_params(right=True, 
                        labelright=True)
        
        ax1.plot(self.OperatingFrequencies,
                 self.RadiationResistance(),
                 "-",
                 markersize=1,
                 label="Radiation Resistance (Calculated)",
                 alpha=0.9)
    
        ax1.plot(self.OperatingFrequencies,
                 self.HighFrequencyResistance(),
                 "-",
                 markersize=1,
                 label="HF Resistance (Calculated)",
                 alpha=0.9)    
        
        ax1.axvline(x=self.IntendedOperatingFrequency,
                    ymin=0.0, 
                    ymax=1.0, 
                    color='black',
                    label='Operating Frequency = %s MHz' %self.IntendedOperatingFrequency,
                    linewidth=0.1,
                    linestyle='--')
        
        ax1.legend(loc=0)
    
        # Saving Plot Output
        plt.tight_layout()
        plt.savefig(path_to_plots+"Radiation_Resistance.pdf")

        # Creating Plot of Calculated Resistance and Reactance
        fig, ax1 = plt.subplots(figsize=(18,10))
        
        ax2 = ax1.twinx()
        
        ax1.set_xlim(xmax=80.0)
        ax2.set_xlim(xmax=80.0)
        
        ax1.set_yscale('log')
        ax2.set_yscale('symlog')
        
        ax1.set_xlabel(r'Operating Frequency [MHz]',
                       fontname="Palatino Linotype",
                       fontsize=12)
    
        ax1.set_ylabel(r'Resistance [Ohms]',
                       fontname="Palatino Linotype",
                       fontsize=12)

        ax2.set_ylabel(r'Reactance [Ohms]',
                       fontname="Palatino Linotype",
                       fontsize=12)
        
        ax1.tick_params(right=False, 
                        labelright=False)
        
        param2 = 0.5*self.SelfResonantFrequency()**2/((1.35*np.pi)**-6*self.OperatingFrequencies**2 + (self.OperatingFrequencies - self.SelfResonantFrequency())**2)

        test = self.InputImpedance().Resistance*param2
        nec =self. Data.InputResistanceNEC[1,:]
        
        diff = np.trapz(test) - np.trapz(nec)
        
        #print(test.size)
        #print(nec.size)
        #print(diff)
        
        # Calculated Resistance
        ax1.plot(self.OperatingFrequencies,
                 self.InputImpedance().Resistance*param2,
                 "-",
                 markersize=1,
                 label="HF + Rad Resistance (Calculated)",
                 color='royalblue',
                 alpha=0.9)
        
        param5 = np.tan((self.SelfResonantFrequency()-self.OperatingFrequencies)/self.SelfResonantFrequency())*param2
        
        # Calculated Reactance
        ax2.plot(self.OperatingFrequencies,
                 self.InputImpedance().Reactance*param5,
                 "-",
                 markersize=1,
                 label="Int + Ext Inductor Reactance (Calculated)",
                 color='orange',
                 alpha=0.9)  
    
        # NEC Simulated Resistance
        ax1.plot(self.Data.InputResistanceNEC[0,:],
                 self.Data.InputResistanceNEC[1,:],
                 "--",
                 markersize=3,
                 label="Total Resistance (NEC2 Simulated)",
                 color='royalblue',
                 alpha=0.9)    
        
        # NEC Simulated Reactance
        ax2.plot(self.Data.InputReactanceNEC[0,:],
                 self.Data.InputReactanceNEC[1,:],
                 "o",
                 markersize=3,
                 label="Total Reactance (NEC2 Simulated)",
                 color='orange',
                 alpha=0.9)

        # EZNEC Pro/2++ 7.0 Simulated Resistance
        ax1.plot(self.Data.InputResistanceReactance_EZNEC[0,:],
                 self.Data.InputResistanceReactance_EZNEC[2,:],
                 "--",
                 markersize=3,
                 label="Total Resistance (EZNEC Pro Simulated)",
                 color='green',
                 alpha=0.9)    
        
        # EZNEC Pro/2++ 7.0 Simulated Reactance
        ax2.plot(self.Data.InputResistanceReactance_EZNEC[0,:],
                 self.Data.InputResistanceReactance_EZNEC[3,:],
                 "o",
                 markersize=3,
                 label="Total Reactance (EZNEC Pro Simulated)",
                 color='brown',
                 alpha=0.9)

        # NanoVNA Measured Resistance w/Conversion
        ax1.plot(self.Convert_S1P().Frequency,
                 self.Convert_S1P().Resistance,
                 ".",
                 markersize=1,
                 label="Total Resistance (NanoVNA)",
                 color='grey',
                 alpha=0.9)
        
        # NanoVNA Measured Reactance w/Conversion
        ax2.plot(self.Convert_S1P().Frequency,
                 self.Convert_S1P().Reactance,
                 ".",
                 markersize=1,
                 label="Total Reactance (NanoVNA)",
                 color='black',
                 alpha=0.9)
        
        ax2.axvline(x=self.IntendedOperatingFrequency,
                    ymin=0.0, 
                    ymax=1.0, 
                    color='black',
                    label='Operating Frequency = %s MHz' %self.IntendedOperatingFrequency,
                    linewidth=0.1,
                    linestyle='--')
        
        ax2.axvline(x=self.SelfResonantFrequency(), 
                    ymin=0.0, 
                    ymax=1.0, 
                    color='black',
                    label='Self Resonant Frequency (NEC2 Simulated) = %s MHz' %self.SelfResonantFrequency(),
                    linewidth=0.1,
                    linestyle='--')
        
        # Combing the Axis Lines and Labels into a single Legend
        lines_ax1, labels_ax1 = ax1.get_legend_handles_labels()
        lines_ax2, labels_ax2 = ax2.get_legend_handles_labels()
        lines = lines_ax1 + lines_ax2
        labels = labels_ax1 + labels_ax2
        
        ax1.legend(lines, 
                   labels, 
                   loc=0)
        
        # Saving Plot Output
        plt.tight_layout()
        plt.savefig(path_to_plots+"Resistance_Reactance.pdf")

        # Creating Plot of Calculated Phase Angle
        fig, ax1 = plt.subplots(figsize=(12,6))
        
        #ax1.set_yscale('symlog')
        #ax1.set_xscale('log')
        
        ax1.set_xlabel(r'Operating Frequency [MHz]',
                       fontname="Palatino Linotype",
                       fontsize=12)
    
        ax1.set_ylabel(r'Phase Angle []',
                       fontname="Palatino Linotype",
                       fontsize=12)
        
        ax1.plot(self.OperatingFrequencies,
                 self.InputImpedance().Phase,
                 "-",
                 markersize=0.1,
                 label="Phase",
                 alpha=0.7)
        
        ax1.plot(self.Data.PhaseNEC[0,:],
                 self.Data.PhaseNEC[1,:],
                 "--",
                 markersize=0.1,
                 label="Phase (NEC2 Simulated)",
                 alpha=0.7)
        
        ax1.plot(self.OperatingFrequencies,
                 np.arctan((self.SelfResonantFrequency() - self.OperatingFrequencies)/self.SelfResonantFrequency()*95)*1.0025*180/np.pi,
                 "--",
                 markersize=0.1,
                 label="Phase (TAN)",
                 alpha=0.7)
        
        ax1.legend(loc=0)
    
        # Saving Plot Output
        plt.tight_layout()
        plt.savefig(path_to_plots+"Phase.pdf")

        # Creating Plot of Medhurst Experimental Data and Equations
        fig, ax1 = plt.subplots(figsize=(12,6))
                    
        ax1.set_yscale('log')
        ax1.set_xscale('log')

        ax1.set_xlabel(r'Length/Diameter [-]',
                       fontname="Palatino Linotype",
                       fontsize=12)
    
        ax1.set_ylabel(r'Self Capacitance [pF/m]',
                       fontname="Palatino Linotype",
                       fontsize=12)
    
        ax1.tick_params(right=True, 
                        labelright=True)
        
        # Converting Data from pF/cm to pF/m
        ax1.plot(self.Data.MedhurstData[0,:],
                 ConvertUnits(self.Data.MedhurstData[1,:],
                              'Centimeters',
                              'Meters'),
                 "o",
                 markersize=3,
                 label="Medhurst Data",
                 alpha=0.9)
    
        ax1.plot(self.LengthDiameterRatios,
                 self.SelfCapacitance().CL_DAE,
                 "-",
                 markersize=1,
                 label="CL-DAE",
                 alpha=0.9)
        
        ax1.plot(self.LengthDiameterRatios,
                 self.SelfCapacitance().Medhurst,
                 "-",
                 markersize=1,
                 label="Medhurst Refitted",
                 alpha=0.9)
                
        ax1.plot(self.LengthDiameterRatios,
                 self.SelfCapacitance().Smythe,
                 "-",
                 markersize=1,
                 label="Smythe",
                 alpha=0.9)
        
        ax1.legend(loc=0)
    
        # Saving Plot Output
        plt.tight_layout()
        plt.savefig(path_to_plots+"Medhurst_Data.pdf")

        # Creating Plot of Resistance and Reactance at Resonance
        fig, ax1 = plt.subplots(figsize=(12,6))
        
        ax2 = ax1.twinx()
        
        ax1.set_yscale('symlog')
          
        ax1.set_xlabel(r'Operating Frequency [MHz]',
                       fontname="Palatino Linotype",
                       fontsize=12)
    
        ax1.set_ylabel(r'Reactance [Ohms]',
                       fontname="Palatino Linotype",
                       fontsize=12)

        ax1.axvline(x=self.SelfResonantFrequency(), 
                    ymin=0.0, 
                    ymax=1.0, 
                    color='black',
                    label='Self Resonant Frequency = %s MHz' %self.SelfResonantFrequency(),
                    linewidth=0.1,
                    linestyle='--')
        
        ax1.axhline(y=0.0, 
                    xmin=0.0, 
                    xmax=1.0, 
                    color='black',
                    linewidth=0.1,
                    linestyle='--')
    
        ax1.plot(self.Data.InputReactanceNEC_AtResonance[0,:],
                 self.Data.InputReactanceNEC_AtResonance[1,:],
                 ".",
                 markersize=1,
                 label="Reactance (NEC2 Simulated)",
                 color='royalblue',
                 alpha=0.9)       

        ax2.set_yscale('log')
        
        ax2.set_ylabel(r'Resistance [Ohms]',
                       fontname="Palatino Linotype",
                       fontsize=12)
        
        ax2.plot(self.Data.InputResistanceNEC_AtResonance[0,:],
                 self.Data.InputResistanceNEC_AtResonance[1,:],
                 ".",
                 markersize=1,
                 label="Resistance (NEC2 Simulated)",
                 color='red',
                 alpha=0.9)
      
        # Combing the Axis Lines and Labels into a single Legend
        lines_ax1, labels_ax1 = ax1.get_legend_handles_labels()
        lines_ax2, labels_ax2 = ax2.get_legend_handles_labels()
        lines = lines_ax1 + lines_ax2
        labels = labels_ax1 + labels_ax2
        
        ax1.legend(lines, 
                   labels, 
                   loc=0)    
    
        # Saving Plot Output
        plt.tight_layout()
        plt.savefig(path_to_plots+"Resistance_Reactance_At_Resonance.pdf")

        return

    def GeometryCheck(self):
        
        print('Diameter (Octagon)          ', self.Diameter)
        print('Circumscribed Radius        ', self.CircumscribedRadius)
        print('Circumference Derived Radius', self.AntennaCircumferance/(2*np.pi))

        return      
        
    # Not needed for single turn loop, HF Resistance is used, if more than one turn, then use this equation.
    def OhmicResistance(self):
        
        return self.NumberTurns*self.AntennaCircumferance/(2*np.pi)*1/(self.ConductorRadius)*self.SurfaceResistance()*(0+1)
    
    def SkinEffectResistance(self):
        
        return self.NumberTurns/(2*np.pi*self.ConductorRadius)*self.SurfaceResistance()