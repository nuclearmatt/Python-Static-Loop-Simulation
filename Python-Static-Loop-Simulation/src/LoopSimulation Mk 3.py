# Importing Required Modules
from headers.StandardModules import *
from headers.DefaultFilePaths import *
from classes.LoopAntenna_Class import LoopAntenna
from classes.LoadData_Class import LoadData
from classes.LocalConstants_Class import LocalConstants

OctagonLoop = LoopAntenna(ConductorDiameter=0.875,
                          AntennaCircumferance=25.0,
                          WidthCorrection=0.0,
                          NumberTurns=1,
                          UpperOperatingFrequency=40.0,
                          LowerOperatingFrequency=1.5,
                          IntendedOperatingFrequency=7.1,
                          Model='EqualSides',
                          LocalConstants=LocalConstants(),
                          Data=LoadData())    

OctagonLoop.Plot()

OctagonLoop.GeometryCheck()  
    
print('Knight         L (uH): ', OctagonLoop.ExternalInductance().Knight)
print('Knight         C (pF): ', OctagonLoop.SelfCapacitance().NEC2)
print('Grover         L (uH): ', OctagonLoop.ExternalInductance().Grover)
print('Rayleigh-Niven L (uH): ', OctagonLoop.ExternalInductance().RayleighNiven)
print('ARRL           L (uH): ', OctagonLoop.ExternalInductance().ARRL)
    