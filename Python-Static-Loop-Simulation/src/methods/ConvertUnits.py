def ConvertUnits(Attribute,
                 ConvertFrom,
                 ConvertTo,
                 **kwargs):
    
    Converted=False
    
    if (ConvertFrom == 'Inches'):
        if (ConvertTo == 'Meters'):
            ConversionFac = 1/12*1200/3937.0
            Converted=True
            
    if (ConvertFrom == 'Feet'):
        if (ConvertTo == 'Meters'):
            ConversionFac = 1200/3937.0
            Converted=True
    
    if (ConvertFrom == 'Centimeters'):
        if (ConvertTo == 'Meters'):
            ConversionFac = 100.0
            Converted=True
            
    if (ConvertFrom == 'Hz'):
        if (ConvertTo == 'MHz'):
            ConversionFac = 10**-6            
            Converted=True
        
    if not Converted:
        
        print('Error: Unknown Conversion Combination: Convert To = %s, Convert From = %s' %(ConvertTo, ConvertFrom))
        ConversionFac = 0.0

    return Attribute*ConversionFac    