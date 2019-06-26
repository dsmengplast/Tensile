# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 15:03:47 2019

@author: 364970
"""

def consolidate(folder):
    
    #user input to name final excel file
    print('What to name the output xlsx?')
    name = input()
    print('What to name the overlay?')
    title = input()

    import os, glob
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    import csv
    from scipy import stats
    #sets working directory to the folder with raw data files
    os.chdir(folder)
    
    # .glob returns a list of filenames matching a pattern 
    # * represents any text. this line returns a list of all .csv files in the directory
    filenames = glob.glob('*.csv')  
    
    #creates a list of dataframes with ALL Instron generated data - one per .csv file 
    dfs = [pd.read_csv(filename, header = 6) for filename in filenames]
    
    #pulls stress and strain values from the dataframes as series and lists them
    strain = [pd.read_csv(filename, header=6).iloc[1: , 3].astype(float) for filename in filenames]
    stress = [pd.read_csv(filename, header=6).iloc[1: , 5].astype(float) for filename in filenames]
    time =  [pd.read_csv(filename, header=6).iloc[1: , 0].astype(float) for filename in filenames]
   
    for i in range(len(list(time))):
        Tcrit = dfs[i]['Time'][1:].astype(float).max() -.100
        tyme = time[i].where(time[i] > Tcrit).dropna()
        breakstress = stress[i][-len(tyme): ].max()
        breakindex = stress[i].where(stress[i] == breakstress).dropna().index[0]
        breakstrain = dfs[i]['Tensile strain (AVE2)'][breakindex]
        
        
    
    
    #finds maximum stress and strain to set x and y axis for overlay 
    stress_df = pd.concat(stress[:])
    max_stress = stress_df.max()
    stress_axis_limit = max_stress*1.1
    strain_df = pd.concat(strain[:])
    max_strain = strain_df.max()
    strain_axis_limit = max_strain*1.1
    longest_strain = len(strain_df)
    
    #creates a list of dataframes for specimen info from the first 6 lines        
    headers = [pd.read_csv(filename, skiprows = range(6, longest_strain), names = ['Property', 'Value', 'Unit']) for filename in filenames]
        
    #creates the overlay
   
    plt.figure()
    for i in range(len(list(filenames))):
        grd = pd.Series(np.gradient(stress[i]))
        grd = grd.where(grd>-1).dropna()
        plt.plot(strain[i][0:len(grd)], stress[i][0:len(grd)], label= 'Specimen ' + str(i+1))
        plt.legend()
        plt.xlabel('Tensile Strain (%)')
        plt.ylabel('Tensile Stress (MPa)')
        plt.ylim(ymin=0, ymax = stress_axis_limit)
        plt.xlim(xmin=0, xmax = strain_axis_limit)
        plt.title(title)
        plt.savefig('Overlay')
    
               
    #creates individal plots and defines stress+strain @ break
    SAB_list = []
    StressAB_list = []
    Stress_yield_list = []
    Strain_yield_list = []
    for i in range(len(list(filenames))):
                strains = dfs[i].iloc[1: , 3].astype(float)
                stresses = dfs[i].iloc[1: , 5].astype(float)
                
                Tcrit = dfs[i]['Time'][1:].astype(float).max() -.100
                tyme = time[i].where(time[i] >= Tcrit).dropna()
                breakstress = stress[i][-len(tyme): ].max().astype(float)
                breakindex = stress[i].where(stress[i] == breakstress).dropna().index[0]
                breakstrain = float(dfs[i]['Tensile strain (AVE2)'][breakindex])
                
                plt.figure()
                plt.plot(strains[0:breakindex], stresses[0:breakindex])
                SAB_list.append(breakstrain)
                StressAB_list.append(breakstress)
                plt.xlabel('Tensile Strain (%)')
                plt.ylabel('Tensile Stress (MPa)')
                plt.ylim(ymin=0)
                plt.xlim(xmin=0)
                plt.title('Specimen' + str(i+1))
                plt.savefig('Plot' + str(i+1))
                plt.close()
                
                num_of_regions = int(100)
                overlap = 0.5
                len_region = (len(strains)/num_of_regions)*(1/(1-overlap))
                for i in range(num_of_regions):
                    strain_region = strains[int((i* len_region)-(i*len_region*overlap)) : int(((i+1)*len_region)-(i*len_region*overlap))]
                    stress_region = stresses[int((i* len_region)-(i*len_region*overlap)) : int(((i+1)*len_region)-(i*len_region*overlap))]
                    r= stats.linregress(strain_region, stress_region)
                    if r.slope <=0:
                        ystress = stress_region.max()
                        print(ystress)
                        Stress_yield_list.append(ystress)
                        yindex = stress_region.where(stress_region == ystress).dropna().index
                        ystrain = float(strain_region[yindex])
                        Strain_yield_list.append(ystrain)
                        break
                        
                

                    
    

    
    #Calculates Modulus
    mod_list =[]
    for i in range(len(list(filenames))):
        strains = dfs[i].iloc[1: , 3].astype(float)
        stresses = dfs[i].iloc[1: , 5].astype(float)
        ss = [strains, stresses]
        ss = pd.DataFrame(ss).T    
        
        #interpolation to find stress at exactly .05, .25% strain
        ss_low = ss.where(ss['Tensile strain (AVE2)'] <= .05).dropna()
        
        ss_mod = ss.where(ss['Tensile strain (AVE2)'] > .05).dropna()
        ss_mod = ss.where(ss['Tensile strain (AVE2)'] < .25).dropna()
        
        ss_high = ss.where(ss['Tensile strain (AVE2)'] >= .25).dropna()
        
        x0A = ss_low['Tensile strain (AVE2)'].max()
        x1A = ss_mod['Tensile strain (AVE2)'].min()
        y0A = ss_low['Tensile stress'].max()
        y1A = ss_mod['Tensile stress'].min()
        
        x0B = ss_mod['Tensile strain (AVE2)'].max()
        x1B = ss_high['Tensile strain (AVE2)'].min()
        y0B = ss_mod['Tensile stress'].max()
        y1B = ss_high['Tensile stress'].min()
        
       
        x1 = .05
        x2 = .25
        y1 = (y0A*(x1A-x1)+y1A*(x1-x0A))/(x1A-x0A)
        y2 = (y0B*(x1B-x2)+y1B*(x2-x0B))/(x1B-x0B)
        rise = (y2-y1)
        run = (x2-x1)/100
        modulus = rise/run
        mod_list.append(modulus)
    
        
    #Combines Stress at yield, Strain at yield, and Modulus     
    labels = []
    for i in range(len(list(filenames))):
        labels.append('Specimen '+str(i+1))
    SSM = [mod_list, Strain_yield_list, Stress_yield_list]
    SSM = pd.DataFrame(SSM, index = ['Modulus (MPa)', 'Strain at Yield (%)', 'Stress at Yield (MPa)'], columns=labels)
    dev = SSM.std(axis=1).round(decimals=3)
    avg= SSM.mean(axis=1).round(decimals=1)
    
    SSM['AVG'] = avg
    SSM['DEV'] = dev
#    print(SSM)
    
    
    
#####################################################################################################################################################
    #Sets an excel file to write all dataframes and inserts plots 
    with pd.ExcelWriter(name) as writer: 
        for i in range(len(list(filenames))):
            headers[i].to_excel(writer, sheet_name='Specimen Properties', index = False, startrow=i*8) #sends headers to the generated excel file
            dfs[i].to_excel(writer, sheet_name = 'Sample' + str(i+1) + 'RawData', index = False) #sends dfs to their own sheet in the generated excel file
            worksheet = writer.sheets['Sample' + str(i+1) + 'RawData']
            worksheet.insert_image('K1', 'Plot' + str(i+1) + '.png')
        SSM.to_excel(writer, sheet_name='Report')
        reportsheet = writer.sheets['Report']
        reportsheet.insert_image('A6', 'Overlay.png')
        reportsheet.set_column(0,15,20)
   
    
    #Writes a CSV file for LIMS upload - might need to be as .erg?
#    print('what is the test #?')
#    tn = input()
#    with open (tn+'.csv', mode = 'w') as csv_file:
#        fieldnames = ['Testnumber', 'E modulus_x', 'Stress at 0.5% Strain_x', 'Stress at 1.0% Strain_x', 'Stress at 2.0% Strain_x', 'Stress at 5% Strain_x', 'Stress at 10% Strain_x', 'Stress at 50% Strain_x', 'Stress at 100% Strain_x', 'Stress at 300% Strain_x', 'Yield stress_x', 'Stress at Break_x', 'Maximum Stress_x', 'Yield Strain_x', 'Strain at Break_x', 'Cross sectional area_x', 'Thickness_x', 'Width_x']
#        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
#        writer.writeheader()
#        writer.writerow({'Testnumber': tn, 'E modulus_x': SSM['AVG']['Modulus (MPa)'], 'Stress at 0.5% Strain_x': ' ', 'Stress at 1.0% Strain_x': ' ', 'Stress at 2.0% Strain_x': ' ', 'Stress at 5% Strain_x': ' ', 'Stress at 10% Strain_x': ' ', 'Stress at 50% Strain_x': ' ', 'Stress at 100% Strain_x': ' ', 'Stress at 300% Strain_x': ' ', 'Yield stress_x': ' ', 'Stress at Break_x': SSM['AVG']['Stress at Break (MPa)'], 'Maximum Stress_x': ' ', 'Yield Strain_x': ' ', 'Strain at Break_x': SSM['AVG']['Strain at Break (%)'], 'Cross sectional area_x': ' ', 'Thickness_x' : ' ', 'Width_x': ' ' })
#        
#                     
consolidate('C:\\Users\\364970\\Documents\\python\\19-0123\\F232-D Natural\\5-17-2019 19-0123 F232-D Natural Lot # 0000178323.is_tens_RawData')