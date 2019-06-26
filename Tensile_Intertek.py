# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 15:03:30 2019

@author: 364970
"""

def consolidate(folder):
    
    #user input to name final excel file
    print('What to name the output xlsx?')
    name = input()

    import os, glob
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    import csv

    #sets working directory to the folder with raw data files
    os.chdir(folder)
    
    # .glob returns a list of filenames matching a pattern 
    # * represents any text. this line returns a list of all .csv files in the directory
    filenames = glob.glob('*.csv')  
    
    #creates a list of dataframes with ALL Instron generated data - one per .csv file 
    dfs = [pd.read_csv(filename, header = 2) for filename in filenames]
    #longest_df = []
    #for i in range(len(list(filenames))):
        
    
    #pulls stress and strain values from the dataframes as series and lists them
    strain = [pd.read_csv(filename, header=2).iloc[: , 1].astype(float) for filename in filenames]
    stress = [pd.read_csv(filename, header=2).iloc[: , 0].astype(float) for filename in filenames]

    
    #finds maximum strain to use for x axis limit of overlay 
    stress_df = pd.concat(stress[:])
    max_stress = stress_df.max()
    stress_axis_limit = max_stress*1.1
    strain_df = pd.concat(strain[:])
    max_strain = strain_df.max()
    strain_axis_limit = max_strain*1.1
    
    
        
    #creates the overlay
    labels = []
    plt.figure(figsize = (20,10))
    for i in range(len(list(filenames))):
        grd = pd.Series(np.gradient(stress[i]))
        #plt.figure(figsize=(20,10))
       # plt.plot(grd, 'o')
        grd = grd.where(grd>-.1).dropna()
        #plt.plot(grd)
#        plt.plot(strain[i], stress[i], 'bo')
#        plt.plot(strain[i][len(grd)], stress[i][len(grd)],'ro')
            
        #plt.figure(figsize = (20,10))
        plt.plot(strain[i][0:len(grd)], stress[i][0:len(grd)], label= 'Specimen '+ str(i+1))
        plt.legend(fontsize = 18)
        plt.xlabel('Tensile Strain (%)', fontsize = 18)
        plt.ylabel('Tensile Stress (MPa)', fontsize = 18)
        plt.ylim(ymin=0, ymax = stress_axis_limit)
        plt.xlim(xmin=0, xmax = strain_axis_limit)
        #plt.title('Specimen 1-' + str(len(filenames))+ ' Overlay')
        plt.title('ForTii MX53T', fontsize = 20)
        plt.savefig('Overlay')
    
               
    #creates individal plots and defines stress+strain @ break
    SAB_list = []
    StressAB_list = []
    
    for i in range(len(list(filenames))):
                strains = dfs[i].iloc[: , 1].astype(float)
                stresses = dfs[i].iloc[: , 0].astype(float)
                grd = pd.Series(np.gradient(stresses))
                grd1 = grd.where(grd>-.1).dropna()
                grd2 = grd.where(grd<-.1).dropna()           
                plt.figure(figsize = (20, 10))
                plt.plot(grd1, 'ro')
                plt.plot(grd2, 'bo')
                
#                grd3 = savgol_filter(grd, window_length = 49, polyorder = 11)
#                plt.figure(figsize = (20,10))
#                plt.plot(grd3, 'go')

                plt.figure()
                plt.plot(strains[0:len(grd1)], stresses[0:len(grd1)])
                SAB = strains[0:len(grd1)+1].max()
                SAB_list.append(SAB)
                StressAB = stresses[0:len(grd)+1].max()
                StressAB_list.append(StressAB)
                plt.xlabel('Tensile Strain (%)')
                plt.ylabel('Tensile Stress (MPa)')
                plt.ylim(ymin=0)
                plt.xlim(xmin=0)
                plt.title('Specimen' + str(i+1))
                plt.savefig('Plot' + str(i+1))
                
    
   
    
    
    #Calculates Modulus
    mod_list =[]
    for i in range(len(list(filenames))):
        strains = dfs[i].iloc[: , 1].astype(float)
        stresses = dfs[i].iloc[: , 0].astype(float)
        ss = [strains, stresses]
        ss = pd.DataFrame(ss).T    
        
        #interpolation to find stress at exactly .05, .25% strain
        ss_low = ss.where(ss['(%)'] <= .05).dropna()
        
        ss_mod = ss.where(ss['(%)'] > .05).dropna()
        ss_mod = ss.where(ss['(%)'] < .25).dropna()
        
        ss_high = ss.where(ss['(%)'] >= .25).dropna()
        
        x0A = ss_low['(%)'].max()
        x1A = ss_mod['(%)'].min()
        y0A = ss_low['(MPa)'].max()
        y1A = ss_mod['(MPa)'].min()
        
        x0B = ss_mod['(%)'].max()
        x1B = ss_high['(%)'].min()
        y0B = ss_mod['(MPa)'].max()
        y1B = ss_high['(MPa)'].min()
        
       # plt.figure()
        #plt.plot(ss_mod['Tensile strain (AVE2)'], ss_mod['Tensile stress'])
        x1 = .05
        x2 = .25
        y1 = (y0A*(x1A-x1)+y1A*(x1-x0A))/(x1A-x0A)
        y2 = (y0B*(x1B-x2)+y1B*(x2-x0B))/(x1B-x0B)
        rise = (y2-y1)
        run = (x2-x1)/100
        modulus = rise/run
        mod_list.append(modulus)
    
        
    #Combines Stress at break, Strain at Break, and Modulus     
    labels = []
    for i in range(len(list(filenames))):
        labels.append('Specimen '+str(i+1))
    SSM = [mod_list, SAB_list, StressAB_list]
    SSM = pd.DataFrame(SSM, index = ['Modulus (MPa)', 'Strain at Break (%)', 'Stress at Break (MPa)'], columns=labels)
    dev = SSM.std(axis=1).round(decimals=3)
    avg= SSM.mean(axis=1).round(decimals=1)
    
    SSM['AVG'] = avg
    SSM['DEV'] = dev
    
    
    
    
#####################################################################################################################################################
    #Sets an excel file to write all dataframes and inserts plots 
    with pd.ExcelWriter(name) as writer: 
        for i in range(len(list(filenames))):
            dfs[i].to_excel(writer, sheet_name = 'Sample' + str(i+1) + 'RawData', index = False) #sends dfs to their own sheet in the generated excel file
            worksheet = writer.sheets['Sample' + str(i+1) + 'RawData']
            worksheet.insert_image('K1', 'Plot' + str(i+1) + '.png')
        SSM.to_excel(writer, sheet_name='Report')
        reportsheet = writer.sheets['Report']
        reportsheet.insert_image('A6', 'Overlay.png')
        reportsheet.set_column(0,15,20)
   
    
    #Writes .erg file for LIMS upload
#    print('what is the test #?')
#    tn = input()
#    with open (tn+'.erg', mode = 'w') as csv_file:
#        fieldnames = ['Testnumber', 'E modulus_x', 'Stress at 0.5% Strain_x', 'Stress at 1.0% Strain_x', 'Stress at 2.0% Strain_x', 'Stress at 5% Strain_x', 'Stress at 10% Strain_x', 'Stress at 50% Strain_x', 'Stress at 100% Strain_x', 'Stress at 300% Strain_x', 'Yield stress_x', 'Stress at Break_x', 'Maximum Stress_x', 'Yield Strain_x', 'Strain at Break_x', 'Cross sectional area_x', 'Thickness_x', 'Width_x']
#        writer = csv.DictWriter(csv_file, fieldnames=fieldnames, delimiter = ';')
#        writer.writeheader()
#        writer.writerow({'Testnumber': tn, 'E modulus_x': SSM['AVG']['Modulus (MPa)'], 'Stress at 0.5% Strain_x': ' ', 'Stress at 1.0% Strain_x': ' ', 'Stress at 2.0% Strain_x': ' ', 'Stress at 5% Strain_x': ' ', 'Stress at 10% Strain_x': ' ', 'Stress at 50% Strain_x': ' ', 'Stress at 100% Strain_x': ' ', 'Stress at 300% Strain_x': ' ', 'Yield stress_x': ' ', 'Stress at Break_x': SSM['AVG']['Stress at Break (MPa)'], 'Maximum Stress_x': ' ', 'Yield Strain_x': ' ', 'Strain at Break_x': SSM['AVG']['Strain at Break (%)'], 'Cross sectional area_x': ' ', 'Thickness_x' : ' ', 'Width_x': ' ' })
#        
                     
consolidate('O:\\Tech_Service\\Analytical - New\\Analytical\\TAR Data\\2019\\19-0075 Spec. development GMMPD for Akulon K-FKG3 BK00001- Deans 3-19-2019\\20190524 Intertek 1572 Raw Data\\1572 Raw Data\\1572 Tensile -40°C.is_tens_RawData')
     
   

#DONE: check break cutoff location is good
#DONE: interpolate for modulus region
#DONE make plots start at origin
#DONE: remove tails at end of S-S curve.
#DONE: plot overlay
#DONE: Place specimen properties as 1st sheet  
#DONE: find max strain in iterative way
#Done: Column width formatting 