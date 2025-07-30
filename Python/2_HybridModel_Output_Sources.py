
# ========================================
# ?? HYBRIDMODEL TOOL: HybridModelSources
# ========================================

import pandas as pd
import numpy as np
import math
import os
import tkinter as tk
from tkinter import messagebox
from scipy.optimize import minimize_scalar

# This section of the code takes the selected combination as input for hazard calculation and displays the graphical interface with the earthquake occurrence model for that combination and source

def HybridModelSources(i):
    # Load input data from CSV files
    Seismic = pd.read_csv("Seismic_Data.csv", delimiter=",")
    Faults = pd.read_csv("FaultCatalog.csv", delimiter=",")
    Output_HM = pd.read_csv("Output_HM.csv", delimiter=",")
    Input = pd.read_csv("InputData.csv", delimiter=",")

    # Extraction key input parameters
    Mmin = float(Input.iloc[0, 1]) # Minimum magnitude
    mu = float(Input.iloc[0, 4]) #Crustal rigidity modulus
	
	#Extract values from Output_HM based on index i:
    MmaxC = float(Output_HM.iloc[i-1, 1]) #Maximum magnitude of completeness
    btf = float(Output_HM.iloc[i-1, 2]) #'?-value' for fault-type sources
    btz = float(Output_HM.iloc[i-1, 3]) #'?-value' for zone-type source
    MmaxZone = float(Output_HM.iloc[i-1, 4]) #Maximum magnitude of zone-type source
	
	#Calculate seismic moment (Mo) using Hanks and Kanamori (1979)
    MoMmin = 10 ** (16.1 + 1.5 * Mmin) #For Minimum magnitude
    MoMmaxC = 10 ** (16.1 + 1.5 * (MmaxC+0.1)) #For Maximum magnitude of completeness. Added 0.1 bin to ensure MmaxC is reached.
    MoMmaxZone = 10 ** (16.1 + 1.5 * (MmaxZone+0.1)) #For Maximum magnitude of zone-type source. Added 0.1 bin to ensure MmaxZone is reached.
	
	#Calculate fault-type sources parameters:
	
	#Calculate fault seismic moment (Mof)
    Faults['Mof'] = Faults['slip_rate'] * Faults['Area'] * mu * 1000 * 10000000 #Brune (1968)
	
	#Calculate seismic moment (Mo) of the MmaxFaults
    Faults['MoMaxFault'] = 10**(16.1 + 1.5 * (Faults['MmaxFault']+0.1)) #Hanks and Kanamori (1979). Added 0.1 bin to ensure MmaxFault is reached.
    
    ID_Selected = pd.DataFrame({'ID':[i],'Mmin': [Mmin], 'MmaxC': [MmaxC], 'btf': [btf], 'btz': [btz], 'Mmax_Zone': [MmaxZone]})
    
	#Calculate cumulative seismic rate (Nmin) for fault-type sources
    Faults['NMmin_MmaxFault'] = Faults['Mof'] * (((1.5 * math.log(10) - btf) * (np.exp(-btf * Mmin) - np.exp(-btf * (Faults['MmaxFault']+0.1))))/(btf * (np.exp(-btf * (Faults['MmaxFault']+0.1)) * Faults['MoMaxFault'] - np.exp(-btf * Mmin) * MoMmin))) #Anderson (1979)
	
	#Calculate cumulative seismic rate for fault-type sources between MmaxC and MmaxFault 
    Faults['NMmaxC_MmaxFault'] = np.where(Faults['MmaxFault'] <= round(MmaxC,1), 0, Faults['NMmin_MmaxFault'] * (((np.exp(-btf * round(MmaxC,1)))-(np.exp(-btf * (Faults['MmaxFault']+0.1))))/((np.exp(-btf * Mmin))-(np.exp(-btf * (Faults['MmaxFault']+0.1)))))) #Cosentino (1977)
	
	#Calculate cumulative seismic rate between Mmin and MmaxC for fault-type sources
    Faults['NMin_MmaxC'] = Faults['NMmin_MmaxFault'] - Faults['NMmaxC_MmaxFault']
	
	#Calculate cumulative moment rate (?o) for fault-type sources between Mmin and MmaxC 
    Faults['MoMin_MmaxC'] = np.where(Faults['MmaxFault'] > round(MmaxC,1), Faults['NMin_MmaxC'] * (btf * (np.exp(-btf * round(MmaxC+0.1,1)) * MoMmaxC - np.exp(-btf * Mmin) * MoMmin))/((1.5 * math.log(10) - btf) * (np.exp(-btf * Mmin) - np.exp(-btf * round(MmaxC+0.1,1)))), Faults['NMin_MmaxC'] * (btf * (np.exp(-btf * (Faults['MmaxFault']+0.1)) * Faults['MoMaxFault'] - np.exp(-btf * Mmin) * MoMmin))/((1.5 * math.log(10) - btf) * (np.exp(-btf * Mmin) - np.exp(-btf * (Faults['MmaxFault']+0.1))))) #Anderson (1979)
	
	#Calculate total seismic moment (?o) for fault-type sources
    Mom_Mmin_MmaxC_Faults = Faults['MoMin_MmaxC'].sum()
	
	#Create a data frame with fault characteristics and rates calculations
    FaultsGR = pd.DataFrame({'ID': Faults['ID_Fault'], 'Name': Faults['Name_Fault'], 'Mmax': Faults['MmaxFault'], 'NMmin_Mmax': Faults['NMmin_MmaxFault'], 'Beta': btf})

	#Seismicity recorded in the catalog (Region)
	
	#Initialize a source model for storing seismic rates
    SourceModel = pd.DataFrame({'ID_Fault': Faults['ID_Fault']})
    SourceModel = SourceModel.transpose() #Transpose for easier manipulation
    SourceModel.columns = FaultsGR.iloc[:, 0]
    SourceModel = SourceModel.drop('ID_Fault') #Set fault ID as column names
	
	#Determine the higher maximum magnitude between fault-type sources, seismic data, and the zone-type source
    maxSource = max(max(Seismic['m']), Faults['MmaxFault'].max(), MmaxZone)
	
	#Create a sequence of magnitude values for source and faults
    VecMSource = np.arange(round(Mmin,1),round(maxSource+0.1,1),0.1)
    VecFault = Faults['ID_Fault']
	
	#Loop through each fault and magnitude to calculate seismic rates
    for j in range(len(VecFault)):
        for k in range(len(VecMSource)):
            if VecMSource[k] < Faults.iloc[j,4]:
                SourceModel.loc[k, j] = Faults.iloc[j, 7] * ((np.exp(-btf * VecMSource[k]) - np.exp(-btf * (Faults.iloc[j, 4]+0.1))) / (np.exp(-btf * Mmin) - np.exp(-btf * (Faults.iloc[j, 4]+0.1))))
            else:
                SourceModel.loc[k, j] = ""
    for j in range(len(VecFault)):
        for k in range(len(VecMSource)):
            if VecMSource[k] < Faults.iloc[j,4]:
                SourceModel.iloc[k, j] = Faults.iloc[j, 7] * ((np.exp(-btf * VecMSource[k]) - np.exp(-btf * (Faults.iloc[j, 4]+0.1))) / (np.exp(-btf * Mmin) - np.exp(-btf * (Faults.iloc[j, 4]+0.1))))
            else:
                SourceModel.iloc[k, j] = ""
    SourceModel = SourceModel.drop(range(len(VecFault)), axis=1)
	
	#Filter the seismic catalog for magnitudes between Mmin and MmaxC
    Seismic_Mmin_MmaxC = Seismic[(Seismic['m'] >= Mmin) & (Seismic['m'] <= MmaxC)]
    Seismic_Mmin_MMax = Seismic[Seismic['m'] >= Mmin]
	
	#Calculate the cumulative seismic rate and the seismic moment rate for the region
    Nmin_Mmin_MmaxC_Reg = Seismic_Mmin_MmaxC['tn'].sum()
    Mom_Mmin_MmaxC_Reg = Seismic_Mmin_MmaxC['tn_Mo'].sum()

    #Calculate Zone-type source parameters:
	
	#Calculate the seismic moment rate for the zone (difference between regional and fault moment rates)
    Mom_Mmin_MmaxC_Zone = Mom_Mmin_MmaxC_Reg - Mom_Mmin_MmaxC_Faults
	
	#Calculate the cumulative seismic rate (Nmin) for the zone-type source between Mmin and MmaxC
    Nmin_Mmin_MmaxC_Zone = np.where(MmaxZone < round(MmaxC,1),Mom_Mmin_MmaxC_Zone * (((1.5 * math.log(10) - btz) * (np.exp(-btz * Mmin) - np.exp(-btz * (MmaxZone+0.1))))/(btz * (np.exp(-btz * (MmaxZone+0.1)) * MoMmaxZone - np.exp(-btz * Mmin) * MoMmin))), Mom_Mmin_MmaxC_Zone * (((1.5 * math.log(10) - btz) * (np.exp(-btz * Mmin) - np.exp(-btz * (MmaxC+0.1))))/(btz * (np.exp(-btz * (MmaxC+0.1)) * MoMmaxC - np.exp(-btz * Mmin) * MoMmin)))) #Anderson (1979)
    
	#Create a data frame for the zone-type source seismic parameters
    Zone = pd.DataFrame({'ID': ['Z'], 'Name': ['Zone'], 'Mmax': [MmaxZone], 'Beta': [btz], 'Nmin_MmaxC': [Nmin_Mmin_MmaxC_Zone]})
	
	#Calculate the seismic rate contribution for the zone-type source between MmaxC and MmaxZone
    Nmin_MmaxC_MmaxZone = Nmin_Mmin_MmaxC_Zone * ((np.exp(-btz * (MmaxZone + 0.1)) - np.exp(-btz * MmaxC))/(np.exp(-btz * Mmin) - np.exp(-btz * MmaxC))) #Cosentino (1977)
	
	#Calculate the final seismic rate for the zone-type source between Mmin and MmaxZone
    Zone['NMmin_MmaxZone'] = np.where(MmaxC < MmaxZone, Nmin_Mmin_MmaxC_Zone - Nmin_MmaxC_MmaxZone, Nmin_Mmin_MmaxC_Zone)
	
	#Create a summary table for the zone-type source seismic parameters
    ZoneGR = pd.DataFrame({'ID': ['Z'], 'Name': ['Zone'], 'Mmax': [MmaxZone], 'NMmin_Mmax': Zone['NMmin_MmaxZone'], 'Beta': btz})

	#Assign seismic rate to the zone-type source model based on magnitude values
    SourceModelZone = pd.DataFrame({'Zone':['Zone']})
    VecMZone = np.arange(round(Mmin,1),round(MmaxZone+0.1,1),0.1)
    VecZone = SourceModelZone.transpose()
    for l in range(len(VecZone)):
        for m in range(len(VecMZone)):
            if VecMSource[m] < MmaxZone:
                SourceModelZone.loc[m, l] = ZoneGR.iloc[0,3] * ((np.exp(-btz * VecMSource[m]) - np.exp(-btz * (MmaxZone+0.1))) / (np.exp(-btz * Mmin) - np.exp(-btz * (MmaxZone+0.1))))
            else:
                SourceModelZone.loc[m, l] = ""
    for l in range(len(VecZone)):
        for m in range(len(VecMZone)):
            if VecMSource[m] < MmaxZone:
                SourceModelZone.iloc[m, l] = ZoneGR.iloc[0,3] * ((np.exp(-btz * VecMSource[m]) - np.exp(-btz * (MmaxZone+0.1))) / (np.exp(-btz * Mmin) - np.exp(-btz * (MmaxZone+0.1))))
            else:
                SourceModelZone.iloc[m, l] = ""
    SourceModelZone = SourceModelZone.drop(range(len(VecZone)), axis=1)
    
	#Creating a serie with the maximum bin of magnitude
    VecMag = pd.Series(np.round(np.arange(round(Mmin,1), round(maxSource+0.1,1),0.1), 1), name="m")
	
    #Combine fault-type sources and zone-type source seismic rates into one table
    SourceGR = pd.concat([FaultsGR, ZoneGR]) # Merge all sources into a single dataset

    #Calculate Gutenberg-Richter "b-value" and "a-value" parameters for each seismic source
    SourceGR['b'] = SourceGR['Beta'] / np.log(10)
    SourceGR['a'] = np.log10(SourceGR['NMmin_Mmax']) + SourceGR['b'] * Mmin

    #Creating the data frame for graphic the seismic sources    
    Final_SourceModel = pd.concat([VecMag,SourceModelZone,SourceModel],axis=1)
    
	# Graphic of the results according to selected ID
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    Graphic_SourceModel = Final_SourceModel.copy()
    Graphic_SourceModel = Graphic_SourceModel.replace('','NaN')
    x_values=Graphic_SourceModel.iloc[:,0]
    y_values = Graphic_SourceModel.iloc[:, 1:]
    plt.figure(figsize=(5.91, 5.12))#inches
    for column in y_values.columns:
        plt.plot(x_values, y_values[column], label=column)

    plt.xlabel('Magnitude (Mw)', fontsize=12, fontweight='bold', fontname='Times New Roman')
    plt.ylabel('Cummulative seismic rate Ṅ(m)', fontsize=12, fontweight='bold', fontname='Times New Roman')
    plt.yscale('log')
    plt.title('MFD GR modified', fontsize=18, fontweight='bold', fontname='Times New Roman')
    plt.xticks(fontsize=10, fontname='Times New Roman')
    plt.yticks(fontsize=10, fontname='Times New Roman')
    ax = plt.gca()
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, subs=[1.0], numticks=10))
    ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=10))
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    ax.grid(which='major', linestyle='-', linewidth=0.2, color='darkgray')
    ax.grid(which='minor', linestyle='--', linewidth=0.2, color='darkgray')
    plt.legend()
	

    # Saving results
    path_actual=os.path.dirname(os.path.abspath(__file__))
    file_name = "\RESULT_ID_" + str(i)
    complete_path = path_actual + file_name
	
    if not os.path.exists(complete_path):
     try:
        # Create the directory if it doesn't exist
        os.mkdir(complete_path)
        print(f"Folder '{file_name}' created successfully in {complete_path}")
     except Exception as e:
        print(f"Error to created the folder: {e}")
     

    # Specify the filename for the saved figure
    file_name = "GR_plot.png"

    # Combine the folder path and the filename
    file_path = os.path.join(complete_path, file_name)
    plt.savefig(file_path)

	
    Final_SourceModel.to_csv(complete_path + "/SourceModel.csv", index=False)
    SourceGR.to_csv(complete_path + "/SourceGR.csv", index=False)
    ID_Selected.to_csv(complete_path + "/ID_Selected.csv", index=False)

#Create GUI window
window = tk.Tk()
window.title("HybridModel Tool Output")

label = tk.Label(window, text="Enter ID number from Output_HM.csv:")
label.pack(padx=10, pady=10)

entry_id = tk.Entry(window)
entry_id.pack(padx=10, pady=10)

def run_model():
    try:
        id_value = int(entry_id.get())
        HybridModelSources(id_value)
    except ValueError:
        messagebox.showerror("Invalid input", "Please enter a valid integer ID.")

button_run = tk.Button(window, text="Run HybridModel Output", command=run_model)
button_run.pack(padx=10, pady=20)

window.mainloop()