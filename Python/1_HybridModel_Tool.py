
# =========================================
# 📌 HYBRIDMODEL TOOL: HybridModelBalancing
# =========================================

# The code performs iterative calculations to display all combinations so that the distribution of the seismic moment rate potential between the zone-type and faults-type sources is in balance

import pandas as pd
import numpy as np
import math
import tkinter as tk
from tkinter import filedialog, ttk, messagebox
import threading
import os
import time

def HybridModelBalancing(LYobs, Mmin, MmaxZone1, MmaxZone2, mu, dec, seismic_file, fault_file, progress_bar):
    try:
        # Load seismic and faults catalog from CSV files
        Seismic = pd.read_csv(seismic_file, delimiter=",")
        Faults = pd.read_csv(fault_file, delimiter=",")

        # Create a DataFrame with the model's input parameters 
        InputData = pd.DataFrame({'LYobs': [LYobs], 'Mmin': [Mmin], 'MmaxZone1': [MmaxZone1], 'MmaxZone2': [MmaxZone2], 'mu': [mu], 'dec': [dec]})
        
		#Save the input parameters to a CSV file for reference or further use
        InputData.to_csv("InputData.csv", index=False)

        # Creating Vectors of MmaxC, btf, btz, MmaxZone:
        SeismicMaxReg = Seismic[Seismic['n'] != 0] # Filter the seismic catalog to exclude entries where 'n' is zero
        MmaxC = np.arange(round(Mmin+1,1),round(SeismicMaxReg['m'].max()+0.1,1), 0.1) # Create a sequence for MmaxC, ranging from (Mmin + 1.0) to the maximum observed magnitude, with increments of 0.1

        vectorBt = np.arange(1,3+dec,dec) # Generate a vector for 'β-value' ranging from 1.0 to 3.0, with step size defined by 'dec'
        btf = vectorBt
        btz = vectorBt

        MmaxZone = np.arange(round(MmaxZone1,1),round(MmaxZone2+0.1, 1), 0.1) # Create a sequence for MmaxZone, ranging from MmaxZone1 to MmaxZone2, with increments of 0.1
		
		# Define the structure of the output DataFrame
        OutputHM = pd.DataFrame(columns=['ID','MmaxC', 'btf', 'btz', 'MmaxZone'])

        # ========================================================================================================================
		# Step 1: Calculation of seismic moment (Mo) and moment rate (tn_Mo) for each bin of magnitude (m) in the seismic catalog
		# ========================================================================================================================
		
        MoMmin = 10**(16.1 + 1.5 * Mmin) # Seismic moment for the minimum magnitude (Mmin) - Hanks and Kanamori (1979)
        Seismic['Mo_m'] = 10**(16.1 + 1.5 * Seismic['m']) # Seismic moment of each magnitude bin (m) - Hanks and Kanamori (1979)
        Seismic['tn'] = Seismic['n'] / (LYobs - Seismic['CYm']) # Seismic rate (tn) for each magnitude bin (m) during the completeness period
        Seismic['tn_Mo'] = Seismic['tn'] * Seismic['Mo_m'] # Seismic moment rate (tn_Mo) for each magnitude bin (m)

        # ==================================================================================
		# Step 2:Calculate the seismic moment (Mo) and moment rate (Mof) for faults in the range Mmin to MmaxFault 
		# ==================================================================================
		
        # Fault seismic moment rate (Mof)
        Faults['Mof'] = Faults['slip_rate'] * Faults['Area'] * mu * 1000 * 10000000 #Brune (1968)
        # Seismic moment for MmaxFault
        Faults['MoMmaxFault'] = 10**(16.1 + 1.5 * (Faults['MmaxFault']+0.1)) #Hanks and Kanamori (1979). A bin (0.1) is added to ensure the maximum magnitude is reached.
        
		#Define loop length
        id = 1
        total_iterations = len(MmaxC) * len(btf) * len(btz) * len(MmaxZone)
        current_iteration = 0
        
		#Loop through different MmaxC values to compute seismic parameters: btf, btz, MmaxZone
        for k in range(len(MmaxC)):
            
			#Calculate the seismic moment for MmaxC
            MoMmaxC = 10**(16.1 + 1.5 * round(MmaxC[k]+0.1,1)) #Hanks and Kanamori (1979). A bin (0.1) is added to ensure the maximum magnitude is reached.

            #Filter seismic events between Mmin and MmaxC
            Seismic_Mmin_MmaxC = Seismic[(Seismic['m'] >= Mmin) & (Seismic['m'] <= round(MmaxC[k],1))]

            # =============================================================================================================================================================
			# Step 3: Calculation of cumulative seismic rate Mmin (Nmin) and moment rate (tn_Mo) for the region (faults+zone) between Mmin - MmaxC recorded in the catalog
			# =============================================================================================================================================================
            Nmin_Mmin_MmaxC_Reg = Seismic_Mmin_MmaxC['tn'].sum()
            Mo_Mmin_MmaxC_Reg = Seismic_Mmin_MmaxC['tn_Mo'].sum()
			
			#Iterate over different 'β-value' parameters for the Fault-type sources (btf)
            for j in range(len(btf)):
                
				# ==========================================================================================================================
				# Step 4: Calculate the cumulative seismic rate [N(m)] and moment rate (Mo) for each fault-type source between Mmin - MmaxC
				# ==========================================================================================================================
				
				#Cumulative seismic rate [N(m)] between Mmin  - MmaxFault for each fault-type source:
                Faults['NMmin_MmaxFault'] = Faults['Mof'] * (((1.5 * math.log(10) - btf[j]) * (np.exp(-btf[j] * Mmin) - np.exp(-btf[j] * (Faults['MmaxFault']+0.1))))/(btf[j] * (np.exp(-btf[j] * (Faults['MmaxFault']+0.1)) * Faults['MoMmaxFault'] - np.exp(-btf[j] * Mmin) * MoMmin))) # Anderson (1979). A bin (0.1) is added to ensure the maximum magnitude is reached.
				
				#Cumulative seismic rate [N(m)] between MmaxC  - MmaxFault for each fault-type source:
                Faults['NMmaxC_MmaxFault'] = np.where(Faults['MmaxFault'] <= round(MmaxC[k],1), 0, Faults['NMmin_MmaxFault'] * (((np.exp(-btf[j] * round(MmaxC[k],1)))-(np.exp(-btf[j] * (Faults['MmaxFault']+0.1))))/((np.exp(-btf[j] * Mmin))-(np.exp(-btf[j] * (Faults['MmaxFault']+0.1)))))) # Cosentino et al. (1977). A bin (0.1) is added to ensure the maximum magnitude is reached.
				
				#Cumulative seismic rate [N(m)] between Mmin  - MmaxC for each fault-type source:
                Faults['NMin_MmaxC'] = Faults['NMmin_MmaxFault'] - Faults['NMmaxC_MmaxFault']
				
				#Cumulative moment rate (Mo) between Mmin  - MmaxC for each fault-type source:
                Faults['MoMin_MmaxC'] = np.where(Faults['MmaxFault'] > round(MmaxC[k],1), Faults['NMin_MmaxC'] * (btf[j] * (np.exp(-btf[j] * round(MmaxC[k]+0.1,1)) * MoMmaxC - np.exp(-btf[j] * Mmin) * MoMmin))/((1.5 * math.log(10) - btf[j]) * (np.exp(-btf[j] * Mmin) - np.exp(-btf[j] * round(MmaxC[k]+0.1,1)))), Faults['NMin_MmaxC'] * (btf[j] * (np.exp(-btf[j] * (Faults['MmaxFault']+0.1)) * Faults['MoMmaxFault'] - np.exp(-btf[j] * Mmin) * MoMmin))/((1.5 * math.log(10) - btf[j]) * (np.exp(-btf[j] * Mmin) - np.exp(-btf[j] * (Faults['MmaxFault']+0.1))))) # Anderson (1979). A bin (0.1) is added to ensure the maximum magnitude is reached.

                # ===============================================================================================================================
				# Step 5: Calculate the total cumulative seismic rate (Nmin) and moment rate (Mo) for all fault-type sources between Mmin - MmaxC
				# ===============================================================================================================================
				
                Nmin_Mmin_MmaxC_Faults = Faults['NMin_MmaxC'].sum()
                Mo_Mmin_MmaxC_Faults = Faults['MoMin_MmaxC'].sum()

                # ========================================================================================================================
				# Step 6: Calculate the cumulative seismic rate (Nmin) and moment rate (Mo) for the zone-type source between Mmin - MmaxC
				# ========================================================================================================================
				
                Nmin_Mmin_MmaxC_Zone = Nmin_Mmin_MmaxC_Reg - Nmin_Mmin_MmaxC_Faults
                Mo_Mmin_MmaxC_Zone = Mo_Mmin_MmaxC_Reg - Mo_Mmin_MmaxC_Faults
				
				# ==================================================================================================================================================================================
				# Step 7: Calculate the theoretical cumulative seismic rate (Nmin theoretical) for the zone-type source between Mmin - MmaxZone, which MmaxZone take values between MmaxZone1 - MmaxZone2
				# ==================================================================================================================================================================================
				
                # Iterate over different 'β-value' parameters for the zone-type source (btz)
                for i in range(len(btz)):
                    
					# Iterate over different values of MmaxZone
                    for h in range(len(MmaxZone)):
					    
						# Calculate the seismic moment (Mo) for MmaxZone 
                        MoMmaxZone = 10**(16.1 + 1.5 * round(MmaxZone[h]+0.1,1)) #  Hanks and Kanamori (1979). A bin (0.1) is added to ensure the maximum magnitude is reached.
						
						# Calculate the theoretical cumulative seismic rate (Nmin theoretical) for the zone-type source 
                        Nmin_Mmin_MmaxC_ZoneT = Mo_Mmin_MmaxC_Zone * np.where(round(MmaxZone[h],1) < round(MmaxC[k],1), (((1.5 * math.log(10) - btz[i]) * (np.exp(-btz[i] * Mmin) - np.exp(-btz[i] * round(MmaxZone[h]+0.1,1))))/(btz[i] * (np.exp(-btz[i] * round(MmaxZone[h]+0.1,1)) * MoMmaxZone - np.exp(-btz[i] * Mmin) * MoMmin))), (((1.5 * math.log(10) - btz[i]) * (np.exp(-btz[i] * Mmin) - np.exp(-btz[i] * round(MmaxC[k]+0.1,1))))/(btz[i] * (np.exp(-btz[i] * round(MmaxC[k]+0.1,1)) * MoMmaxC - np.exp(-btz[i] * Mmin) * MoMmin)))) # Anderson (1979). A bin (0.1) is added to ensure the maximum magnitude is reached.
						
						# ====================================================================================
						# Step 8: Acceptance tolerance and required conditions for obtaining a coherent result
						# ====================================================================================
						
                        DifZone = round(abs(Nmin_Mmin_MmaxC_ZoneT - Nmin_Mmin_MmaxC_Zone),3)

                        if DifZone >= 0.001 or Nmin_Mmin_MmaxC_Zone < 0 or Mo_Mmin_MmaxC_Zone < 0:
                            continue
                        else: 
                            # Save only coherent results that meet the defined conditions 
                            OutputHM = OutputHM.append({'ID':id,'MmaxC': round(MmaxC[k],1), 'btf': btf[j], 'btz': btz[i], 'MmaxZone': round(MmaxZone[h],1)}, ignore_index=True)
                            id = id + 1

                        current_iteration += 1
                        progress = int((current_iteration / total_iterations) * 100)
                        progress_bar["value"] = progress
                        window.update_idletasks()

        # Save the results to CSV files 
        Seismic.to_csv("Seismic_Data.csv", index=False)
        OutputHM.to_csv("Output_HM.csv", index=False)

        messagebox.showinfo("Success", "Hybrid model executed successfully. Results saved.")
    except Exception as e:
        messagebox.showerror("Error", f"An error occurred: {str(e)}")

def select_seismic_file():
    global lbl_seismic_file
    seismic_file_path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
    if seismic_file_path:
        seismic_file_name = os.path.basename(seismic_file_path)
        lbl_seismic_file.config(text=f"Seismic file: {seismic_file_name}")

def select_fault_file():
    global lbl_fault_file
    fault_file_path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
    if fault_file_path:
        fault_file_name = os.path.basename(fault_file_path)
        lbl_fault_file.config(text=f"Faults file: {fault_file_name}")

def run_hybrid_model():
    global progress_bar
    try:
        seismic_file_path = lbl_seismic_file.cget("text").split(": ")[1]
        fault_file_path = lbl_fault_file.cget("text").split(": ")[1]

        if not seismic_file_path or not fault_file_path:
            messagebox.showwarning("Warning", "Please select both Seismic and Faults files.")
            return

        LYobs = float(entry_CY.get())
        Mmin = float(entry_Mmin.get())
        MmaxZone1 = float(entry_MmaxZ1.get())
        MmaxZone2 = float(entry_MmaxZ2.get())
        mu = float(entry_n.get())
        dec = float(entry_dec.get())

        # Disable the button while processing
        btn_run_hybrid_model.config(state="disabled")

        # Run HybridModel function in a separate thread
        thread = threading.Thread(target=HybridModelBalancing, args=(LYobs, Mmin, MmaxZone1, MmaxZone2, mu, dec, seismic_file_path, fault_file_path, progress_bar))
        thread.start()

        # Update progress bar periodically
        while thread.is_alive():
            window.update()
            progress_bar.update_idletasks()
            progress_bar.step(10)  # Adjust the step value as needed
            time.sleep(0.1)  # Adjust sleep time as needed

        # Re-enable the button after processing
        btn_run_hybrid_model.config(state="normal")

    except Exception as e:
        messagebox.showerror("Error", f"An error occurred: {str(e)}")

def show_help():
    help_text = (
        "Parameters explanation:\n"
		"\n"
        "LYobs: Last year of registration in the seismic catalog.\n"
        "Mmin: Minimum magnitude for seismic distribution (Mw).\n"
        "MmaxZone1: Lower limit for the maximum magnitude of the zone (Mw).\n"
        "MmaxZone2: Upper limit for the maximum magnitude of the zone (Mw).\n"
        "mu: rigidity modulus (Pa) (it is recommended to use 3.00E+10).\n"
        "dec: precision for the b-value (NOTE: for small values of precision the calculation times will take longer).\n\n"
        "\n"
		"CSV file format: delimited by coma (,)\n"
		"\n"
        "Seismic: file with header m, CYm, n \n"
		"Where:\n"
		"m = each magnitude bin from Mmin to the maximum magnitude observed in the catalog (MmaxObs) \n"
		"CYm = the year of completeness of the seismic catalog for each magnitude bin \n"
		"n = number of earthquakes for each magnitude bin in the completeness period \n"
		"\n"
        "Faults: file with header ID_Fault, Name_Fault, slip_rate, Area, MmaxFault \n"
		"Where:\n"
		"ID_Fault = Numeric fault identifier \n"
		"Name_Fault = Name of the fault segment \n"
		"Area = Failure segment area (Km^2) \n"
		"slip_rate = average fault slip rate (mm/yr) \n"
		"MmaxFault = maximum magnitude that the fault segment can generate (Mw) (NOTE: a magnitude scaling relationship must be applied) \n"
    )

    # Create the help window
    help_window = tk.Toplevel()
    help_window.title("Help")
    
    # Add a label with the help information
    help_label = tk.Label(help_window, text=help_text, justify="left", padx=10, pady=10)
    help_label.pack()

# Create GUI window
window = tk.Tk()
window.title("HybridModel Tool")

# Create input fields and labels
labels = ["LYobs", "Mmin", "MmaxZone1", "MmaxZone2", "mu", "dec"]
for i, label in enumerate(labels):
    tk.Label(window, text=label).grid(row=i, column=0, padx=5, pady=5)
    entry = tk.Entry(window)
    entry.grid(row=i, column=1, padx=5, pady=5)
    if label == "LYobs":
        entry_CY = entry
    elif label == "Mmin":
        entry_Mmin = entry
    elif label == "MmaxZone1":
        entry_MmaxZ1 = entry
    elif label == "MmaxZone2":
        entry_MmaxZ2 = entry
    elif label == "mu":
        entry_n = entry
    elif label == "dec":
        entry_dec = entry

# Progress bar
progress_bar = ttk.Progressbar(window, orient="horizontal", length=300, mode="determinate")
progress_bar.grid(row=len(labels)+1, column=0, columnspan=2, pady=10)

# Buttons to select files
btn_select_seismic = tk.Button(window, text="Select Seismic Catalog File", command=select_seismic_file)
btn_select_seismic.grid(row=len(labels)+2, column=0, padx=5, pady=5)

btn_select_fault = tk.Button(window, text="Select Faults Catalog File", command=select_fault_file)
btn_select_fault.grid(row=len(labels)+2, column=1, padx=5, pady=5)

# Labels to display selected files
lbl_seismic_file = tk.Label(window, text="Seismic file: ")
lbl_seismic_file.grid(row=len(labels)+3, column=0, padx=5, pady=5)

lbl_fault_file = tk.Label(window, text="Faults file: ")
lbl_fault_file.grid(row=len(labels)+3, column=1, padx=5, pady=5)

# Button to run Hybrid Model
btn_run_hybrid_model = tk.Button(window, text="Run HybridModel", command=run_hybrid_model)
btn_run_hybrid_model.grid(row=len(labels)+4, column=0, columnspan=2, padx=5, pady=5)

# Button to show the Help 
btn_help = tk.Button(window, text="Help", command=show_help)
btn_help.grid(row=len(labels)+6, column=0, columnspan=2, padx=5, pady=5)

# Run GUI
window.mainloop()