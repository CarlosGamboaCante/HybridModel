
# ===========================================================
# 📌 HYBRIDMODEL TOOL: HybridModelBalancing + HybridModelSources
# ===========================================================

# The first section of the code performs iterative calculations to display all combinations so that the distribution of the seismic moment ratepotential between the zone-type sources and faults-type sources is in balance

  HybridModelBalancing <- function(LYobs,Mmin,MmaxZone1,MmaxZone2,mu,dec) {
    
# Load seismic and faults catalogs from CSV files  

Seismic <- read.csv("SeismicCatalog.csv", header = TRUE, sep = ";")
Faults <- read.csv("FaultCatalog.csv", header = TRUE, sep = ";")

# Create a DataFrame with the model's input parameters 
InputData <- data.frame(LYobs = LYobs, Mmin = Mmin, MmaxZone1 = MmaxZone1, MmaxZone2 = MmaxZone2, mu = mu, dec = dec)

#Save the input parameters to a CSV file for reference or further use 
write.csv(InputData, "InputData.csv", row.names = FALSE)

# Creating Vectors of MmaxC, btf, btz, MmaxZone:

SeismicMaxReg <- Seismic[Seismic$n != 0,] # Filter the seismic catalog to exclude entries where 'n' is zero
MmaxC <- seq(from = round(Mmin + 1, 1), to = round(max(SeismicMaxReg$m), 1), by = 0.1) # Create a sequence for MmaxC, ranging from (Mmin + 1) to the maximum observed magnitude, with increments of 0.1  

vectorBt <- seq(from = 1, to = 3, by = dec) # Generate a vector for 'b-value' ranging from 1 to 3, with step size defined by 'dec'
btf <- vectorBt
btz <- vectorBt

MmaxZone <- seq(from = round(MmaxZone1, 1), to = round(MmaxZone2, 1), by = 0.1) # Create a sequence for MmaxZone, ranging from MmaxZone1 to MmaxZone2, with increments of 0.1  

# Define the structure of the output DataFrame

OutputHM<- data.frame(MmaxC=as.integer(),btf=as.integer(),btz=as.integer(),MmaxZone=as.integer())

# Step 1: Calculate the seismic moment (Mo) and moment rate (tn_Mo) for each magnitude bin (m) in the seismic catalog  

MoMmin <- 10^(16.1 + 1.5 * Mmin) # Seismic moment for the minimum magnitude (Mmin) - Hanks and Kanamori (1979)

Seismic$Mo_m <- 10^(16.1 + 1.5 * Seismic$m) # Seismic moment (Mo) for each magnitude bin (m) in the seismic catalog - Hanks and Kanamori (1979)
Seismic$tn <- Seismic$n / (LYobs - Seismic$CYm) # Seismic rate (tn) for each magnitude bin (m) during the completeness period  
Seismic$tn_Mo <- Seismic$tn * Seismic$Mo_m # Seismic moment rate (tn_Mo) for each magnitude bin (m)

# Step 2: Calculate the seismic moment (Mo) for faults in the range Mmin to MmaxFault  

# Fault seismic moment (Mof)
Faults$Mof <- Faults$slip_rate * Faults$Area * mu * 1000 * 10000000 # Brune (1968)

#  Seismic moment for MmaxFault
Faults$MoMmaxFault <- 10^(16.1 + 1.5 * (Faults$MmaxFault + 0.1)) #   Hanks and Kanamori (1979).  A bin (0.1) is added to ensure the maximum magnitude is reached.

# Loop through different MmaxC values to compute seismic parameters  
for (k in 1:length(MmaxC)[1]) {

  # Calculate the seismic moment for MmaxC
    MoMmaxC <- 10^(16.1 + 1.5 * round(MmaxC[k] + 0.1, 1)) # Hanks and Kanamori (1979). A bin (0.1) is added to ensure the maximum magnitude is reached.
  
  # Filter seismic events between Mmin and MmaxC  
    Seismic_Mmin_MmaxC <- subset(Seismic, Seismic$m >= Mmin & Seismic$m <= MmaxC[k] )
  
  # Step 3: Calculation of cumulative seismic rate Mmin (Nmin) and moment rate (tn_Mo) for the region (faults+zone) between Mmin - MmaxC recorded in the catalog
  Nmin_Mmin_MmaxC_Reg <- sum(Seismic_Mmin_MmaxC$tn)
  Mo_Mmin_MmaxC_Reg <- sum(Seismic_Mmin_MmaxC$tn_Mo)
  
  # Iterate over different  bt-value parameters for the Fault-sources (btf) 
  for (j in 1:length(btf)[1]) {
  
    # Step 4: Calculate the cumulative seismic rate Mmin (Nmin) and moment rate (tn_Mo) for each faults between Mmin - Mmax and Mmin - MmaxC
    
	# Fault cumulative seismic rate N (Mmin  - MmaxFault):
	Faults$NMmin_MmaxFault <- Faults$Mof * (((1.5 * log(10) - btf[j]) * (exp(-btf[j] * Mmin) - exp(-btf[j] * (Faults$MmaxFault + 0.1)))) / (btf[j] * (exp(-btf[j] * (Faults$MmaxFault + 0.1)) * Faults$MoMmaxFault - exp(-btf[j] * Mmin) * MoMmin))) # Anderson (1979). A bin (0.1) is added to ensure the maximum magnitude is reached.
    
	# Fault cumulative seismic rate N (MmaxC - MmaxFault):
    Faults$NMmaxC_MmaxFault = ifelse(Faults$MmaxFault <= round(MmaxC[k], 1), 
                             0, 
                             Faults$NMmin_MmaxFault * (((exp(-btf[j] * round(MmaxC[k], 1))) - (exp(-btf[j] * (Faults$MmaxFault + 0.1)))) / ((exp(-btf[j] * Mmin)) - (exp(-btf[j] * (Faults$MmaxFault + 0.1)))))) # Cosentino et al. (1977). A bin (0.1) is added to ensure the maximum magnitude is reached.
    
    # Fault cumulative seismic rate N (Mmin - MmaxC):
	Faults$NMin_MmaxC = Faults$NMmin_MmaxFault - Faults$NMmaxC_MmaxFault
    
	# Fault cumulative moment rate Mo (Mmin - MmaxC):
    Faults$MoMin_MmaxC <- ifelse(Faults$MmaxFault > round(MmaxC[k], 1), 
                              Faults$NMin_MmaxC * (btf[j] * (exp(-btf[j] * round(MmaxC[k] + 0.1, 1)) * MoMmaxC   - exp(-btf[j] * Mmin) * MoMmin)) / ((1.5 * log(10) - btf[j]) * (exp(-btf[j] * Mmin) - exp(-btf[j] * round(MmaxC[k] + 0.1, 1)))), 
                              Faults$NMin_MmaxC * (btf[j] * (exp(-btf[j] * (Faults$MmaxFault + 0.1)) * Faults$MoMmaxFault - exp(-btf[j] * Mmin) * MoMmin)) / ((1.5 * log(10) - btf[j]) * (exp(-btf[j] * Mmin) - exp(-btf[j] * (Faults$MmaxFault + 0.1))))) # Anderson (1979). A bin (0.1) is added to ensure the maximum magnitude is reached.
    
    
	# Step 5: Calculate the total cumulative seismic rate (Nmin) and moment rate (tn_Mo) for all faults between Mmin and MmaxC  
    
	Nmin_Mmin_MmaxC_Faults <- sum(Faults$NMin_MmaxC)
    Mo_Mmin_MmaxC_Faults <- sum(Faults$MoMin_MmaxC)
    
	# Step 6: Calculate the cumulative seismic rate (NMmin) and moment rate (tn_Mo) for the zone between Mmin and MmaxC
	
    Nmin_Mmin_MmaxC_Zone <- Nmin_Mmin_MmaxC_Reg - Nmin_Mmin_MmaxC_Faults
    Mo_Mmin_MmaxC_Zone <- Mo_Mmin_MmaxC_Reg - Mo_Mmin_MmaxC_Faults
    
    #Step 7: Calculate the theoretical cumulative seismic rate (Nmin) for the zone between Mmin and MmaxC  
    
	# Iterate over different bt-value parameters for the zone-source (btz) 
    for (i in 1:length(btz)[1]) {
      
      # Iterate over different values of MmaxZone
      for (h in 1:length(MmaxZone)[1]) {
        
        # Calculate the seismic moment for MmaxZone 
		MoMmaxZone <- 10^(16.1 + 1.5 * round(MmaxZone[h] + 0.1, 1)) #  Hanks and Kanamori (1979). A bin (0.1) is added to ensure the maximum magnitude is reached.
        
        # Calculate the theoretical cumulative seismic rate (Nmin) for the zone 
		Nmin_Mmin_MmaxC_ZoneT <- Mo_Mmin_MmaxC_Zone * ifelse(round(MmaxZone[h], 1) < round(MmaxC[k], 1), 
                                                          (((1.5 * log(10) - btz[i]) * (exp(-btz[i] * Mmin) - exp(-btz[i] * round(MmaxZone[h] + 0.1, 1)))) / (btz[i] * (exp(-btz[i] * round(MmaxZone[h] + 0.1, 1)) * MoMmaxZone - exp(-btz[i] * Mmin) * MoMmin))), 
                                                          (((1.5 * log(10) - btz[i]) * (exp(-btz[i] * Mmin) - exp(-btz[i] * round(MmaxC[k] + 0.1, 1))))      / (btz[i] * (exp(-btz[i] * round(MmaxC[k] + 0.1, 1)) * MoMmaxC           - exp(-btz[i] * Mmin) * MoMmin))))
        
        
		# Step 8: Acceptance tolerance and required conditions for obtaining a coherent result
        DifZone <- round(abs(Nmin_Mmin_MmaxC_ZoneT - Nmin_Mmin_MmaxC_Zone), 3)
        
		if (DifZone >= 0.001) next
        
        if (Nmin_Mmin_MmaxC_Zone < 0) next
        
        if (Mo_Mmin_MmaxC_Zone < 0) next
        
		# Save only coherent results that meet the defined conditions   
		
        OutputHM[nrow(OutputHM) +1,] <- c(MmaxC[k],btf[j],btz[i],MmaxZone[h])
        
      }
    }
    
  }
}
       
	# Save the results to CSV files  
		write.csv(Seismic, file = "Seismic_Data.csv", row.names = FALSE)
		write.csv(OutputHM, file = "Output_HM.csv", row.names = FALSE)
			
	# Return the final output data frame  
			
		return(OutputHM)

  }


# The second section takes the selected combination as input for hazard calculation and displays the graphical interface with the earthquake occurrence model for that combination and source

  HybridModelSources <- function(i) {
    
    library(readr)
    # Load input data from CSV files  
    Seismic<- read_delim("Seismic_Data.csv", ",",escape_double = FALSE, trim_ws = TRUE) # Seismic catalog
    Faults<- read_delim("FaultCatalog.csv", ";",escape_double = FALSE, trim_ws = TRUE) # Fault catalog
    Output_HM<- read_delim("Output_HM.csv", ",",escape_double = FALSE, trim_ws = TRUE) # Model output from previous step
    Input<- read_delim("InputData.csv", ",",escape_double = FALSE, trim_ws = TRUE) # Input parameters
    
    # Extract key input parameters  
    Mmin = as.numeric(Input [1,2])  # Minimum magnitude
    mu =  as.numeric(Input [1,5])  # Crustal rigidity modulus 
    
	# Extract values from Output_HM based on index i :
    MmaxC = as.numeric(Output_HM [i,1]) # maximum completeness magnitude 
    btf = as.numeric(Output_HM [i,2]) # bt-value for faults
    btz = as.numeric(Output_HM[i,3]) # bt-value for zone
    MmaxZone = as.numeric(Output_HM [i,4]) # Maximum magnitude for the zone
    
    # Calculate seismic moments using Hanks and Kanamori (1979)  
	MoMmin=10^(16.1+1.5*Mmin) # Seismic moment for Mmin
    MoMmaxC=10^(16.1+1.5*(MmaxC + 0.1)) # Seismic moment for MmaxC. Added 0.1 bin to ensure max magnitude is reached
    MoMmaxZone=10^(16.1+1.5*(MmaxZone + 0.1)) # SeismicMmaxZone. Added 0.1 bin to ensure max magnitude is reached

    #Calculate fault-sources parameters:
	
		# Calculate fault seismic moment    
		Faults$Mof <- Faults$slip_rate * Faults$Area * mu * 1000 * 10000000 # Brune (1968)
		
		# Calculate seismic Moment of the MmaxFault
		Faults$MoMmaxFault <- 10^(16.1 + 1.5 * (Faults$MmaxFault + 0.1)) #  Hanks and Kanamori (1979). Added 0.1 bin to ensure max magnitude is reached
		
		# Calculate cumulative seismic rate (Nmin) for faults
		Faults$NMmin_MmaxFault = Faults$Mof * ((1.5*log(10)-btf)*(exp(-btf*Mmin)-exp(-btf*(Faults$MmaxFault+0.1)))/(btf*(exp(-btf*(Faults$MmaxFault+0.1))*Faults$MoMmaxFault-exp(-btf*Mmin)*MoMmin)))
		
		# Calculate cumulative seismic rate for faults between MmaxC and MmaxFault  
		Faults$NMmaxC_MmaxFault = ifelse(Faults$MmaxFault <= MmaxC,
								 0, 
								 Faults$NMmin_MmaxFault * ((exp(-btf*MmaxC)-exp(-btf*(Faults$MmaxFault+0.1)))/(exp(-btf*Mmin)-exp(-btf*(Faults$MmaxFault+0.1)))))
		
		# Calculate cumulative seismic rate between Mmin and MmaxC  
		Faults$NMin_MmaxC = Faults$NMmin_MmaxFault - Faults$NMmaxC_MmaxFault
		
		# Calculate cumulative moment rate for faults between Mmin and MmaxC  
		Faults$MoMin_MmaxC = ifelse(Faults$MmaxFault > MmaxC, 
								 Faults$NMin_MmaxC*(btf*(exp(-btf*(MmaxC+0.1))*MoMmaxC-exp(-btf*Mmin)*MoMmin))/((1.5*log(10)-btf)*(exp(-btf*Mmin)-exp(-btf*(MmaxC+0.1)))), 
								 Faults$NMin_MmaxC * (btf*(exp(-btf*(Faults$MmaxFault+0.1))*Faults$MoMmaxFault-exp(-btf*Mmin)*MoMmin))/((1.5*log(10)-btf)*(exp(-btf*Mmin)-exp(-btf*(Faults$MmaxFault+0.1)))))
		
		# Calculate total seismic moment for faults
		Mo_Mmin_MmaxC_Faults <- sum(Faults$MoMin_MmaxC)
		
		# Create a data frame with fault characteristics and rates calculations 
		FaultGR <- data.frame(ID=Faults$ID_Fault,Name=Faults$Name_Fault,Mmax=Faults$MmaxFault,NMmin_Mmax=Faults$NMmin_MmaxFault,Beta=btf)
    
    
	# Seismicity recorded in the catalog (Region)
   
		# Initialize a source model for storing seismic rates 
		
		SourceModel <- data.frame(ID_Fault=Faults$ID_Fault)
		SourceModel <- data.frame(t(SourceModel[-1])) # Transpose for easier manipulation  
		colnames(SourceModel) <- FaultGR[, 1]  # Set fault ID as column names 
		
		# Determine the maximum magnitude from faults, seismic data, and the zone  
		maxSource <- max(max(Seismic$m),Faults$MmaxFault,MmaxZone)
		
		# Create a sequence of magnitude values for source and faults  
		VecMSource <- ((Mmin*10):(maxSource *10))/10
		VecFault <- Faults$ID_Fault
		
		# Loop through each fault and magnitude to calculate seismic rates
		for (j in 1:length(VecFault)[1]) {
		  for (k in 1:length(VecMSource)[1]) {
			
			SourceModel[k,j]= ifelse(VecMSource[k] <= Faults[j,5], Faults[j,8]  * ((exp(-btf*(VecMSource[k]))-exp(-btf*(Faults[j,5]+0.1)))/(exp(-btf*Mmin)-exp(-btf*(Faults[j,5]+0.1)))), "")
			
		  }
		}
		
		# Filter the seismic catalog for magnitudes between Mmin and MmaxC  
		Seismic_Mmin_MmaxC <- subset(Seismic, Seismic$m >= Mmin & Seismic$m <= MmaxC)
		Seismic_Mmin_MMax<- subset(Seismic, Seismic$m >= Mmin)
		
		# Calculate the cumulative seismic rate and the seismic moment rate for the region (recorded in the seismic catalog).
		Nmin_Mmin_MmaxC_Reg <- sum(Seismic_Mmin_MmaxC$tn)
		Mo_Mmin_MmaxC_Reg <- sum(Seismic_Mmin_MmaxC$tn_Mo)
		
    #Calculate Zone-Source parameters:
    
		# Calculate the seismic moment rate for the zone (difference between regional and fault moment rates)
		Mo_Mmin_MmaxC_Zone <- Mo_Mmin_MmaxC_Reg - Mo_Mmin_MmaxC_Faults
    
		# Calculate the cumulative seismic rate (Nmin) for the zone between Mmin and MmaxC
		Nmin_Mmin_MmaxC_Zone <- Mo_Mmin_MmaxC_Zone * ifelse(round(MmaxZone, 1) < round(MmaxC, 1), 
														  (((1.5 * log(10) - btz) * (exp(-btz * Mmin) - exp(-btz * round(MmaxZone + 0.1, 1)))) / (btz * (exp(-btz * round(MmaxZone + 0.1, 1)) * MoMmaxZone - exp(-btz * Mmin) * MoMmin))), 
														  (((1.5 * log(10) - btz) * (exp(-btz * Mmin) - exp(-btz * round(MmaxC + 0.1, 1))))    / (btz * (exp(-btz * round(MmaxC + 0.1, 1)) * MoMmaxC       - exp(-btz * Mmin) * MoMmin))))
		
		# Create a data frame for the zone-source seismic parameters
		Zone <- data.frame(ID="Z",Name="Zone",Mmax=MmaxZone,Beta=btz,Nmin_MmaxC= Nmin_Mmin_MmaxC_Zone)
		
		# Calculate the seismic rate contribution for the zone-source between MmaxC and MmaxZone
		Nmin_MmaxC_MMax_Zone = Nmin_Mmin_MmaxC_Zone  * ((exp(-btz*(MmaxZone + 0.1))-exp(-btz*(MmaxC)))/(exp(-btz*(Mmin))-exp(-btz*(MmaxC))))
		
		# Calculate the final seismic rate for the zone-source
		Zone$NMmin_Mmax=ifelse(MmaxC<MmaxZone,Nmin_Mmin_MmaxC_Zone - Nmin_MmaxC_MMax_Zone,Nmin_Mmin_MmaxC_Zone)
		
		# Create a summary table for the zone-sopurce seismic parameters
		ZoneGR <- data.frame(ID="Z",Name="Zone",Mmax=MmaxZone,NMmin_Mmax= Zone$NMmin_Mmax,Beta=btz)  
    
		# Assign seismic rate to the zone-source model based on magnitude values
		SourceModel$Zone = ifelse(VecMSource <= MmaxZone, Zone$NMmin_Mmax * ((exp(-btz*(VecMSource))-exp(-btz*(MmaxZone+0.1)))/(exp(-btz*(Mmin))-exp(-btz*(MmaxZone+0.1)))), "")
		
		# Combine fault-spurces and zone-source seismic parameters into one table
		SourceGR <- rbind(FaultGR,ZoneGR) # Merge all sources into a single dataset
		
		# Calculate Gutenberg-Richter "b" and "a" parameters for each seismic source
		SourceGR$b = SourceGR$Beta/log(10)
		SourceGR$a = log10(SourceGR$NMmin_Mmax) + SourceGR$b * Mmin
		
    # Plot the Gutenberg-Richter model for fault-sources and the zone-source. 
	
		# Plot fault-cources GR models
		matplot(VecMSource, SourceModel, xlab = "Magnitude", ylab = "N(m)", log = "y", type = "l", col = 1:(length(Faults) + 2), lty = 2, main = "Gutenberg-Richter Model")  
		
		# Add zone-source GR model
		matlines (VecMSource, SourceModel$Zone, type = "l", lwd=2, col = 2)

    # Save results to CSV files for further analysis

    write.csv(SourceGR, "SourceGR.csv") 
    write.csv(SourceModel, "SourceModel.csv")
    
  }
  
  
