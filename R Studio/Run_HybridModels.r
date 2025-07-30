
# ===========================================================
# 📌 EXECUTING THE CODE: RUNNING HYBRID MODELS
# ===========================================================

# Set working directory for Guatemala (update this path to match your local folder)

setwd("C:/Users/.../Guatemala/")  

# Run the `HybridModelBalancing` function with the input parameters:
# - Year: 2020
# - Minimum magnitude (Mmin): 4.5
# - Maximum magnitudes (MmaxZone1, MmaxZone2): 5.0, 6.5
# - Crustal rigidity modulus (mu): 3.00E+10
# - Increment step for bt-value (dec): 0.02

HybridModelBalancing(2020,4.5,5.0,6.5,3.00E+10,0.02)

# Run the `HybridModelSources` function for each dataset (one at a time)
HybridModelSources(1)
HybridModelSources(2)
HybridModelSources(3)
HybridModelSources(4)
HybridModelSources(5)
HybridModelSources(6)
HybridModelSources(7)
HybridModelSources(8)

# Set working directory for Granada (update this path to match your local folder)

setwd("C:/Users/.../Granada/") 

# Run the same process for the Granada study region

HybridModelBalancing(2023, 4.0, 4.5, 5.5, 3.00E+10, 0.2)  

# Run the `HybridModelSources` function for each dataset (one at a time)
HybridModelSources(1)
HybridModelSources(2)
HybridModelSources(3)
HybridModelSources(4)
HybridModelSources(5)
HybridModelSources(6)
HybridModelSources(7)
HybridModelSources(8)
