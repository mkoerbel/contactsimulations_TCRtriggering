""" 
Automatically start simulations and loop over parameters

"""

import os
import numpy as np
import pandas as pd
import random
import concurrent.futures
from itertools import repeat

# IO Parameters

smoldyn_path = "/usr/local/bin/smoldyn" #exe file to call smoldyn
root_path = "/Users/koerbem/Documents/smoldyn/" 
model_file = "model_1_9Phos_kphos.txt"
results_folder = "/Users/koerbem/Documents/smoldyn/Model_1_9Phos_scan-k/"
simnum = 100

# Model Parameters
pMHC_mdens=[0,10]                           # mol/um^2 
k_phos=[0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5]    # s-1; phosphorylation rate
randomint = 42                              # answer to everything else

# Function Definitions

def simulate(smoldyn_path, root_path, model_file, results_folder, simnum, pMHC_dens, k_phos):
    os.chdir(root_path)
    string = smoldyn_path + ' ' + model_file + " --define python=42 --define results_folder="+results_folder+" --define sim_num={:.0f} --define randomint={:.0f} --define pMHC_mdens={:.0f} --define k_phos={:.1f} -tqw".format(simnum, random.randint(1,1e10), pMHC_dens, k_phos)
    print(string)
    os.system(string)
    return 0

# Main

if __name__ == '__main__':
    # Run simulations
    random.seed()
    print("Chosen model: "+root_path+model_file)
    print("Starting Simulation")

    for i_phos in k_phos:
        for i_mhc in pMHC_mdens:
            print('Simulating k_phops {} and pMHC_mdens {}'.format(i_phos, i_mhc))
            with concurrent.futures.ProcessPoolExecutor() as executor:
                results = executor.map(simulate, repeat(smoldyn_path), repeat(root_path), repeat(model_file), repeat(results_folder), range(1, simnum+1), repeat(i_mhc), repeat(i_phos))
                

    print("\r Well Done!")