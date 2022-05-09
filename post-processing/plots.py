""" 
PLOTS
======================

Reads .h5 files and plots the results 
    
(c) Matteo M. 2022

"""

import os, fnmatch
import numpy as np
import matplotlib.pyplot as plt
import h5py


##########   SETTINGS  #############
results_dir = '../results/'    # Results files location
plot_profiles = True           # Plot profiles along the column
skip_profiles = 0              # Skip n profiles
ds = 1                         # Plot every ds samples
units = 2                      # 1: ppb, 2:ug/m3
####################################


# List of files
def listfiles(directory):
    listOfFiles = os.listdir(directory)
    pattern = "*.h5"
    l = []
    for entry in listOfFiles:
        if fnmatch.fnmatch(entry, pattern):
            l.append(entry)
    if len(l) < 1:
        print('Error. No .h5 file found.')
    return l



# Print the list of files and select one from the list
l_files = listfiles(results_dir) 
num_files = len(l_files)
for i in range(num_files):
    print(f'[{i}] {l_files[i]}')
print(f'Select a file (default: last file) [{num_files-1}]:', end = '')
num = input()
if len(num) > 0:
    sel = int(num)
else:
    sel = num_files-1

filename = results_dir + l_files[sel]



########### READ .H5 FILES ############
hf = h5py.File(filename, 'r')

model = {}
model['t'] = hf['/results/output/t'][:]
model['c'] = hf['/results/output/c'][:]


use_ts = hf['parameters/use_ts']
if use_ts[()]:
    inlet = {}
    inlet['t'] = hf['data/inlet/t'][:]
    inlet['c'] = hf['data/inlet/c'][:]

    outlet = {}
    outlet['t'] = hf['data/outlet/t'][:]
    outlet['c'] = hf['data/outlet/c'][:]

export_profiles = hf['parameters/export_profiles']
if export_profiles[()]:
    profiles = {}
    profiles['t'] = hf['/results/profiles/t'][:]
    profiles['z'] = hf['/results/profiles/z'][:]
    profiles['o'] = hf['/results/profiles/o']
    profiles['c'] = hf['/results/profiles/c']

mixing_model = hf['parameters/mixing_model']
if mixing_model[()]:
    model['c_mix'] = hf['/results/output/c_mix'][:]
	
# Print simulation filename (.yml)
print(hf['parameters/filename'].asstr()[()])

	
################ PLOTS ################

# Read molar weight
mw = hf['parameters/mw'][()]

# Conversion factor
if units == 2:
    cf = mw/24.45
else:
    cf = 1


# Output concentrations
plt.figure()

if use_ts[()]:
    plt.plot(inlet['t']/3600,inlet['c']*cf,linestyle='-',marker='',fillstyle='none', label='Inlet')
    plt.plot(outlet['t']/3600,outlet['c'],linestyle='-',marker='',fillstyle='none', label='Outlet')

plt.plot(model['t'][::ds]/3600,model['c'][::ds]*cf,linestyle='none',marker='o',fillstyle='none', label='Column outlet')

if mixing_model[()]:
    plt.plot(model['t'][::ds]/3600,model['c_mix'][::ds]*cf,linestyle='none',marker='o',fillstyle='none', label='Inlet air - mixed')

plt.legend()
plt.xlabel('Time (h)')
if units == 2:
    plt.ylabel('Concentration ($\mu$g/m$^3$)')
else:
    plt.ylabel('Concentration (ppb)')


if (export_profiles[()] == 1) and plot_profiles:
    # Profiles
    # Oxygen
    plt.figure()

    for i in range(len(profiles['t'])):
        if i > skip_profiles:
            plt.plot(profiles['z'],profiles['o'][i][:], linestyle='-',marker='',fillstyle='none', label=f"t = {profiles['t'][i]/3600:.1f} h")
    plt.legend()
    plt.xlabel('Distance from inlet (m)')
    plt.ylabel('Oxygen concentration (ppm)')

    # VOC
    plt.figure()
    for i in range(len(profiles['t'])):
        if i > skip_profiles:
            plt.plot(profiles['z'],profiles['c'][i][:]*cf, linestyle='-',marker='',fillstyle='none', label=f"t = {profiles['t'][i]/3600:.1f} h")
    plt.legend()
    plt.xlabel('Distance from inlet (m)')
    if units == 2:
        plt.ylabel('Concentration ($\mu$g/m$^3$)')
    else:
        plt.ylabel('Concentration (ppb)')

plt.show()




