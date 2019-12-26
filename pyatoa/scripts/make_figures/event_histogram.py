"""
Create a histogram of event depths for a given Catalog read in from an
event xml file. Colors coordinated with Tape et al. 2010
"""
from obspy import read_events
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# Set some plotting parameters
mpl.rcParams['font.size'] = 17
mpl.rcParams['font.weight'] = "normal"
mpl.rcParams['lines.linewidth'] = 1.25
mpl.rcParams['axes.linewidth'] = 3 

# Grab depth information from the Catalog
cat = read_events('/path/to/catalog_xml')
depths = []
for event in cat:
    depths.append(event.preferred_origin().depth * 1E-3)
depths = np.array(depths)

# Colors and binsizes for depth following Tape et al. 2010
colors = ["red", "gold", "yellowgreen", "royalblue", "royalblue", "deepskyblue", 
          "deepskyblue", "blueviolet", "blueviolet", "fuchsia"]
bin_step = 5
binsize = len(np.arange(0, bin_step, 2.5))

# Create individual histograms w/ respective colors for each given depth bin
f, ax = plt.subplots(1, figsize=(5.25, 6))
for i, bin_start in enumerate(range(0, 60, bin_step)):
    bin_end = bin_start + bin_step
    binned = np.where((depths > bin_start) & (depths < bin_end))
    depths_binned = depths[binned]
    
    if i > len(colors)-1:
        color = colors[-1]
    else:
        color = colors[i]

    plt.hist(x=depths_binned, bins=binsize, color=color, edgecolor='k',
             linewidth=1., histtype='bar', orientation='horizontal',
             range=(bin_start, bin_end), rwidth=1)

# Finish up with some labels
plt.ylabel("Depth (km)")
plt.xlabel("Number")
plt.tick_params(which='both', direction='in', top=True, right=True)
plt.tick_params(which='major', length=8, width=1.5)
plt.gca().invert_yaxis()
plt.show()

