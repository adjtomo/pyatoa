"""
Simple scatterplot to show misfit convergence based on the misfits.json
output that is created during a Seisflows plugin run
"""
import os
import json
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['font.size'] = 14
mpl.rcParams['lines.linewidth'] = 2.
mpl.rcParams['lines.markersize'] = 6
mpl.rcParams['axes.linewidth'] = 3

def parse_json(fid):
    """
    Parse the misfits json file and 
    """
    misfits = json.load(open(fid))

    misfit_array, window_array, model = [], [], []
    for i, key in enumerate(misfits):
        try:
            misfit = misfits[key]['s00']['misfit']
            windows = misfits[key]['s00']['windows']
            # if s00 is available but is the same as the previous misfit,
            # force the exception to get s01
            if misfit_array and misfit == misfit_array[-1]:
                raise KeyError
        # if thrifty inversion, s00 is taken from the previous iteration
        except KeyError:
            misfit = misfits[key]['s01']['misfit']
            windows = misfits[key]['s01']['windows']
        
        misfit_array.append(misfit)
        window_array.append(windows)
        model.append(i)
    
    return model, misfit_array, window_array


if __name__ == "__main__":
    # For comparisons of various misfit criteria
    basepath = ("/scale_wlg_nobackup/filesets/nobackup/nesi00263/bchow/tomo/"
                "seisflows/checkerboard/30event_75e1785/")
    fid_list = [os.path.join(basepath, "cc/pyatoa.io/misfits.json"),
                os.path.join(basepath, "mtm/pyatoa.io/misfits.json"),
                os.path.join(basepath, "both/pyatoa.io/misfits.json")
                ]    
    color_list = ["mediumpurple", "darkorange", "mediumturquoise"]
    label_list = ["traveltime_cc", "multitaper misfit", "both"]

    # set up plots
    f, ax1 = plt.subplots(figsize=(8,6))
    ax2 = ax1.twinx()

    for i, fid in enumerate(fid_list):
        models, misfits, windows = parse_json(fid)
        ax1.plot(models, misfits, 'o-', c=color_list[i], label=label_list[i])
        ax2.plot(models, windows, 'd-.', c=color_list[i], alpha=0.5)
    
    ax1.legend(loc="center right")
    plt.title('Checkerboard Convergence')
    ax1.set_xticks(models)
    ax1.set_xlabel('Model Number')
    ax1.set_ylabel('Total Misfit (solid)')
    ax2.set_ylabel('Number of Windows (dashed)', rotation=270, labelpad=15.)
    ax1.grid(True, alpha=0.5, linestyle='--', linewidth=1.)
    
    plt.savefig('convergence.png') 
    plt.show()
        
