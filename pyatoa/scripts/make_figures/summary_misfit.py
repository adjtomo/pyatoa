"""
Generate a summary figure with multiple histograms using the Inspector class
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pyatoa import Inspector


def gen_fig(insp, choice="3x3"):
    """
    A larger summary misfit that includes separation by component and starttime
    """
    # The choices for the 3 rows and 3 columns of the figure
    if choice == "3x3":
        time_split = 100
        choices = [
          ["cc_shift_in_seconds", "dlnA", "max_cc_value"],
          ["cc_shift_in_seconds", "cc_shift_in_seconds", "cc_shift_in_seconds"],
          ["relative_starttime", "cc_shift_in_seconds", "cc_shift_in_seconds"]
          ]
        binsizes = [[2, None, 0.1], 
                    [2, 2, ,2]
    elif choice == "2x2":
        choices = [["cc_shift_in_seconds", "dlnA"],
                   ["max_cc_value", "relative_starttime"]]
    else:
        raise NotImplementedError


    # Instantiate the 3x3 plot
    f = plt.figure()
    gs = mpl.gridspec.GridSpec(len(choices), len(choices[0]), 
                               wspace=0.75, hspace=0.75)

    windows = insp.windows.copy()
    comps = windows.component.unique()

    for row in range(0, gs.get_geometry()[0]):
        for col in range(0, gs.get_geometry()[1]):
            ax = plt.subplot(gs[row, col])

            if choice == "3x3":
                # Plot middle row by component if 3x3 figure
                if row == 1:
                   insp.windows = windows.loc[
                           windows["component"] == comps[col]]
                # Plot last row by splitting on a certain time value
                elif row == 2:
                    if col == 0:
                        insp.windows = windows
                    elif col == 1:
                        insp.windows = windows.loc[
                                windows["relative_starttime"] < time_split]
                    elif col == 2:
                        insp.windows = windows.loc[
                                windows["relative_starttime"] >= time_split]
                else:
                    insp.windows = windows

            insp.hist(choice=choices[row][col], 
                      # binsize=binsizes[i], xlabel=xlabels[i], ymax=ymaxes[i], 
                      f=f, ax=ax,  show=False, figsize=(5.0, 6.5), fontsize=8)
    plt.show()



if __name__ == "__main__":
    fid = "sixty"
    insp = Inspector(fid)
    gen_fig(insp)


    # comps = ["N", "E", "Z"]
    # comp_choice = "cc_shift_in_seconds"
    # cb = 1.75
    # 
    # choices = ["cc_shift_in_seconds", "dlnA", "max_cc_value", 
    #            comp_choice, comp_choice, comp_choice,
    #            "relative_starttime", comp_choice, comp_choice
    #            ]
    # binsizes = [cb, 0.25, 0.15, cb, cb, cb, 40, cb, cb]
    # ymaxes = [1400, 1400, 1400, 550, 550, 550, 3200, 1000, 1000]
    # xlabels = [None, None, None, 
    #            "R Time Shift (s)", "T Time Shift (s)", "Z Time Shift (s)",
    #            None, "<100s Time Shift (s)", ">100s Time Shift (s)"]
    #           
    # 
    # plt.show()

