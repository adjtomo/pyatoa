"""
Generate a large summary figure with multiple histograms using the Inspector
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
from pyatoa import Inspector


fid = "nz1d"
model = "m00"
step = "s00"
comps = ["R", "T", "Z"]
comp_choice = "cc_shift_in_seconds"
cb = 1.75

insp = Inspector(fid)
windows = insp.windows.copy()
choices = ["cc_shift_in_seconds", "dlnA", "max_cc_value", 
           comp_choice, comp_choice, comp_choice,
           "relative_starttime", comp_choice, comp_choice
           ]
binsizes = [cb, 0.25, 0.15, cb, cb, cb, 40, cb, cb]
ymaxes = [1400, 1400, 1400, 550, 550, 550, 3200, 1000, 1000]
xlabels = [None, None, None, 
           "R Time Shift (s)", "T Time Shift (s)", "Z Time Shift (s)",
           None, "<100s Time Shift (s)", ">100s Time Shift (s)"]
          

f = plt.figure()
gs = mpl.gridspec.GridSpec(3, 3, wspace=0.75, hspace=0.75)
i = 0
for row in range(0, gs.get_geometry()[0]):
    for col in range(0, gs.get_geometry()[1]):
        ax = plt.subplot(gs[row, col])
        if row == 1:
           insp.windows = windows.loc[windows["component"] == comps[col]]
        elif row == 2:
            if col == 0:
                insp.windows = windows
            elif col == 1:
                insp.windows = windows.loc[windows["relative_starttime"] < 100]
            elif col == 2:
                insp.windows = windows.loc[windows["relative_starttime"] > 100]
        insp.hist(model, step, choice=choices[i], binsize=binsizes[i],
                  f=f, ax=ax, show=False, figsize=(4.5,6), fontsize=8,
                  ymax=ymaxes[i], xlabel=xlabels[i])
        i += 1
plt.show()

