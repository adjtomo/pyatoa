"""
Each window represents a pixel in a figure, supposed to be used to get quick
estimations of comparisons between evaluations but is quite hard to make sense
of. Left here incase I want to revisit the idea
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from pyatoa import Inspector


iteration = "i01"
step_count = "s00"
key = "cc_shift_in_seconds"
insp = Inspector('sixty')
events = insp.events[:10]
stations = insp.stations[:10]

max_val = 0
data = []
for event in events:
    event_data = []
    for station in stations:
        df = insp.isolate(iteration, step_count, event=event, station=station,
                          keys=key)
        val = float(df.abs().max())
        event_data.append(val)
        if val > max_val:
            max_val = val

    data.append(event_data)

# Set any NaN values (unavailable windows) to Red
current_cmap = mpl.cm.get_cmap("Greys")
current_cmap.set_bad(color="red")

# Plot attributes
plt.imshow(np.array(data), interpolation="nearest", cmap=current_cmap)

plt.title(f"Maximum {key}: {max_val}")
ax = plt.gca()
# plt.xticks(np.arange(0.5, len(stations), 1.), rotation=90)
# plt.yticks(np.arange(0, len(events), 1), rotation=0)
ax.set_xticklabels(stations)
ax.set_yticklabels(events)

for tick in ax.xaxis.get_major_ticks():
    tick.label.set_horizontalalignment('right')

for tick in ax.yaxis.get_major_ticks():
    tick.label.set_horizontalalignment('right')

plt.show()

