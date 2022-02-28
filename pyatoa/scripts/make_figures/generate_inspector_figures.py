"""
Make all available Inspector plots and save for use in the docs page
"""
import os
from pyatoa import Inspector

# User-defined parameters
tag = "forest"
path_out = "./figures"
insp = Inspector(tag="forest")
insp.read(path="../../tests/test_data/docs_data/inspector_doc/")
show = False

if not os.path.exists(path_out):
    os.mkdir(path_out)

# Let the plotting begin
insp.map(show=show, save=f"{tag}_map.png")

insp.event_depths(xaxis="longitude", show=show, save=f"{tag}_event_depths.png")

insp.raypaths(iteration="i01", step_count="s00", show=show, 
              save=f"{tag}_raypaths.png")

insp.raypath_density(iteration="i01", step_count="s00", show=show,
                     save=f"{tag}_raypath_density.png")

insp.event_hist(choice="magnitude", show=show, save=f"{tag}_event_hist.png")
insp.travel_times(t_offset=-20, constants=[2, 4, 6, 8, 10], show=show,
                  save=f"{tag}_travel_times.png")

insp.plot_windows(iteration="i01", step_count="s00", show=show, 
                  save=f"{tag}_windows.png")

insp.convergence(windows="nwin", show=show, save=f"{tag}_convergence.png")

insp.hist(choice="cc_shift_in_seconds", show=show, save=f"{tag}_cc_hist.png")

insp.scatter(x="relative_starttime", y="max_cc_value", show=show,
             save=f"{tag}_scatter.png")

insp.measurement_hist(iteration="i01", step_count="s00", choice="station",
                      show=show, save=f"{tag}_station_hist.png")

insp.measurement_hist(iteration="i01", step_count="s00", choice="event",
                      show=show, save=f"{tag}_event_hist.png")

insp.station_event_misfit_map(station="BFZ", iteration="i01", step_count="s00",
                              choice="misfit", show=show,
                              save=f"{tag}_BFZ_misfit_map.png")

insp.event_station_misfit_map(event="2018p130600", iteration="i01", 
                              step_count="s00", show=show, choice="misfit",
                              save=f"{tag}_2018p130600_misfit_map.png")

insp.event_misfit_map(choice="misfit", show=show, 
                      save=f"{tag}_event_misfit.png")






