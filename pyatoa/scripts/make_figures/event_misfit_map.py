"""
Create an event misfit map from the auxiliary data in a Pyasdf
Dataset.
"""
import os
import glob
import pyasdf
from pyatoa.core.config import Config
from pyatoa.utils.visuals.mapping import event_misfit_map

cfg = Config()

path_to_ds = os.getcwd()

for dataset in glob.glob(os.path.join(path_to_ds, '*.h5')):
    with pyasdf.ASDFDataSet(dataset) as ds:
        event_id = ds.events[0].resource_id.id.split('/')[-1]
        for model in ds.auxiliary_data.Statistics.list():
            # Normalize by the first model
            if model == "m00":
                fidout = "{eid}_{mdl}_misfitmap.png".format(eid=event_id,
                                                            mdl=model)
                event_misfit_map(map_corners=cfg.map_corners,
                                 ds=ds, model=model,
                                 annotate_station_info='simple',
                                 contour_overlay=True, filled_contours=True,
                                 show=False, save=fidout)

