import pyatoa
config = pyatoa.Config()
mgmt = pyatoa.Manager(config=config, empty=True)
mgmt.plot_map(show=True, show_faults=False, annotate_names=True,
              figsize=(10, 9.4), dpi=100)
