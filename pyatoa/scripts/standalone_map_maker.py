import pyatoa
config = pyatoa.Config()
mgmt = pyatoa.Manager(config=config, empty=True)
mgmt.plot_map(show=True, show_faults=False, color_by_network=True,
              figsize=(10, 9.4), dpi=100,
              map_corners=[-42.65, -36.8, 172.95, 179.5])
# [-42.5007, -36.9488, 172.9998, 179.5077]

