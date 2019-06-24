#!/usr/bin/env python3
"""
Preconfigured maps. Calls on plot_map.py to generate basemaps and creates ontop
e.g. catalog maps, ray density maps. etc.

"""
import matplotlib.pyplot as plt
import pyatoa.utils.visuals.plot_map as plot_map


def ray_density(cat, inv, map_corners=[-42.5007, -36.9488, 172.9998, 179.5077],
                **kwargs):
    """
    Take all events in a catalog and all stations in an inventory and connects
    them by straight lines to show ray coverage
    :param cat:
    :param inv:
    :return:
    """
    figsize = kwargs.get("figsize", (10, 9.4))
    dpi = kwargs.get("dpi", 100)
    show = kwargs.get("show", True)
    save = kwargs.get("save", None)

    f = plt.figure(figsize=figsize, dpi=dpi)
    m = plot_map.initiate_basemap(map_corners=map_corners, scalebar=True,
                                  coastline_zorder=101)

    for event in cat:
        # plot_map.event_beachball(m, event, type="focal_mechanism", width=2.6E4)
        for net in inv:
            for sta in net:
                plot_map.connect_source_receiver(m, event, sta, linestyle="-",
                                                 markercolor="w", linewidth=0.5,
                                                 linecolor="k", zorder=100)





def event_catalog(cat)
    """
    """

    inputs = np.load(fid)
    lats = inputs["lats"]
    lons = inputs["lons"]
    depths = inputs["depths"]
    mags = inputs["mags"]
    dates = inputs["dates"]
    ids = inputs["ids"]
    depth_color = build_colormap(depths)
    import ipdb;ipdb.set_trace()
    # ADDITION TO JANKILY GET MOMENT TENSORS
    import pyatoa
    for i, eid in enumerate(ids):
        if eid in ['2013p613824',
                    '2016p858951',
                    '2015p768477',
                    '2016p858260',
                    '2017p735876',
                    '2017p819775',
                    '2017p795065',
                    '2018p130600',
                    '2017p708412',
                    '2017p861155']:
            try:
                config = pyatoa.Config(event_id=eid)
                mgmt = pyatoa.Manager(config)
                width = mgmt.event.preferred_magnitude().mag**2*6E2
                color = depth_color.to_rgba(depths[i])
                event_beachball(m, mgmt.event, width=width, facecolor=color)
            except Exception as e:
                print(e)
                continue


    # ADDITION STOP

    # mags = [m**2.6 for m in mags]
    #
    # depth_color = build_colormap(depths)
    # x, y = m(lons, lats)
    # m.scatter(x, y, marker='o', edgecolor='k', s=mags, zorder=999,
    #           linewidth=0.5, c=depth_color.to_rgba(depths),
    #           )
    # for M, D in zip([4.5, 6.0], [0, 80]):
    #     m.scatter(0, 0, marker='o', edgecolor='k', s=M**2.6, zorder=1,
    #               linewidth=0.5, label="M{}, z={}km".format(M, D),
    #               c=depth_color.to_rgba(D)
    #               )

