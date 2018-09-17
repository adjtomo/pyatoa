"""
General utility functions used in plotting scripts
"""


def align_yaxis(ax1, ax2):
    """
    adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1

    :type ax?: matplotlib axis
    :param ax?: axes to adjust
    """
    ymin_a1, ymax_a1 = ax1.get_ylim()
    ymin_a2, ymax_a2 = ax2.get_ylim()

    _, y1 = ax1.transData.transform((0, (ymax_a1+ymin_a1)/2))
    _, y2 = ax2.transData.transform((0, (ymax_a2+ymin_a2)/2))
    inv = ax2.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, y1-y2))
    ax2.set_ylim(ymin_a2+dy, ymax_a2+dy)


def pretty_grids(input_ax):
    """
    standard plot skeleton formatting, thick lines and internal tick marks etc.

    :type input_ax: matplotlib axis
    :param input_ax: axis to prettify
    """
    input_ax.set_axisbelow(True)
    input_ax.tick_params(which='both', direction='in', top=True, right=True)
    input_ax.minorticks_on()
    input_ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    for axis_ in ['major', 'minor']:
        input_ax.grid(which=axis_, linestyle=':', linewidth='0.5', color='k',
                      alpha=0.25)


def format_axis(input_ax):
    """
    Sit the tick marks away from the plot edges to prevent overlapping when
    multiple subplots are stacked atop one another, and for general gooood looks
    will check if the plot is two sided (e.g. waveforms) or only positive

    :type input_ax: matplotlib axis
    :param input_ax: axis to prettify
    """
    ymin, ymax = input_ax.get_ylim()
    maxvalue = max([abs(_) for _ in input_ax.get_ylim()])
    percentover = maxvalue * 0.125
    if abs(round(ymin/ymax)) != 0:
        bounds = (-1 * (maxvalue+percentover), (maxvalue+percentover))
    else:  # elif abs(round(ymin/ymax)) == 0:
        bounds = (-0.05,1.05)
    input_ax.set_ylim(bounds)