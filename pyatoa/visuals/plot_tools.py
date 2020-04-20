"""
General utility functions used in plotting scripts
"""
from PIL import Image


def align_yaxis(ax1, ax2):
    """
    adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1

    :type ax1: matplotlib axis
    :param ax1: axes to adjust
    :type ax2: matplotlib axis
    :param ax2: axes to adjust
    """
    ymin_a1, ymax_a1 = ax1.get_ylim()
    ymin_a2, ymax_a2 = ax2.get_ylim()

    _, y1 = ax1.transData.transform((0, (ymax_a1+ymin_a1)/2))
    _, y2 = ax2.transData.transform((0, (ymax_a2+ymin_a2)/2))
    inv = ax2.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, y1-y2))
    ax2.set_ylim(ymin_a2+dy, ymax_a2+dy)


def pretty_grids(input_ax, twax=False, grid=False):
    """
    standard plot skeleton formatting, thick lines and internal tick marks etc.

    :type input_ax: matplotlib axis
    :param input_ax: axis to prettify
    :type twax: bool
    :param twax: If twax (twin axis), do not set grids
    """
    input_ax.set_axisbelow(True)
    input_ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    input_ax.tick_params(which='both', direction='in', top=True, right=True)

    # Set the grids 'on' only if main axis
    if not twax:
        input_ax.minorticks_on()
        if grid:
            for axis_ in ['major', 'minor']:
                input_ax.grid(which=axis_, linestyle=':', linewidth='0.5',
                              color='k', alpha=0.25)


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
        bounds = (-0.05, 1.05)
    input_ax.set_ylim(bounds)


def tile_imgs(fids, fid_out):
    """
    A function to combine two figures horizontally using the Pillow library.
    Tiles them in the order they are specified in the input list.
    So far this is only tested with .png files

    :type fids: list
    :param fids: list of file ids with full pathnames to be tiled
    :type fid_out: str
    :param fid_out: the name of the file to be saved with full pathname
    """
    # .png files require conversion to properly get the alpha layer
    images = []
    for fid in fids:
        images.append(Image.open(fid).convert("RGBA"))

    widths, heights = zip(*(i.size for i in images))
    total_width = sum(widths)
    max_height = max(heights)

    # Create the new image that will be returned
    im_out = Image.new(mode="RGBA", size=(total_width, max_height))
    x_offset = 0
    for im in images:
        im_out.paste(im=im, box=(x_offset, 0))
        x_offset += im.size[0]

    im_out.save(fid_out)


def imgs_to_pdf(fids, fid_out):
    """
    A function to combine a list of image files into a single PDF document

    :type fids: list
    :param fids: list of file ids with full pathnames to be combined
    :type fid_out: str
    :param fid_out: the name of the file to be saved with full pathname
    """
    images = []
    for fid in fids:
        images.append(Image.open(fid).convert("RGB"))

    image_main = images[0]
    images = images[1:]

    image_main.save(fp=fid_out, format="PDF", resolution=100., save_all=True,
                    append_images=images)




