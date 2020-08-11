"""
Utility functions for manipulating image files such as .png and .pdfs, using the
Python Pillow package
"""
from PIL import Image


def imgs_to_pdf(fids, fid_out):
    """
    Combine a list of image files into a single PDF document

    :type fids: list
    :param fids: list of file ids with full pathnames to be combined
    :type fid_out: str
    :param fid_out: the name of the file to be saved with full pathname
    """
    images = []
    for fid in fids:
        # PNGs need to be converted to RGB to get alpha to play nice
        images.append(Image.open(fid).convert("RGB"))

    image_main = images[0]
    images = images[1:]

    image_main.save(fp=fid_out, format="PDF", resolution=100., save_all=True,
                    append_images=images)


def tile_imgs(fids, fid_out):
    """
    Combine a list of images into a single, horizontally tiled image.

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


