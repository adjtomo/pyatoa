"""
Utility functions for manipulating image files such as .png and .pdfs.
Makes use use of the Python Pillow package if working with .png files, and 
the PyPDF2 package if manipulating pdf files. Internal imports for all functions
to remove Pyatoa-wide dependencies on these packages for short functions.
"""


def merge_pdfs(fids, fid_out):
    """
    Merge a list of pdfs into a single output pdf using the PyPDF2 package.
    Any desired order to the pdfs should be set in the list of input fids.

    :type fids: list
    :param fids: list of paths to .pdf files
    :type fid_out: str
    :param fid_out: path and name of the resulting output .pdf file
    """
    if not fids:
        return

    from PyPDF2 import PdfFileMerger
    
    merger = PdfFileMerger()
    for fid in fids:
        merger.append(fid)

    merger.write(fid_out)
    merger.close()


def imgs_to_pdf(fids, fid_out):
    """
    Combine a list of .png files into a single PDF document

    :type fids: list
    :param fids: list of file ids with full pathnames to be combined
    :type fid_out: str
    :param fid_out: the name of the file to be saved with full pathname
    """
    from PIL import Image

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
    from PIL import Image

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


