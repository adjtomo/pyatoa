:py:mod:`pyatoa.utils.images`
=============================

.. py:module:: pyatoa.utils.images

.. autoapi-nested-parse::

   Utility functions for manipulating image files such as .png and .pdfs.
   Makes use use of the Python Pillow package if working with .png files, and
   the PyPDF2 package if manipulating pdf files. Internal imports for all functions
   to remove Pyatoa-wide dependencies on these packages for short functions.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   pyatoa.utils.images.merge_pdfs
   pyatoa.utils.images.imgs_to_pdf
   pyatoa.utils.images.tile_imgs
   pyatoa.utils.images.tif_to_array



.. py:function:: merge_pdfs(fids, fid_out)

   Merge a list of pdfs into a single output pdf using the PyPDF2 package.
   Any desired order to the pdfs should be set in the list of input fids.

   :type fids: list
   :param fids: list of paths to .pdf files
   :type fid_out: str
   :param fid_out: path and name of the resulting output .pdf file


.. py:function:: imgs_to_pdf(fids, fid_out)

   Combine a list of .png files into a single PDF document

   :type fids: list
   :param fids: list of file ids with full pathnames to be combined
   :type fid_out: str
   :param fid_out: the name of the file to be saved with full pathname


.. py:function:: tile_imgs(fids, fid_out)

   Combine a list of images into a single, horizontally tiled image.

   :type fids: list
   :param fids: list of file ids with full pathnames to be tiled
   :type fid_out: str
   :param fid_out: the name of the file to be saved with full pathname


.. py:function:: tif_to_array(fid)

   Convert GeoTiff images (e.g., ETOPO1 topography) to a numpy array for
   conversion and processing

   :type fid: str
   :param fid: .tif(f) file
   :rtype: np.array
   :return: array of data contained within the tiff file


