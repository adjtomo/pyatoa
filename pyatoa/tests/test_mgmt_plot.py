"""
Test plotting capabilities by comparing baseline images for example data
"""
import os
import pytest
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.testing.compare import compare_images
from pyatoa import logger, Manager

# Turn off the logger for tests
logger.propogate = False
logger.setLevel("CRITICAL")

IMGDIR = "./test_data/baseline_images"


@pytest.fixture
def mgmt():
    """
    Fully realized manager object which will be used to make plots
    """
    mgmt = Manager()
    mgmt.load()
    mgmt.flow(remove_response=True, output="DISP")
    return mgmt


def reset_matplotlib():
    """
    Reset matplotlib to a common default. Copied from Pyflex
    """
    # Set all default values.
    mpl.rcdefaults()
    # Force agg backend.
    plt.switch_backend('agg')
    # These settings must be hardcoded for running the comparision tests and
    # are not necessarily the default values.
    mpl.rcParams['font.family'] = 'Bitstream Vera Sans'
    mpl.rcParams['text.hinting'] = 'none'  # used to be False, invalid

    # Not available for all matplotlib versions.
    try:
        mpl.rcParams['text.hinting_factor'] = 8
    except KeyError:
        pass

    # Force classic style.
    mpl.style.use("classic")

    import locale
    locale.setlocale(locale.LC_ALL, str('en_US.UTF-8'))


def images_are_identical(created_img):
    """
    Modified from Pyflex which was partially copied from ObsPy.
    Used to check images for equality.
    """
    baseline_img = os.path.join(IMGDIR, os.path.basename(created_img))

    assert os.path.exists(created_img)
    assert os.path.exists(baseline_img)

    # Use a reasonably high tolerance to get around difference with different
    # freetype and possibly agg versions.
    # 25 is fairly high but this should work on matplotlib 1 and 2.
    result = compare_images(created_img, baseline_img, 25, in_decorator=True)
    if result is not None:
        print(result)
    assert result is None


def test_waveform_plot(tmpdir, mgmt):
    """
    Test that plotting waveforms, windows, adjsrc etc. by themselves works
    """
    # reset_matplotlib()
    save_fid = os.path.join(tmpdir, "default_manager_plot_wav.png")
    mgmt.plot(choice="wav", show=False, save=save_fid)
    images_are_identical(save_fid)


def test_map_plot(tmpdir, mgmt):
    """
    Test that map plotting with moment tensor works
    """
    save_fid = os.path.join(tmpdir, "default_manager_plot_map.png")
    mgmt.plot(choice="map", show=False, save=save_fid)
    images_are_identical(save_fid)


def test_combined_plot(tmpdir, mgmt):
    """
    Test that a combined gridspec waveform + map figure works
    """
    save_fid = os.path.join(tmpdir, "default_manager_plot_both.png")
    mgmt.plot(choice="both", show=False, save=save_fid)
    images_are_identical(save_fid)
