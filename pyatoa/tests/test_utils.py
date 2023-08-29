"""
Test all the various utility functions except for the processing functions,
which get their own test suite.
"""
import os
import shutil
import pytest
import numpy as np
from glob import glob
from obspy import UTCDateTime, read_inventory
from obspy.core.util.testing import streams_almost_equal
from pyasdf import ASDFDataSet
from pyatoa.utils import (adjoint, calculate, form, images, srcrcv, window, 
                          write)


@pytest.fixture
def station_fid():
    """Re-used test data pointing to STATIONS file"""
    return "./test_data/test_STATIONS_NZ_BFZ"


@pytest.fixture
def sem_fid(tmpdir):
    """
    Re-used test data pointing to a two-column ascii file representing
    SPECFEM3D seismograms
    """
    src = "./test_data/test_adjoint_source_NZ_BFZ_N_2018p130600.adj"
    dst = os.path.join(tmpdir, "NZ.BFZ.BXZ.semd")
    shutil.copy(src, dst)
    return dst


@pytest.fixture
def ds():
    """Test ASDFDataSet"""
    fid = "./test_data/test_ASDFDataSet.h5"
    return ASDFDataSet(fid)


# ============================= TEST ADJOINT UTILS =============================
def test_traveltime_adjoint_source():
    """
    !!! TO DO
    """
    pass

# ============================= TEST CALCULATE UTILS ===========================
def test_calc_utils():
    """
    Wrapping all the calculate utilities into a single function because they're
    all pretty basic, but it's good to know they function the same each time
    """
    # Test abs_max
    assert(calculate.abs_max(np.array([-100, 10])) == -100)

    # Test myround
    value = 34.501
    base = 5
    assert(calculate.myround(value, base, "near") == 35)
    assert(calculate.myround(value, base, "up") == 35)
    assert(calculate.myround(value, base, "down") == 30)

    # Test overlapping_adys
    otime = UTCDateTime("2000-01-01T00:00:00")  # midnight
    daylist = calculate.overlapping_days(otime, start_pad=20, end_pad=200)
    assert(daylist == [365, 1])

    daylist = calculate.overlapping_days(otime, start_pad=0, end_pad=0)
    assert(daylist == [1])

    # Test normalize_a_to_b
    arr = np.arange(0, 101, 10)
    norm_arr = calculate.normalize_a_to_b(arr, a=0, b=1)
    assert((norm_arr == (arr / arr.max())).all())

    # Test amplitude_anomaly
    # !!! This is only checking that the function returns, not output values
    dx = 0.5
    x = np.arange(-2 * np.pi, 2 * np.pi, dx)
    y1 = np.sin(x)
    y2 = np.sin(x + np.pi / 4)  # shift by quarter pi to get some offset
    aa = calculate.amplitude_anomaly(y1, y2, dx)
    assert(aa == pytest.approx(0.0023, .1))

    # Test logarithmic variance reduction - reuse amplitude anomaly arrays
    output = calculate.vrl(y1, y2, y2)
    assert(output == pytest.approx(6.735, .1))

# ============================= TEST FORMAT UTILS ==============================
def test_form_utils():
    """
    Test the string formatting utilities which ensures that IO is standardized
    throughout the package. Again, place them all in the same function.
    """
    # Test format_iter
    for check in ["i09", 9, "9"]:
        assert(form.format_iter(check) == "i09")

    assert(form.format_iter(9.) is None)

    # Test format step
    for check in ["s09", 9, "9"]:
        assert(form.format_step(check) == "s09")

    assert(form.format_step(9.) is None)

    # Test format event name with various test cases
    eid = "abc123"  # test event id
    test_cases = [
        "smi:nz.org.geonet/{eid}",  # GEONET
        "smi:service.iris.edu/fdsnws/event/1/query?eventid={eid}",  # IRIS
        "smi:local/ndk/{eid}/event",  # SPUD
        "smi:local/cmtsolution/{eid}/event",  # GCMT
        "quakeml:earthquake.usgs.gov/fdsnws/event/1/query?eventid={eid}"
                                                    "&format=quakeml",  # USGS
        "quakeml:us.anss.org/event/{eid}",  # USGS COMCAT,
                  ]
    for test_case in test_cases:
        test_case = test_case.format(eid=eid)
        assert(form.format_event_name(test_case) == eid)



# ============================= TEST READ/ WRITE UTILS =========================


def test_write_utils(tmpdir, station_fid, sem_fid, ds):
    """
    Test write utilities
    """
    # Test write_inv_seed with edited dir structure
    # Mess with the default template naming schema and see if it returns goodly
    inv = read_inventory()
    write.write_inv_seed(inv, path=tmpdir)

    check_fids = ["FUR.GR/RESP.GR.FUR..BHE", "RJOB.BW/RESP.BW.RJOB..EHZ",
                  "WET.GR/RESP.GR.WET..LHN"]
    for check_fid in check_fids:
        assert(os.path.exists(os.path.join(tmpdir, check_fid)))

    # Test write_misfit
    event_id = form.format_event_name(ds)
    write.write_misfit(ds, iteration=1, step_count=0, path=tmpdir)
    misfit = np.loadtxt(os.path.join(tmpdir, event_id), dtype=float)
    assert(misfit == pytest.approx(0.2512, 3))

    # Test write_stations_adjoint
    # This will only write out NZ.BFZ because STATION file only contains BFZ
    write.write_stations_adjoint(ds=ds, iteration="i01",
                                 specfem_station_file=station_fid,
                                 pathout=tmpdir
                                 )
    sta_adj = np.loadtxt(os.path.join(tmpdir, "STATIONS_ADJOINT"), dtype=str)
    assert(sta_adj[0] == "BFZ")

    # Test write_adj_src_to_ascii
    sta = "NZ_BFZ_BXN"
    fid = sta.replace("_", ".") + ".adj"

    write.write_adj_src_to_ascii(ds, iteration="i01", pathout=tmpdir)
    max_val = ds.auxiliary_data.AdjointSources.i01.s00[sta].data[()].max()

    path_check = os.path.join(tmpdir, fid)
    check_max = np.loadtxt(path_check).max()
    assert(max_val == pytest.approx(check_max, .1))

# ============================= TEST SRCRCV UTILS ==============================
# srcrcv functions are all pretty short, likely don't need a test

# ============================= TEST WINDOW UTILS ==============================
# not enough window utils to warrant writing tests

