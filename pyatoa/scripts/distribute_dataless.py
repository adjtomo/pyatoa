"""
Write an ObsPy Inventory object, which can be comprised of various 

Response files written through Obspy come out as a single object, but Pyatoa
will look for response information from individual components and individual
stations. Distrubute this dataless information into the necessary components
"""
import os
import sys
from obspy import read_inventory

def distribute_dataless(inv, path="./"):
    """
    Note: The template here is hardcoded with SEED convention

    :type inv: obspy.core.inventory.Inventory
    :param inv: inventory with response to be delinieated into separate objects
    :type path: str
    :param path: path to save the new response files to
    """
    # This follows SEED convention
    inner_folder = "{STA}.{NET}"
    fid_template = "RESP.{NET}.{STA}.{LOC}.{CHA}"

    full_template = os.path.join(path, inner_folder, fid_template)
    for net in inv:
        for sta in net:
            try:
                # Create the container directory unless it exists
                os.mkdir(os.path.join(path, inner_folder.format(STA=sta.code, 
                                                                NET=net.code))
                )
            except FileExistsError:
                pass

            for cha in sta:
                # Write the individual channel as a STATIONXML file
                inv_temp = inv.select(network=net.code, station=sta.code,
                                      location=cha.location_code,
                                      channel=cha.code)
                inv_temp.write(full_template.format(
                    STA=sta.code, NET=net.code, LOC=cha.location_code,
                    CHA=cha.code),
                    format="STATIONXML"
                    )

if __name__ == "__main__":
    path = "./"
    inv = read_inventory(sys.argv[1])
    distribute_dataless(inv, path)     
