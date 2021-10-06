"""
Given ObsPy Catalog and Inventory objects, create a raypath density plot to 
provide a more deatiled illustration of raypath gradients, which may be 
interpreted alongside tomographic inversion results as a preliminary resolution
test.

The idea behind this script is that it will partition a single raypath line into 
discrete points (adjustable), divide the entire 2D domain into a grid of 
adjustable squares, and then count the number of descrete point within each
grid square. A contour plot is then generated using these 'pixels'
"""
from obspy import read_events, read_stations
