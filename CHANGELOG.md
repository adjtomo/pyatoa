# Current (moving to v0.2.0)

- Replaced Basemap mapping dependency with Cartopy due to Basemap EOL
- Refactored Gatherer class to simplify code structure. API remains the same.
- New Executive class takes care of standalone parallel processing 
- Pyaflowa per-station parallel processing using concurrent.futures
