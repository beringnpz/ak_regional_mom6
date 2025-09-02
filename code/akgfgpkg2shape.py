import os
import geopandas as gpd

akgffol = "/Users/kelly.kearney/Documents/Repos/Other/akgfmaps/inst/extdata/"

fbase = ["afsc_bts_strata", "all_crab_from_akgfmaps_grid"]

for ff in fbase:

    fname = os.path.join(akgffol, f"{ff}.gpkg")
    layers = gpd.list_layers(fname)

    for lay in layers.name:
        tmp = gpd.read_file(fname, layer=lay)
        tmp.to_file(f"../supporting_data/akgfmaps_shapefiles/{ff}_{lay}.shp", driver="ESRI Shapefile")