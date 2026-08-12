# cefi.py

import xarray as xr
import pandas as pd
from xoak import SklearnGeoBallTreeAdapter
from arraylake import Client
import time
import pathlib as p
from urllib.parse import urljoin, urlparse
import math
import re
import fnmatch

class cefiportalopts:
    def __init__(self, 
        portalpath="http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/",
        region="nep",
        subdomain="full",
        exptype="hcast",
        freq="daily",
        grid="raw",
        release="latest",
        yyyymmdd="19930101"):

        self._regions     = ""
        self._regionl     = ""
        self._exptypes    = ""
        self._exptypel    = ""
        self._freqs       = ""
        self._freql       = ""
        self._subdomains  = ""
        self._subdomainl  = ""
        
        self._portalpath  = portalpath
        self._release     = release
        self._yyyymmdd    = yyyymmdd
        self._grid        = grid

        self._region      = region
        self._subdomain   = subdomain
        self._exptype     = exptype
        self._freq        = freq
        
        self.region      = self._region
        self.freq        = self._freq
        self.subdomain   = self._subdomain
        self.exptype     = self._exptype
        self.portalpath  = self._portalpath
        self.release     = self._release
        self.yyyymmdd    = self._yyyymmdd
        self.grid        = self._grid

    def __tbllists(self, prop):
        if prop == "region":
            return ["nwa", "northwest_atlantic",
                    "nep", "northeast_pacific",
                    "arc", "arctic",
                    "pci", "pacific_islands",
                    "glk", "great_lakes"]
        elif prop == "exptype":
            return ["hcast",         "hindcast",
                    "ss_fcast",      "seasonal_forecast",
                    "ss_fcast_init", "seasonal_forecast_initialization",
                    "ss_refcast",    "seasonal_reforecast",
                    "dc_fcast_init", "decadal_forecast_initialization",
                    "dc_fcast",      "decadal_forecast",
                    "ltm_proj",      "long_term_projection"]
        elif prop == "freq":
            return ["day", "daily",
                    "mon", "monthly",
                    "ann", "yearly",
                    "sta", "static"]
        elif prop == "subdomain":
              return ["full",     "full_domain",
                      "aksvyreg", "ak_surveyregion_avg"]
        
    @property
    def portalpath(self):
        return self._portalpath
    
    @portalpath.setter
    def portalpath(self,value):
        self._portalpath = value

    @property
    def release(self):
        return self._release
    
    @release.setter
    def release(self,value):
        self._release = value

    @property
    def yyyymmdd(self):
        return self._yyyymmdd
    
    @yyyymmdd.setter
    def yyyymmdd(self,value):
        self._yyyymmdd = value

    @property
    def grid(self):
        return self._grid
    
    @grid.setter
    def grid(self,value):
        self._grid = value

    @property
    def region(self):
        return self._region
    
    @property
    def regions(self):
        return self._regions
    
    @property
    def regionl(self):
        return self._regionl
    
    @region.setter
    def region(self,value):
        tbl = self.__tbllists("region")
        if value in tbl:
            idx = math.floor(tbl.index(value)/2)*2
            self._region = value
            self._regions = tbl[idx]
            self._regionl = tbl[idx+1]
        else:
            print("invalid value for region")

    @property
    def exptype(self):
        return self._exptype
    
    @property
    def exptypes(self):
        return self._exptypes
    
    @property
    def exptypel(self):
        return self._exptypel
    
    @exptype.setter
    def exptype(self,value):
        tbl = self.__tbllists("exptype")
        if value in tbl:
            idx = math.floor(tbl.index(value)/2)*2
            self._exptype = value
            self._exptypes = tbl[idx]
            self._exptypel = tbl[idx+1]
        else:
            print("invalid value for exptype")

    @property
    def freq(self):
        return self._freq
    
    @property
    def freqs(self):
        return self._freqs
    
    @property
    def freql(self):
        return self._freql
    
    @freq.setter
    def freq(self,value):
        tbl = self.__tbllists("freq")
        if value in tbl:
            idx = math.floor(tbl.index(value)/2)*2
            self._freq = value
            self._freqs = tbl[idx]
            self._freql = tbl[idx+1]
        else:
            print("invalid value for freq")

    @property
    def subdomain(self):
        return self._subdomain
    
    @property
    def subdomains(self):
        return self._subdomains
    
    @property
    def subdomainl(self):
        return self._subdomainl
    
    @subdomain.setter
    def subdomain(self,value):
        if isinstance(value, str):
            tbl = self.__tbllists("subdomain")
            if value in tbl:
                idx = math.floor(tbl.index(value)/2)*2
                self._subdomain = value
                self._subdomains = tbl[idx]
                self._subdomainl = tbl[idx+1]
            elif re.fullmatch("iq[0-9]*-[0-9]*jq[0-9]*-[0-9]*", value):
                self._subdomain = value
                self._subdomains = value
                self._subdomainl = value
            else:
                print("invalid value for subdomain")
        else: # TODO: validate 4-element integer-value numeric list
            self._subdomains = f"iq{value[0]}-{value[1]}jq{value[2]}-{value[3]}"
            self._subdomainl = self._subdomains
            self._subdomain = value

    def printprops(self):
        print(f"region:    {self.regions} / {self.regionl}")
        print(f"subdomain: {self.subdomains} / {self.subdomainl}")
        print(f"exp type:  {self.exptypes} / {self.exptypel}")
        print(f"frequency: {self.freqs} / {self.freql}")
        print(f"portal path: {self.portalpath}")
        print(f"release: {self.release}")
        print(f"grid: {self.grid}")

    def cefifolder(self):
        tmp = urlparse(self.portalpath)
        if tmp.scheme in ['http', 'https']:
            return urljoin(self.portalpath, "/".join([self.regionl, self.subdomainl, self.exptypel, self.freql, self.grid, self.release]))
        else:
            return p.Path(self.portalpath) / self.regionl / self.subdomainl / self.exptypel / self.freql / self.grid / self.release

    def cefifilename(self, varname, datestr, flag=False):
        if flag:
            return f"{varname}.{self.regions}.{self.subdomains}.{self.exptypes}.{self.freql}.{self.release}.{datestr}.nc" # temporary error(?) in local scheme, no grid in filename
        else:
            return f"{varname}.{self.regions}.{self.subdomains}.{self.exptypes}.{self.freql}.{self.grid}.{self.release}.{datestr}.nc"
    
    def cefifilelist(self, varname, datestr):
        tmp = urlparse(self.portalpath)
        if tmp.scheme in ['http', 'https']: # parse from OpenDAP catalog

            tmp = self.cefifolder()
            catalog_url = tmp.replace('dodsC', 'catalog')
            
            tbl = pd.read_html("/".join([catalog_url, "catalog.html"]))
            tbl[0].drop(0,inplace=True)

            fglob = self.cefifilename(varname, datestr)
            fmatch = fnmatch.filter(tbl[0].Dataset.to_list(), fglob)

            return ["/".join([self.cefifolder(),x]) for x in fmatch]
        else: # parse from local file listing
            fol = self.cefifolder().expanduser().resolve()
            fglob = self.cefifilename(varname, datestr, True)

            return sorted(list(fol.glob(fglob)))
        
def resamplemodel(fsurvey, fout, model, release, vname, location, 
                  vnamein="", label="", zname="", zsel=0):

    if not label:
        label = f"{model}_{release}_{vname}_{location}"
    
    # Parse model, release, location, and variable name to open model output 
    # as xarray dataset

    if model == "romsb10k":

        # Location of postprocessed ROMS output and grid file
    
        if location == "TDS":
            simbase = "https://data.pmel.noaa.gov/aclim/thredds/dodsC"
            grdpath = "/".join([simbase, "ancillary", "Bering10K_extended_grid.nc"])
        elif location == "klone":
            simbase = p.Path("~/Documents/hpcmount/klone_cicoes/").expanduser().resolve()
            grdpath = simbase / "Bering10K_extended_grid.nc"

        # 5-year blocks used for file names (these vary by release)

        if   release == "K20_CORECFS":
            yr = range(1970,2024,5)
        elif release ==  "K20nobio_CORE":
            yr = range(1970,2004,5)
        elif release ==  "K20nobio_CFS":
            yr = range(1980,2019,5)
        elif release == "K20nobio_CORECFS_daily":
            yr = range(1980,2024,5)
        elif release == "K20P19_CORECFS":
            yr = range(1970,2024,5)
        elif release == "H16_CORECFS":
            yr = range(1970,2019,5)

        # Determine if variable is Level 1 or 2 based on variable name, and set 
        # internal variable name as needed
        
        if vname.endswith(('_bottom5m', '_surface5m', '_integrated', '_depthavg')):
            lev = 2
            vnamein = vname.replace('_bottom5m', "").replace('_surface5m',"").replace('_integrated',"").replace('_depthavg', "")
        elif (vname in ['uEast', 'vNorth']):
            lev = 2
        else:
            lev = 1

        # Coordinate dimension names
            
        ltname = "lat_rho"
        lnname = "lon_rho"
        tname = "ocean_time"

        # Build data file names

        fname = [f"{yyyy}-{yyyy+4}/B10K-{release}_{yyyy}-{yyyy+4}_average_{vname}.nc" for yyyy in yr]
        
        tmp = urlparse(simbase)
        if tmp.scheme in ['http', 'https']:
            fdata = ["/".join([simbase, f"B10K-{release}", f"Level{lev}", x]) for x in fname]
        else:
            fdata = [simbase / f"B10K-{release}" / f"Level{lev}" / x for x in fname]

        # Open dataset

        ds = xr.open_mfdataset(fdata, data_vars='all')
        ds_coord = xr.open_dataset(grdpath)

        ds['mask'] = ds_coord['mask_rho'] == 1
        ds_coord.close()

    elif model == "mom6nep":

        # Coordinate dimension names

        ltname = "geolat"
        lnname = "geolon"
        tname = "time"

        # Read dataset and coordinates

        if location == "TDS": 
            
            # Read from CEFI Data Portal. Use daily dataset if possible.  
            # If variable not found there, move to monthly

            C = cefiportalopts(portalpath='http://psl.noaa.gov/thredds/dodsC/Projects/CEFI/regional_mom6/cefi_portal/',
                               region="nep",
                               subdomain="full",
                               release=release)
            fname = C.cefifilelist(varname=vname, datestr="*")
            if not fname:
                C.freq = "monthly"
                fname = C.cefifilelist(varname=vname, datestr="*")

            sname = "/".join([C.cefifolder(), "ocean_static.nc"])

            ds_coord = xr.open_dataset(sname)
            if len(fname) > 1:
                ds_data  = xr.open_mfdataset(fname)
            else:
                ds_data = xr.open_dataset(fname[0])

            # Merge static variables into main dataset, label coordinate 
            # variables, and set the water mask

            ds_coord = ds_coord.set_coords(['geolon','geolat'])
            ds = xr.merge([ds_coord, ds_data])
            ds['mask'] = ds['wet'] == 1

            # Close files we're done with

            ds_coord.close()
            ds_data.close()

        elif location == "local":

            # Read from local copy of CEFI Data Portal (subsets for PEEC, ESR)

            C = cefiportalopts(portalpath='~/Documents/Data/CEFI/ak_cefiportal',
                               region="nep",
                               subdomain=[0, 342, 446, 743], 
                               release=release)
            fname = C.cefifilelist(varname=vname, datestr="*")
            C.freq="static"
            sname = C.cefifilelist(varname="ocean_static", datestr="*")

            ds_coord = xr.open_dataset(sname)
            ds_data  = xr.open_dataset(fname)

            # Merge static variables into main dataset, label coordinate 
            # variables, and set the water mask

            ds_coord = ds_coord.set_coords(['geolon','geolat'])
            ds = xr.merge([ds_coord, ds_data])
            ds['mask'] = ds['wet'] == 1

            # Close files we're done with

            ds_coord.close()
            ds_data.close()

        elif location == "EM_opendap":

            # Use daily dataset if possible.  If variable not found there, move 
            # to monthly

            fdaily = "https://compute.earthmover.io/v1/services/dap2/NOAA-PMEL/cefi-nep-hindcast-daily/add-statics/raw/main/opendap"
            fmonth = "https://compute.earthmover.io/v1/services/dap2/NOAA-PMEL/cefi-nep-hindcast-monthly/add-statics/raw/main/opendap"

            ds = xr.open_dataset(fdaily)
            if vname not in ds:
                ds.close()
                ds = xr.open_dataset(fmonth)

            # Label coordinate variables, set water mask
                
            ds = ds.set_coords(['geolon','geolat'])
            ds['mask'] = ds['wet'] == 1


        elif location == "EM_arraylake":

            # Arraylake protocol.  Use daily dataset if possible.  If 
            # variable not found there, move to monthly.

            client = Client()
            repo = client.get_repo("NOAA-PMEL/cefi-nep-hindcast-daily")
            session = repo.readonly_session(branch="add-statics")
            ds = xr.open_zarr(session.store, zarr_format=3, group="raw/main")

            if vname not in ds:
                ds.close()
                repo = client.get_repo("NOAA-PMEL/cefi-nep-hindcast-monthly")
                session = repo.readonly_session(branch="add-statics")
                ds = xr.open_zarr(session.store, zarr_format=3, group="raw/main")

            # Label coordinate variables, set water mask
                
            ds = ds.set_coords(['geolon','geolat'])
            ds['mask'] = ds['wet'] == 1


        # A few options use different variables in the file than in the file 
        # name.  For the rest, set to the same
            
        if vnamein == "":
            vnamein = vname
        
    # Read survey data (.csv) in dataset
            
    svy = pd.read_csv(fsurvey)
    ds_svy = xr.Dataset.from_dataframe(svy)

    # Mask out land in regional model output (we want to match closest 
    # ocean grid cell)

    ds[ltname] = ds[ltname].where(ds.mask)
    ds[lnname] = ds[lnname].where(ds.mask)

    ds[ltname] = ds[ltname].fillna(0) # works since 0,0 is far from region, but not universal
    ds[lnname] = ds[lnname].fillna(0)

    # Create index for quick distance searching using haversine distance
    #
    # Note: I'd rather be able to do this in projected coordinates, but I 
    # haven't figured out how to adapt the NDPointIndex for that.  Also, 
    # possibly consider a distance-with-barriers metric if we apply this to 
    # the Aleutians?

    ds = ds.set_xindex(
        (ltname, lnname), 
        xr.indexes.NDPointIndex, 
        tree_adapter_cls=SklearnGeoBallTreeAdapter)


    # Select closet points in horizontal space

    dimsel = {ltname: ds_svy['latitude'], 
              lnname: ds_svy['longitude'],
               tname: ds_svy['start_time'],
            'method': 'nearest'}
    ds_nearest = ds.sel(**dimsel)

    # Select specified vertical layer or value

    if zname:
        if isinstance(zsel,int): # index of layer
            if (zsel==-1): # special case, bottom-most non-NaN value
                ds_nearest[vnamein] = ds_nearest[vnamein].dropna(**{'dim': zname}).isel(**{zname: -1})
            elif zsel >= 0:
                ds_nearest[vnamein] = ds_nearest[vnamein].isel(**{zname: zsel})
        else: # depth value
            ds_nearest[vnamein] = ds_nearest[vnamein].sel(**{zname: zsel})

    # Extract values of closest point, save to file

    svy[label] = ds_nearest[vnamein].values
    svy.to_csv(fout, index=False)

    # Close dataset
    
    ds.close()