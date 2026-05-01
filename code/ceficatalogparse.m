function ceficatalogparse(pth)

arguments
    pth {mustBeTextScalar} ="http://psl.noaa.gov/thredds/catalog/Projects/CEFI/regional_mom6/cefi_portal/"
end

tree = parseXML(fullfile(pth, 'catalog.xml');

