library(akgfmaps)
library(ncdf4)
 
shapefol <- "./akgfmaps_shapefiles"
if (! dir.exists(shapefol)) {
  dir.create(shapefol)
}

# Groundfish bottom trawl survey regions

sebs <- akgfmaps::get_base_layers(select.region = "bs.south", set.crs = "auto")
nebs <- akgfmaps::get_base_layers(select.region = "bs.north", set.crs = "auto")
goa  <- akgfmaps::get_base_layers(select.region = "goa",      set.crs = "auto")
ecs  <- akgfmaps::get_base_layers(select.region = "ecs",      set.crs = "auto")
ai   <- akgfmaps::get_base_layers(select.region = "ai",       set.crs = "auto")
  
sebs$survey.area <- st_simplify(sebs$survey.area, preserveTopology = TRUE, dTolerance = 1000)
nebs$survey.area <- st_simplify(nebs$survey.area, preserveTopology = TRUE, dTolerance = 1000)
goa$survey.area  <- st_simplify(goa$survey.area,  preserveTopology = TRUE, dTolerance = 1000)
ecs$survey.area  <- st_simplify(ecs$survey.area,  preserveTopology = TRUE, dTolerance = 1000)
ai$survey.area   <- st_simplify(ai$survey.area,   preserveTopology = TRUE, dTolerance = 1000)

if (! file.exists(file.path(shapefol, 'surveyarea_goa.shp'))) {
  st_write(goa$survey.area,  file.path(shapefol, 'surveyarea_goa.shp'))
}
if (! file.exists(file.path(shapefol, 'surveyarea_sebs.shp'))) {
  st_write(sebs$survey.area,  file.path(shapefol, 'surveyarea_sebs.shp'))
}
if (! file.exists(file.path(shapefol, 'surveyarea_nebs.shp'))) {
  st_write(nebs$survey.area,  file.path(shapefol, 'surveyarea_nebs.shp'))
}
if (! file.exists(file.path(shapefol, 'surveyarea_ecs.shp'))) {
  st_write(ecs$survey.area,  file.path(shapefol, 'surveyarea_ecs.shp'))
}
if (! file.exists(file.path(shapefol, 'surveyarea_ai.shp'))) {
  st_write(ai$survey.area,  file.path(shapefol, 'surveyarea_ai.shp'))
}

# Bering Sea survey strata

if (! file.exists(file.path(shapefol, 'surveystrata_ebs.shp'))) {
  st_write(sebs$survey.strata,  file.path(shapefol, 'surveystrata_ebs.shp'))
}

# ESR subareas

esr_subareas <- get_esr_regions(select.region = "esr_subarea", set.crs = "auto")
esr_subareas <- st_simplify(esr_subareas, preserveTopology = TRUE, dTolerance = 1000)

if (! file.exists(file.path(shapefol, 'esr_subareas.shp'))) {
  st_write(esr_subareas, file.path(shapefol, 'esr_subareas.shp'))
}



