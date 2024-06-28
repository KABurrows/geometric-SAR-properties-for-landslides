//With this script you can calculate geometric properties of landslide properties and how they relate to Sentinel-1 SAR images in Google Earth Engine
//Landslides that face directly towards the SAR sensor may be susceptible to geometric decorrelation when forming interferograms with long perpendicular baselines (see e.g. Lee and Liu, 1999, DOI: IEEE.10.1109/IGARSS.1999.773541)
//Landslides in areas of steep topography may be distorted in SAR amplitude images
//This script does not currently estimate areas of foreshortening and layover, which may also be present in SAR images acquired over areas of steep topography

//read in landslide inventory
//the landslide inventory used here (Nepal_monsoon_2019_2000) is available at https://zenodo.org/records/7970874
var landslides = ee.FeatureCollection('path/to/Nepal_monsoon_2019_2000')
Map.addLayer(landslides)
Map.centerObject(landslides.first(),9);

//polygon outlining extent of landslide inventory (can just draw in GEE)
var AOI=geometry

//read in ALOS DEM
var dataset = ee.ImageCollection('JAXA/ALOS/AW3D30/V3_2');
var elevation = dataset.select('DSM');
//visualise elevation dataset
var elevationVis = {
  min: 0,
  max: 5000,
  palette: ['0000ff', '00ffff', 'ffff00', 'ff0000', 'ffffff']
};
Map.addLayer(elevation, elevationVis, 'Elevation');

// Calculate slope and aspect
// Reproject an image mosaic using a projection from one of the image tiles, rather than using the default projection returned by .mosaic().
// This is necessary for ALOS and Copernicus DEMs, but for SRTM you can skip this and load the DEM directly.
var proj = elevation.first().select(0).projection();
var slopeReprojected = ee.Terrain.slope(elevation.mosaic()
                             .setDefaultProjection(proj));
Map.addLayer(slopeReprojected, {min: 0, max: 45}, 'Slope')
var aspectReprojected = ee.Terrain.aspect(elevation.mosaic()
                             .setDefaultProjection(proj));
Map.addLayer(aspectReprojected, {min: 0, max: 360}, 'Aspect')

//read in stack of Sentinel-1 images
var image_s1 = ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterBounds(AOI)
    .filterDate('2018-04-01','2018-05-01')//Filter by date (although we only need one SAR image for the calculation so we will just select the first one in the series at line 44 of the script)
    .filter(ee.Filter.eq('relativeOrbitNumber_start',19))//Filter by orbit number. If you don't know this, test which images cover your AOI by loading the SAR data without this line
    .filter(ee.Filter.and(
      ee.Filter.eq('instrumentMode', 'IW')
    )).first()//for now, we only need one SAR image to get the look direction etc.
    
// get Sentinel-1 image geometry and projection
var geom_s1 = image_s1.geometry();
var proj_s1 = image_s1.select(1).projection();
print(geom_s1)
// get look direction angle
var look_dir = (ee.Terrain.aspect(image_s1.select('angle')).reduceRegion(ee.Reducer.mean(), geom_s1, 1000).get('aspect'));//
print(look_dir)

// Radar geometry
var theta_iRad = image_s1.select('angle').multiply(Math.PI/180).clip(geom_s1);//incidence angle of SAR
var phi_iRad = ee.Image.constant(look_dir).multiply(Math.PI/180);//look direction in radians (flat image)
Map.addLayer(theta_iRad)

// Cut topography images and reproject to match SAR images
var alpha_sRad = slopeReprojected.multiply(Math.PI/180).setDefaultProjection(proj_s1).clip(geom_s1);//topographic slope in radians
var phi_sRad = aspectReprojected.multiply(Math.PI/180).setDefaultProjection(proj_s1).clip(geom_s1);//topographic aspect in radians

//calculate slope steepnesses in range (and azimuth) directions
//slope aspect relative to satellite look direction
var phi_rRad = phi_iRad.subtract(phi_sRad);//satellite look direction in radians - local aspect
Map.addLayer(phi_rRad)

// slope steepness in range
var alpha_rRad = (alpha_sRad.tan().multiply(phi_rRad.cos())).atan();
Map.addLayer(alpha_rRad)

//slope steepness in azimuth - commented out because it's not quite as useful.
//var alpha_azRad = (alpha_sRad.tan().multiply(phi_rRad.sin())).atan();    

//Difference between SAR incidence angle and local slope in LOS. Geometric decorrelation is highest as this value approaches zero.
var theta_iRad_loc = theta_iRad.subtract(alpha_rRad)
Map.addLayer(theta_iRad_loc,{min:0,max:1.57})

//Next step - extract the geometric properties for every landslide so we can download the factor and calculate it for every landslide and every Bperp and thus do a correction on every 
//get and export the SAR incidence angle for every landslide polygon
var ls_theta = theta_iRad.reduceRegions({
  collection: landslides.select(['object_id']),//object_id is the name of an attribute in the shapefile labelling every landslide (landslide1, landslide2 etc.)
  reducer: ee.Reducer.median(),
  scale:10
})
Export.table.toDrive({
  collection: ls_theta,
  description: 'ls_theta_T019D',
  fileFormat: 'CSV',
  folder: 'foldername'
});

//get and export the slope in Line of Sight direction for every landslide polygon
var ls_alpha = alpha_rRad.reduceRegions({
  collection: landslides.select(['object_id']),
  reducer: ee.Reducer.median(),
  scale:10
})

Export.table.toDrive({
  collection: ls_alpha,
  description: 'ls_slope_in_LOS_T019D',
  fileFormat: 'CSV',
  folder: 'foldername'
});

//export topographic slope of every landslide polygon
var ls_slope=alpha_sRad.reduceRegions({
  collection: landslides.select(['object_id']),
  reducer: ee.Reducer.median(),
  scale:10
});
Export.table.toDrive({
  collection: ls_slope,
  description: 'ls_slopes',
  fileFormat: 'CSV',
  folder: 'foldername'
})

//export aspect of every landslide polygon
var ls_aspect=phi_sRad.reduceRegions({
  collection: landslides.select(['object_id']),
  reducer: ee.Reducer.median(),
  scale:10
});
Export.table.toDrive({
  collection: ls_slope,
  description: 'ls_aspects',
  fileFormat: 'CSV',
  folder: 'foldername'
})

//Export the difference between incidence angle and local slope in LOS for every landslide
var ls_theta_alpha_difference = theta_iRad_loc.reduceRegions({
  collection: landslides.select(['object_id']),//object_id is the name of an attribute in the shapefile labelling every landslide (landslide1, landslide2 etc.)
  reducer: ee.Reducer.median(),
  scale:10
})
Export.table.toDrive({
  collection: ls_theta,
  description: 'ls_theta_alpha_difference_T019D',
  fileFormat: 'CSV',
  folder: 'foldername'
});

//Estimate geometric decorrelation for a given perpendicular baseline
//Note that this will almost certainly UNDERESTIMATE decorrelation particularly in vegetated areas where volume decorrelation also needs to be taken into account (See e.g. Hoen and Zebker, 2000 DOI: 10.1109/36.885204)
//decorrelation is also caused by temporal factors such as movement/growth of vegetation, changes in soil moisture which are not included in this estimate.

//generate flat images of constant values
var Bperp=150//perpendicular baseline in m (150 m chosen as an example to test the result)
var c=299792458//Speed of light
var c_im = ee.Image.constant(c).setDefaultProjection(proj_s1).clip(geom_s1);//generate image of c
var lambda = ee.Image.constant(0.05547).setDefaultProjection(proj_s1).clip(geom_s1)//generate image of Sentinel-1 wavelength (in m)
var r = ee.Image.constant(693000).setDefaultProjection(proj_s1).clip(geom_s1)//generate image of earth-satellite distance
var Bw = ee.Image.constant(48300000).setDefaultProjection(proj_s1).clip(geom_s1)//chirp bandwidth - 42.8 MHz correct for IW3
var Bp = ee.Image.constant(Bperp).setDefaultProjection(proj_s1).clip(geom_s1)
var ones = ee.Image.constant(1.0).setDefaultProjection(proj_s1).clip(geom_s1)
//geometric coherence can be estimated from 1-c.Bperp/(lambda.r.Bw.tan(theta-alpha)) - Lee and Liu (1999) DOI: 10.1109/IGARSS.1999.773541
var coh_geom = ones.subtract(c_im.multiply(Bp).divide(lambda.multiply(r).multiply(theta_iRad_loc.tan().abs()).multiply(Bw)))
Map.addLayer(coh_geom,{min:0,max:1.0},'modelled geometric coherence')//

//Export the modelled geometric coherence for every landslide
//Note - be careful how you use this, especially if you are working in a vegetated area!
var ls_modelled_coh_geom = coh_geom.reduceRegions({
  collection: landslides.select(['object_id']),//object_id is the name of an attribute in the shapefile labelling every landslide (landslide1, landslide2 etc.)
  reducer: ee.Reducer.median(),
  scale:10
})
Export.table.toDrive({
  collection: ls_modelled_coh_geom,
  description: 'ls_modelled_coh_geom_T019D',
  fileFormat: 'CSV',
  folder: 'foldername'
});