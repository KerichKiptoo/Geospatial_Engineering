/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var aoi = ee.FeatureCollection("projects/my-project-cnn-334706/assets/Victoria");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
/*
========================================================================================
5. COMPUTE OTHER RS METRICS FROM MASKED DATA SETS
========================================================================================
*/

var s2 = ee.ImageCollection("COPERNICUS/S2_SR");


Map.centerObject(aoi)
// Write a function for Cloud masking
function maskS2clouds(image) {
  var qa = image.select('QA60')
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0).and(
             qa.bitwiseAnd(cirrusBitMask).eq(0))
  return image.updateMask(mask).divide(10000)
      .select("B.*")
      .copyProperties(image, ["system:time_start"])
}
var nocloud =  s2
  .filter(ee.Filter.date('2021-01-01', '2021-12-31'))
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 30))
  .filter(ee.Filter.bounds(aoi))
  .map(maskS2clouds);


function addIndices(image){
  var rangeBands = {
                      'blue': image.select('B2'),
                      'green':image.select('B3'),
                      'red': image.select('B4'),
                      'nir': image.select('B8'),
                      'swir1': image.select('B11')
                    };
                   
  //NDVI
  var ndvi = image.expression({
    expression: '(nir - red) / (nir + red)',
    map: rangeBands
  });
  
  //chlorophyll
  var chl = image.expression({
    expression: '(nir/green)-1',
    map: rangeBands
  });
  
  //SAVI
 var savi = image.expression({
   expression: '((nir-red)/(nir+red+0.5))*1.5',
   map: rangeBands
 });
 //NDWI
 var ndwi = image.expression({
   expression: '(green - nir) / (green + nir)',
   map: rangeBands
 });
 
  image = image.addBands(ndvi.clip(aoi).rename('ndvi'))
               .addBands(savi.clip(aoi).rename('savi'))
               .addBands(chl.clip(aoi).rename('chl'))
               .addBands(ndwi.clip(aoi).rename('ndwi'));
  return image.float();
}


//Compute indices and add to image median image composite
var image = nocloud.map(addIndices);
var indices = image.select('ndvi','savi','ndwi', 'chl');
print("Indices ", indices.first().bandNames());

var savi = indices.select('savi').median()
var ndvi = indices.select('ndvi').median()
var ndwi = indices.select('ndwi').median()
var chl = indices.select('chl').median()
print(chl)


Map.addLayer(chl)
// Map.addLayer(savi)
// Map.addLayer(ndwi)

Export.image.toDrive({
    image: ndvi,
    description: 'ndvi',
    folder: 'earthengine',
    fileNamePrefix: 'ndvi',
    region: aoi,
    maxPixels: 1e13,
    scale: 30,
})

Export.image.toDrive({
    image: ndvi,
    description: 'savi',
    folder: 'earthengine',
    fileNamePrefix: 'savi',
    region: aoi,
    maxPixels: 1e13,
    scale: 30,
})

Export.image.toDrive({
    image: ndvi,
    description: 'ndwi',
    folder: 'earthengine',
    fileNamePrefix: 'ndwi',
    region: aoi,
    maxPixels: 1e13,
    scale: 30,
})
