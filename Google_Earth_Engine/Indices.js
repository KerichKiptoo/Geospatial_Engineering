/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var aoi = ee.FeatureCollection("projects/my-project-cnn-334706/assets/Victoria");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
/*
========================================================================================
5. COMPUTE OTHER RS METRICS FROM MASKED DATA SETS
========================================================================================
*/
//var aoi = ee.FeatureCollection("users/kimeli/Bbound_Busia");
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

// generte 15 day medinan composites
/*var interval = 10;
var increment = 'day';
var start = '2019-03-01';
// make a list of start years
var startDate = ee.Date(start);
var secondDate = startDate.advance(interval, increment).millis();
var increase = secondDate.subtract(startDate.millis());
var list = ee.List.sequence(startDate.millis(), ee.Date('2019-10-31').millis(), increase);

// ten days  composite
var tenDaysmed =  ee.ImageCollection.fromImages(list.map(function(date){
  return nocloud.select('B2','B3','B4','B8','B11').filterDate(ee.Date(date), ee.Date(date).advance(interval, increment))
           .mean().set('system:time_start',ee.Date(date).millis());
           
}));
//print(ui.Chart.image.series(tenDaysmed,aoi,ee.Reducer.mean(),30))
//print (tenDaysmed)

*/
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
 var savi = image.expression({
   expression: '((nir-red)/(nir+red+0.5))*1.5',
   map: rangeBands
 });
 
 var ndwi = image.expression({
   expression: '(green - nir) / (green + nir)',
   map: rangeBands
 });
 

  //Green Normalized difference Vegetation Index (GNDVI)
   var gdnvi =  image.expression({
     expression: '(nir - green) / (nir + green)',
     map: rangeBands
   });

  //EVI
  var evi = image.expression({
    expression: '2.5 * ((nir - red) / (nir + (6 * red) - (7.5 * blue) + 1))',
    map: rangeBands
  });
 
  //Normalize Differennce Mositure Index (NDMI)
  var ndmi = image.expression({
    expression: '(nir - swir1) / (nir + swir1)',
     map: rangeBands
  });
  // normalized special index (NSI)
  var nsi = image.expression({
    expression: '100 * ((swir1 - nir) / (nir - red))',
     map: rangeBands
  });
  //Modified Soil Adjusted Vegetation Index MSAVI
  var msavi =image.expression({
    expression: '(2 * nir + 1 - sqrt(pow((2 * nir + 1), 2) - 8 * (nir - red)) ) / 2',
    map: rangeBands
  });
  //Green Chlorophyll Index (GCI)
  var gci =image.expression({
    expression: '(nir/green) - 1',
    map: rangeBands
  });
  //Green Leaf Index (GLI)
  var gli =image.expression({
    expression: '((green - red) + (green - blue)) / ((2 * green) + red + blue)',
    map: rangeBands
  });
  //Structure Insensitive Pigment Index (sipi)
  var sipi =image.expression({
    expression: '(nir - blue) / (nir - red)',
    map: rangeBands
  });
  image = image.addBands(ndvi.clip(aoi).rename('ndvi'))
               .addBands(gdnvi.clip(aoi).rename('gdnvi'))
               .addBands(evi.clip(aoi).rename('evi'))
               .addBands(ndmi.clip(aoi).rename('ndmi'))
               .addBands(nsi.clip(aoi).rename('nsi'))
               .addBands(msavi.clip(aoi).rename('msavi'))
               .addBands(gci.clip(aoi).rename('gci'))
               .addBands(gli.clip(aoi).rename('gli'))
               .addBands(savi.clip(aoi).rename('savi'))
               .addBands(ndwi.clip(aoi).rename('ndwi'))
               .addBands(sipi.clip(aoi).rename('sipi'));
  return image.float();
}


//Compute indices and add to image median image composite
var image = nocloud.map(addIndices);
var indices = image.select('ndvi','gdnvi','evi','ndmi','nsi','msavi','gci','gli','sipi','savi','ndwi');
print("Indices ", indices.first().bandNames());

var savi = indices.select('savi').median()
var ndvi = indices.select('ndvi').median()
var ndwi = indices.select('ndwi').median()
print(ndwi)


Map.addLayer(ndvi)

                                    // Moving-Window Smoothing

// Specify the time-window
/*var days = 30

// We use a 'join' to find all images that are within
// the time-window
// The join will add all matching images into a
// new property called 'images'
var join = ee.Join.saveAll({
  matchesKey: 'images'
});

// This filter will match all images that are captured
// within the specified day of the source image
var diffFilter = ee.Filter.maxDifference({
  difference: 1000 * 60 * 60 * 24 * days,
  leftField: 'system:time_start', 
  rightField: 'system:time_start'
});

// Select the 'ndvi' band
var ndviCol = indices.select('ndvi')

var joinedCollection = join.apply({
  primary: ndviCol, 
  secondary: ndviCol, 
  condition: diffFilter
});

// Each image in the joined collection will contain
// matching images in the 'images' property
// Extract and return the mean of matched images
var smoothedCollection = ee.ImageCollection(joinedCollection.map(function(image) {
  var collection = ee.ImageCollection.fromImages(image.get('images'));
  return ee.Image(image).addBands(collection.median().rename('moving_median'));
}))
  
// Display a time-series chart
var chart = ui.Chart.image.series({
  imageCollection: smoothedCollection.select(['ndvi', 'moving_median']),
  region: aoi,
  reducer: ee.Reducer.mean(),
  scale: 20
}).setOptions({
      lineWidth: 1,
      title: 'NDVI Time Series',
      interpolateNulls: true,
      vAxis: {title: 'NDVI'},
      hAxis: {title: '', format: 'YYYY-MMM'},
      series: {
        1: {color: 'gray', lineDashStyle: [1, 1]}, // Original NDVI
        0: {color: 'red', lineWidth: 2 }, // Smoothed NDVI
      },

    })
print(chart);
var ndvi_med_ts =smoothedCollection.select("moving_median")
var list = ee.List([0, 3, 6, 9, 12, 15,18,21,24,27,30]);
print(ndvi_med_ts)
print(ee.List.sequence(0, 30,3));
*/
