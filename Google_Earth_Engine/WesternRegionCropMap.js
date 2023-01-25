/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var s2 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2"),
    ref = ee.FeatureCollection("users/Kerich/Kenya_crops"),
    aoi = ee.FeatureCollection("users/Kerich/Kenya_shapefile");
/***** End of imports. If edited, may not auto-convert in the playground. *****/



Map.addLayer(aoi, {color: 'red'}, 'Farm')
Map.centerObject(aoi)
var rgbVis = {min: 0.0, max: 3000, bands: ['B4', 'B3', 'B2']};
var minLabel = 1;
var maxLabel = 10;
var resolution = 30;

var Class = "Landcover1"; //Land-cover attribute in the shapefile
var filtered = s2
  .filter(ee.Filter.date('2015-03-01', '2015-10-31'))
  //.filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 30))
  .filter(ee.Filter.bounds(aoi));
  var filtered1 = s2
  .filter(ee.Filter.date('2016-03-01', '2016-10-31'))
  //.filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 30))
  .filter(ee.Filter.bounds(aoi));

// Write a function for Cloud masking
//Cloud Mask
function CloudMask(image) {
  // Bits 3 and 5 are cloud shadow and cloud, respectively.
  var cloudShadowBitMask = (1 << 3);
  var cloudsBitMask = (1 << 5);
  // Get the pixel QA band.
  var qa = image.select('QA_PIXEL');
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
                 .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  return image.updateMask(mask);
}


var filtered = filtered.map(CloudMask)
print(filtered,'filt')
var filtered2 = filtered1.map(CloudMask)


function addIndices(image){
  var rangeBands = {
                      'blue': image.select('SR_B2'),
                      'green':image.select('SR_B3'),
                      'red': image.select('SR_B4'),
                      'nir': image.select('SR_B5'),
                      'swir1': image.select('SR_B6')
                    };
                   
  //NDVI
  var ndvi = image.expression({
    expression: '(nir - red) / (nir + red)',
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
  //RVI
  var rvi = image.expression({
    expression : '(red - nir)',
    map: rangeBands
  });
  
  //IPVI
  var ipvi = image.expression({
    expression : '(nir/(nir+red)',
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
               .addBands(sipi.clip(aoi).rename('sipi'));
  return image.float();
}

var image = filtered.map(addIndices);
var image_2016 = filtered1.map(addIndices);
// Moving-Window Smoothing

// Specify the time-window
var days = 32

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
var ndviCol = image.select('ndvi','gli','msavi','ndmi','gci','evi')
var ndviCol2 = image_2016.select('ndvi','gli','msavi','ndmi','gci','evi')

var joinedCollection = join.apply({
  primary: ndviCol, 
  secondary: ndviCol, 
  condition: diffFilter
});

var joinedCollection2 = join.apply({
  primary: ndviCol2, 
  secondary: ndviCol2, 
  condition: diffFilter
});
// Each image in the joined collection will contain
// matching images in the 'images' property
// Extract and return the mean of matched images
var NDVIsmoothedCollection = ee.ImageCollection(joinedCollection.map(function(image) {
  var collection = ee.ImageCollection.fromImages(image.get('images'));
  return ee.Image(image).addBands(collection.median().rename('median_nd','median_gli','median_msavi','median_ndmi','median_gci','median_evi'));
}))
  
var NDVIsmoothedCollection2016 = ee.ImageCollection(joinedCollection2.map(function(image) {
  var collection2 = ee.ImageCollection.fromImages(image.get('images'));
  return ee.Image(image).addBands(collection2.median().rename('median_nd','median_gli','median_msavi','median_ndmi','median_gci','median_evi'));
}))
  
var median = NDVIsmoothedCollection.select('median.*');
var median2016 = NDVIsmoothedCollection2016.select('median.*');


var mediancomp =median.toBands()
var mediancomp2016 =median2016.toBands()

var names = ee.List.sequence(0,mediancomp.bandNames().size().subtract(1))
  .map(function(element){
    return ee.String('indice').cat(ee.Number(element).toInt())
  });
  
var names2 = ee.List.sequence(0,mediancomp2016.bandNames().size().subtract(1))
  .map(function(element){
    return ee.String('indice').cat(ee.Number(element).toInt())
  });
var BandsRenamed = mediancomp.rename(names);
var BandsRenamed2016 = mediancomp2016.rename(names2);

print(BandsRenamed2016,'renamed2');

var dates = NDVIsmoothedCollection.map(function(image) {
      return ee.Feature(null, {'date': image.date().format('YYYY-MM-dd')})
    })
    .distinct('date')
    .aggregate_array('date');
    
print (dates,'dates')


var bandnamesy = mediancomp.bandNames()
print(bandnamesy,'MEED')

var interval = 16;
var increment = 'day';
var start2 = '2015-03-01';
// make a list of start years
var startDate2 = ee.Date(start2);
var secondDate2 = startDate2.advance(interval, increment).millis();
var increase2 = secondDate2.subtract(startDate2.millis());
var list2 = ee.List.sequence(startDate2.millis(), ee.Date('2015-10-31').millis(), increase2);

// 16 day composite
var tenDaysmed2 =  ee.ImageCollection.fromImages(list2.map(function(date){
  return image.select('ndvi','evi','msavi','ndmi','gli').filterDate(ee.Date(date), ee.Date(date).advance(interval, increment))
           .median().set('system:time_start',ee.Date(date).millis());
           
}));

var inputFeatures = ['ndvi','evi','msavi','ndmi','gli','gci']
var inputSeries = tenDaysmed2.select(inputFeatures);
var inputImage = (tenDaysmed2.toArray().arrayFlatten([
  ee.List.sequence(1, inputSeries.size()).map(function (it) { return ee.String(ee.Number(it).int()) }),
  inputFeatures
]))

;

var tenDaysmed3 =tenDaysmed2.toBands();
//print(tenDaysmed3,'tenDaysmed3')
var inputBands = inputImage.bandNames();
var inputImage = inputImage.select(inputBands);
//print(inputImage,'ag')

//print(tenDaysmed2,'ten')
//print(ui.Chart.image.series(tenDaysmed2,aoi,ee.Reducer.mean(),30))


// classification 
// var datawithColumn = ref.randomColumn('random');
// print(datawithColumn,'datawith')
// var split = 0.7; // separate 70% for training, 30% for validation
// var trainingData = datawithColumn.filter(ee.Filter.lt('random', split));
// var validationData = datawithColumn.filter(ee.Filter.gte('random', split));
// print(trainingData)
var bandnames = BandsRenamed.bandNames(); //All bands on included here
//print(bandnames,'bnanes')
var points = BandsRenamed.select(bandnames).sampleRegions({
  collection: ref,
  properties: ['Landcover1'],
  scale: 30 
  }).randomColumn();
  
// var split = 0.7; // separate 70% for training, 30% for validation
// var trainingData = datawithColumn.filter(ee.Filter.lt('random', split));
// var validationData = datawithColumn.filter(ee.Filter.gte('random', split));
var training = points.filter(ee.Filter.lt('random', 0.7)); 
var validation = points.filter(ee.Filter.gte('random', 0.7));
print(points.size());
print(training.size());
print(validation.size());
//print(training,'training')
var classifier = ee.Classifier.smileRandomForest(10)
    .train({
    features: training,
    classProperty: Class,
    inputProperties: bandnames
    });


var classified = BandsRenamed.select(bandnames).classify(classifier);
Map.addLayer(classified, {min: minLabel, max: maxLabel, palette: ['008080', '54FF9F', '40870A', 'FFFF00', 'ddba22', '76B485', 'FF00FF', '#E96B27', '0000FF', 'BEBEBE']},'classified');

//Accuracy Assessment
// Get a confusion matrix representing resubstitution accuracy.
var trainAccuracy = classifier.confusionMatrix();
print('Resubstitution error matrix: ', trainAccuracy);
print('Training overall accuracy: ', trainAccuracy.accuracy());

// Classify the validation data.

// var validations = tenDaysmed3.select(bandnames).sampleRegions({
//   collection: validation,
//   properties: ['class'],
  
//   scale: 30 
//   });
var validated = validation.classify(classifier);

// Get a confusion matrix representing expected accuracy.
var testAccuracy = validated.errorMatrix('Landcover1', 'classification');
print('Validation error matrix: ', testAccuracy);
print('Validation overall accuracy: ', testAccuracy.accuracy());

var OA = testAccuracy.accuracy()
var CA = testAccuracy.consumersAccuracy()
var Kappa = testAccuracy.kappa()
var Order = testAccuracy.order()
var PA = testAccuracy.producersAccuracy()
print(OA,'Overall Accuracy')
print(CA,'Consumers Accuracy')
print(Kappa,'Kappa')
print(Order,'Order')
print(PA,'Producers Accuracy')


//2020 CLASSIFICATION
var image2 = filtered2.map(addIndices);

var interval = 16;
var increment = 'day';
var start3 = '2016-03-01';
// make a list of start years
var startDate3 = ee.Date(start3);
var secondDate3 = startDate3.advance(interval, increment).millis();
var increase3 = secondDate3.subtract(startDate3.millis());
var list3 = ee.List.sequence(startDate3.millis(), ee.Date('2016-10-31').millis(), increase3);

// 16 day composite
var tenDaysmed4 =  ee.ImageCollection.fromImages(list3.map(function(date){
  return image2.select('ndvi','evi','msavi','ndmi','gli').filterDate(ee.Date(date), ee.Date(date).advance(interval, increment))
           .median().set('system:time_start',ee.Date(date).millis());
           
}));

print(ui.Chart.image.series(tenDaysmed4,aoi,ee.Reducer.mean(),30))
var tenDaysmed5 =tenDaysmed4.toBands();
var bandnames2 = tenDaysmed5.bandNames(); //All bands on included here
print(bandnames2,'bnanes2')
print(tenDaysmed5,'tenDaysmed5')

var classified_2016 = BandsRenamed2016.select(bandnames).classify(classifier);
Map.addLayer(classified_2016, {min: minLabel, max: maxLabel, palette: ['008080', '54FF9F', '40870A', 'FFFF00', 'ddba22', '76B485', 'FF00FF', '#E96B27', '0000FF', 'BEBEBE']},'classified 2016');


var points_2016 = tenDaysmed5.select(bandnames).sampleRegions({
  collection: ref,
  properties: ['Landcover1'],
  scale: 30 
  }).randomColumn();


var training_2016 = points_2016.filter(ee.Filter.lt('random', 0.7)); 
var validation_2016 = points_2016.filter(ee.Filter.gte('random', 0.7));

// prepare data for DTW

var NDVI_TS_2015 =  ee.ImageCollection.fromImages(list2.map(function(date){
  return image.select('ndvi').filterDate(ee.Date(date), ee.Date(date).advance(interval, increment))
           .median().set('system:time_start',ee.Date(date).millis());
           
}));
var NDVI_TS_20155 =NDVI_TS_2015.toBands();
//print (NDVI_TS_20155,'NDVI_TS')

var EVI_TS_2015 =  ee.ImageCollection.fromImages(list2.map(function(date){
  return image.select('evi').filterDate(ee.Date(date), ee.Date(date).advance(interval, increment))
           .median().set('system:time_start',ee.Date(date).millis());
           
}));
var EVI_TS_20155 =EVI_TS_2015.toBands();
//print (EVI_TS_20155,'EVI_TS')

var MSAVI_TS_2015 =  ee.ImageCollection.fromImages(list2.map(function(date){
  return image.select('msavi').filterDate(ee.Date(date), ee.Date(date).advance(interval, increment))
           .median().set('system:time_start',ee.Date(date).millis());
           
}));
var MSAVI_TS_20155 =MSAVI_TS_2015.toBands();
//print (MSAVI_TS_20155,'MSAVI_TS')

var NDMI_TS_2015 =  ee.ImageCollection.fromImages(list2.map(function(date){
  return image.select('ndmi').filterDate(ee.Date(date), ee.Date(date).advance(interval, increment))
           .median().set('system:time_start',ee.Date(date).millis());
           
}));
var NDMI_TS_20155 =NDMI_TS_2015.toBands();
//print (NDMI_TS_20155,'NDMI_TS')

var GLI_TS_2015 =  ee.ImageCollection.fromImages(list2.map(function(date){
  return image.select('gli').filterDate(ee.Date(date), ee.Date(date).advance(interval, increment))
           .median().set('system:time_start',ee.Date(date).millis());
           
}));
var GLI_TS_20155 =GLI_TS_2015.toBands();
print (GLI_TS_20155,'GLI_TS')

Export.image.toDrive({ 
  image: classified,
  description: 'classified_2015',
  scale: 30, 
  maxPixels: 1e13, 
  region: aoi 
});



