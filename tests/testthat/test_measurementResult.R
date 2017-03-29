################################################################
# test detection                                               #
################################################################


#!#!#!# ERROR with subsetSimulations, run that test first (multiple, unsorted locations)
# it seems: multiple plumes are returned multiple times, multiple locations ar not (as subsetSDF does not do it)
test_that("detection", {
  data(radioactivePlumes_area)

  # sample locations
  sampleLocations1 = sample.int(nLocations(radioactivePlumes_area), 1000)
  locations1 = c(0, 3000, NA, sampleLocations1, rev(sampleLocations1))
  i = sample.int(nPlumes(radioactivePlumes_area), 1)
  # generate intermediate results
  ## layer: detection
  threshold = 1e-7
  radioactivePlumes_area@values$detectable = calc(
    radioactivePlumes_area@values$maxdose,
    fun = function(x){x >= threshold})
  ## plumes:
  ### total dose
  radioactivePlumes_area@plumes$totalDose =
    summaryPlumes(radioactivePlumes_area, fun = sum, kinds = "finaldose")[[2]] # not area-weighted
  ### number of locations where it can be detected
  radioactivePlumes_area@plumes$nDetectable =
    summaryPlumes(radioactivePlumes_area, fun = sum, kinds = "detectable")[[2]]
  ### earliest possible detection of plume inside the area
  radioactivePlumes_area@plumes$earliestDetection =
    summaryPlumes(radioactivePlumes_area, fun = min, kinds = "time", na.rm = TRUE)[[2]]

  expect_error(# simulationsApply: locations invalid
    singleDetection(
      simulations = radioactivePlumes_area,
      locations = locations1)
  )

  # single detection
  ## no sensors
  singleDetection0 = singleDetection(
    simulations = radioactivePlumes_area,
    locations = integer(0))

  expect_equal(
    singleDetection0[["cost"]], 1)
  expect_equivalent(
    singleDetection0[["costPlumes"]], matrix(1, nrow = nPlumes(radioactivePlumes_area), ncol = 1))

  ## all sensors
  singleDetectionAll = singleDetection(
    simulations = radioactivePlumes_area,
    locations = 1:nLocations(radioactivePlumes_area))

  expect_equal(
    singleDetectionAll[["cost"]], 0)
  expect_equivalent(
    singleDetectionAll[["costPlumes"]], matrix(0, nrow = nPlumes(radioactivePlumes_area), ncol = 1))

  ## some sensors
  singleDetection1 = singleDetection(
    simulations = radioactivePlumes_area,
    locations = sampleLocations1)

  expect_equal(
    singleDetection1[["costPlumes"]][i],
    prod(1 - getValues(subset(radioactivePlumes_area,
               locations = sampleLocations1, plumes = i, kinds = "detectable", valuesOnly = TRUE)))
  )
  expect_equal(
    singleDetection1[["cost"]],
    mean(singleDetection1[["costPlumes"]] *
           radioactivePlumes_area@plumes$totalDose)/
      mean(radioactivePlumes_area@plumes$totalDose)
  )

  # multiple detection
  ## no sensors
  multipleDetection0 = multipleDetection(
    simulations = radioactivePlumes_area,
    locations = integer(0))

  expect_equal(
    multipleDetection0[["cost"]], 1)
  expect_equivalent(
    multipleDetection0[["costPlumes"]],
    matrix(0, nrow = nPlumes(radioactivePlumes_area), ncol = 1)) # 0: number of detections, not cost for single plume

  ## all sensors
  multipleDetectionAll = multipleDetection(
    simulations = radioactivePlumes_area,
    locations = 1:nLocations(radioactivePlumes_area))

  expect_equal(
    multipleDetectionAll[["cost"]], 0)
  expect_equivalent(
    multipleDetectionAll[["costPlumes"]],
    matrix(radioactivePlumes_area@plumes$nDetectable,
           nrow = nPlumes(radioactivePlumes_area), ncol = 1))

  ## some sensors
  multipleDetection1 = multipleDetection(
    simulations = radioactivePlumes_area,
    locations = sampleLocations1)

  expect_equal(
    multipleDetection1[[2]][i],
    sum(getValues(subset(radioactivePlumes_area,
                         locations = sampleLocations1,
                         plumes = i, kinds = "detectable", valuesOnly = TRUE)))
  )
  expect_equal(
    multipleDetection1[[1]],
    mean((radioactivePlumes_area@plumes$nDetectable - multipleDetection1[[2]])/
      radioactivePlumes_area@plumes$nDetectable)
  )


 # early detection
 ## no sensors
 earlyDetection0 = earlyDetection(
   simulations = radioactivePlumes_area,
   locations = integer(0))

 expect_equal(
   earlyDetection0[[1]], 1)
 expect_equivalent(
   earlyDetection0[[2]], matrix(Inf, nrow = nPlumes(radioactivePlumes_area), ncol = 1))# detection time, not cost of single plume

 ## all sensors
 earlyDetectionAll = earlyDetection(
   simulations = radioactivePlumes_area,
   locations = 1:nLocations(radioactivePlumes_area))

 expect_equal(
   earlyDetectionAll[[1]], 0)
 expect_equivalent(
   earlyDetectionAll[[2]], matrix(radioactivePlumes_area@plumes$earliestDetection,
                                  nrow = nPlumes(radioactivePlumes_area), ncol = 1))# detection time, not cost of single plume

 ## some sensors
 earlyDetection1 = earlyDetection(
   simulations = radioactivePlumes_area,
   locations = sampleLocations1)

 expect_equal(
   earlyDetection1[[2]][i],
   min(getValues(subset(radioactivePlumes_area,
                        locations = sampleLocations1, plumes = i,
                        kinds = "time", valuesOnly = TRUE)),
       na.rm = TRUE)
 )
 fun_Rp_detectionTime = function(x, nout = 1, weight = 1,
                   notDetected = 6.0480e+05 * 2){
   xFinite = x - radioactivePlumes_area@plumes$earliestDetection
   xFinite[is.infinite(x)] = notDetected
   out = mean(xFinite)/notDetected
 }
 expect_equal(
   earlyDetection1[[1]],
   fun_Rp_detectionTime(x = earlyDetection1[[2]])
 )
})
