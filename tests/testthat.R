library(testthat)
library(sensors4plumes)
library(sensors4plumesData)

source("tests/testthat/test_SpatialPolygridDataFrame_prepare.R") # ok
source("tests/testthat/test_SpatialPolygridDataFrame_SpatialPolygridDataFrame.R") # ok
source("tests/testthat/test_SpatialPolygridDataFrame_coerce.R") # ok
source("tests/testthat/test_points2polygrid.R") # ok, some tests disabled as expect_error seems not to do what I expect
source("tests/testthat/test_polygrid2grid.R") # ok
source("tests/testthat/test_SpatialPolygridDataFrame_coordinates_proj4string_bbox.R") # ok
source("tests/testthat/test_SpatialPolygridDataFrame_subsetSDF.R") # ok
source("tests/testthat/test_SpatialPolygridDataFrame_cbind.R") # ok
source("tests/testthat/test_SpatialPolygridDataFrame_spplot.R") # ok
source("tests/testthat/test_SpatialPolygridDataFrame_testData.R") # ok
source("tests/testthat/test_testDataArtificial.R") # ok
source("tests/testthat/test_subsetSDF.R") # ok
source("tests/testthat/test_areaSDF.R") # ok area on SPolygons does work differently
source("tests/testthat/test_SpatialIndexDataFrame.R") # ok, tests removed of deleted function rbind.SpatialIndexDataFrame
source("tests/testthat/test_SpatialDataFrame.R") # ok
source("tests/testthat/test_testDataTelescopic.R") # ok
source("tests/testthat/test_spplotLog.R") # ok,does only make sense 'by hand'
source("tests/testthat/test_Simulations.R") # ok
source("tests/testthat/test_cbindSimulations.R") # ok
source("tests/testthat/test_extractSpatialDataFrame.R") # ok
source("tests/testthat/test_SDF2simulations.R") # ok
path1 = "/home/kristina/sensors4plumesTests" # path to save files from tests
source("tests/testthat/test_loadSimulations.R") # ok
source("tests/testthat/test_changeSimulationsPath.R") # ok
source("tests/testthat/test_copySimulations.R") # ok
source("tests/testthat/test_radioactivePlumes.R") # ok
source("tests/testthat/test_summaryPlumes.R") # ok
source("tests/testthat/test_summaryLocations.R") # ok
source("tests/testthat/test_subsetSimulations.R") # ok
source("tests/testthat/test_replaceDefault.R") # ok (to be updated if all types have been defined (?))
source("tests/testthat/test_simulationsApply.R") # ok (don't run like this: files generated...)
source("tests/testthat/test_interpolationError.R") # ok
source("tests/testthat/test_interpolation.R") # ok (don't run: generates files, takes long)
source("tests/testthat/test_measurementResult.R") # ok (does not test 'measurementResult' but the derived functions)
source("tests/testthat/test_costFunctions.R") # ok
source("tests/testthat/test_optimiseSD_ssa.R") # ok
source("tests/testthat/test_optimiseSD_genetic.R") # ok
source("tests/testthat/test_optimiseSD_greedy.R") # ok
source("tests/testthat/test_optimiseSD_global.R") # ok
source("tests/testthat/test_optimiseSD_manual.R") # ok (except plotting/clicking for SIndexDF, SLinesDF) interactive, not to be run automatically
source("tests/testthat/test_plotSD.R") # ok (do by hand and visually check result)
source("tests/testthat/test_optimisationCurve.R") # ok (do by hand and visually check result)
source("tests/testthat/test_postprocessing.R") # ok (do by hand and visually check result)
