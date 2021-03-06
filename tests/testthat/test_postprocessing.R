##############################################
#   test    postprocessing                   #
##############################################

## prepare data
# simulations with SpatialPolygridDataFrame
data(radioactivePlumes_local)

# simulations with SpatialPixelsDataFrame
data(meuse.grid)
coordinates(meuse.grid) = ~ x + y
gridded(meuse.grid) = TRUE
meuseGSim = subset(radioactivePlumes_local, 
                  locations = length(meuse.grid), plumes = 1:10)
meuseGSim@locations = meuse.grid

# simulations with rectangular SpatialPixelsDataFrame
data(SPixelsDF)
rectSim = subset(radioactivePlumes_local, 
                 locations = 1:length(SPixelsDF), plumes = 1:10)
rectSim@locations = SPixelsDF

# simulations with SpatialPointsDataFrame
data(meuse)
coordinates(meuse) = ~ x + y
meuseSim = subset(radioactivePlumes_local, 
                  locations = length(meuse), plumes = 1:10)
meuseSim@locations = meuse

# SD: vector, matrix, list
# for meuseGSim, radioactivePlumes_local
SDs0 = sample.int(2001, 10)
SDs1 = matrix(sample.int(2001, 40), nrow = 8)
SDs2 = list(
  sample.int(2001, 5),
  sample.int(2001, 10),
  sample.int(2001, 15),
  sample.int(2001, 20)
)
# for meuseSim
sds0 = c(10,  25,  42,  84,  90,  92,  94,  97, 120, 153)#sample.int(155, 10)
sds1 = matrix(c(97,   79,   40,   68,  131,
                57,   18,   38,  118,   14,
                23,   71,   22,   94,   27,
                125,  108,    4,   80,  129,
                130,   96,  101,  137,   19,
                77,  138,   32,   95,   88,
                140,   73,   43,  153,   12,
                8,  141,   92,   35,  102), 
              byrow = TRUE, nrow = 8)#matrix(sample.int(155, 40), nrow = 8)
sds2 = list(
  c(28, 59, 64, 78, 81),
  c(5,  13,  21,  31,  45,  91,  92, 122, 130, 141),
  c(1,   2,  19,  36,  50,  51,  58,  59,  90, 103, 105, 107, 112, 123, 132),
  c(34,  48,  49,  50,  71,  76,  77,  86,  92,  97, 100, 103, 104, 106, 108, 110, 113, 127, 134, 142)
)

#------------ test returnLonLat --------------------
test_that("SDLonLat", {
  LL0 = SDLonLat(simulations = radioactivePlumes_area, 
                 SD = SDs0)
  LL1 = SDLonLat(simulations = radioactivePlumes_area, 
                 SD = SDs1)
  LL2 = SDLonLat(simulations = radioactivePlumes_area, 
                 SD = SDs2)
  LL = spTransform(
    SpatialPoints(coordinates(radioactivePlumes_area@locations),
                  CRS(proj4string(radioactivePlumes_area@locations))),
    CRS("+init=epsg:4326"))

  expect_equal(
    LL0[[1]],
    LL[SDs0,]
    )
  i = sample.int(nrow(SDs1),1)
  expect_equal(
    LL1[[i]],
    LL[SDs1[i,],]
  )  
  i = sample.int(length(SDs2),1)
  expect_equal(
    LL2[[i]],
    LL[SDs2[[i]],]
  )
})

# ------------ test similaritySD -------------------
test_that("similaritySD", {
  # "equal"
  sim_1_1 = similaritySD(# no simulations needed
    SD = SDs1, referenceSD = SDs0, type = "equal")  
  i = sample.int(nrow(SDs1), 1)
  expect_equal(
    sim_1_1[i],
    sum(is.element(SDs1[i,], SDs0))
  )
  sim_1_2 = similaritySD(# no simulations needed
    SD = SDs2, referenceSD = SDs0, type = "equal")  
  i = sample.int(length(SDs2),1)
  expect_equal(
    sim_1_2[i],
    sum(is.element(SDs2[[i]], SDs0))
  )

  # "Kneighbours
  ## SPointsDF
#   plot(coordinates(meuseSim@locations), pch = ".")
#   points(coordinates(meuseSim@locations)[sds0,], cex = 0.5)
#   for (i in seq(along = sds2)){
#     points(coordinates(meuseSim@locations)[sds2[[i]],], 
#            col = i + 1)#, pch = as.character(i))  
#   }
  sim_2_1_1 = similaritySD(simulations = meuseGSim,
                              SD = sds2, referenceSD = sds0, 
                              type = "Kneighbours", k = 4)  
  expect_equal(sim_2_1_1, c(0,4,2,8))
  sim_2_1_2 = similaritySD(simulations = meuseGSim,
                         SD = sds2, referenceSD = sds0, 
                         type = "Kneighbours", k = 9) 
  expect_equal(sim_2_1_2, c(1,4,6,12))# bigger neighbourhood -> more neighbours

  ## SPixelsDF
#   plot(coordinates(meuseGSim@locations), pch = ".")
#   points(coordinates(meuseGSim@locations)[SDs0,], cex = 0.5)
#   for (i in seq(along = SDs2)){
#     points(coordinates(meuseGSim@locations)[SDs2[[i]],], col = i + 1)#, pch = as.character(i))  
#   }
  sim_2_2_1 = similaritySD(simulations = meuseGSim,
                              SD = SDs2, referenceSD = SDs0, 
                              type = "Kneighbours", k = 9)

  sim_2_2_2 = similaritySD(simulations = meuseGSim,
                              SD = SDs2, referenceSD = SDs0, 
                              type = "Kneighbours", k = 25)   
  
#   plot(coordinates(radioactivePlumes_local@locations)[1:2001,], pch = ".")
#   
#   points(coordinates(radioactivePlumes_local@locations)[SDs0,], cex = 0.5, pch = 2)
#   for (i in 1:nrow(SDs1)){
#     points(coordinates(radioactivePlumes_local@locations)[SDs1[i,],], 
#            col = i, cex = 0.5)#, pch = as.character(i))  
#   }  
  sim_2_3_1 = similaritySD(simulations = radioactivePlumes_local,
                              SD = SDs1, referenceSD = SDs0, 
                              type = "Kneighbours", k = 13)  
  
  # "EarthMoversDistance"
  ## SpatialPixelsDataFrame
  expect_error(similaritySD(simulations = meuseSim,
                                 SD = SDs1, referenceSD = SDs0, 
                                 type = "EarthMoversDistance") 
  )
  sim_3_1_1 = similaritySD(simulations = rectSim,
                                SD = c(2), referenceSD = c(5), 
                                type = "EarthMoversDistance")
  expect_equal(
    sim_3_1_1, 1)
  sim_3_1_2 = similaritySD(simulations = rectSim,
                                SD = c(2,9), referenceSD = c(5,6), 
                                type = "EarthMoversDistance")
  expect_equal(
    sim_3_1_2, 1) # average distance  
  sim_3_1_3 = similaritySD(simulations = rectSim,
                                SD = c(2,9), referenceSD = c(3,8), 
                                type = "EarthMoversDistance")
  expect_equal(
    sim_3_1_3, 1.5) # rectangular grid taken into account
  sim_3_1_4 = similaritySD(simulations = rectSim,
                                SD = c(3,8), referenceSD = c(4,5), 
                                type = "EarthMoversDistance")
  sim_3_1_5 = similaritySD(simulations = rectSim,
                                SD = c(3,8), referenceSD = c(4,5), 
                                type = "EarthMoversDistance", dist = "manhattan")  
  expect_lt(sim_3_1_4, sim_3_1_5) # forwarding parametes to emd works
  
  ## SpatialPolygridDataFrame
  sim_3_1_5 = similaritySD(simulations = radioactivePlumes_local,
                                SD = c(220, 503), referenceSD = c(224, 544), 
                                type = "EarthMoversDistance")  
  expect_equal(
    sim_3_1_5, 1500)
})




# compare SDs
# spatial similarity
# goodness on plumes for calibration and test
