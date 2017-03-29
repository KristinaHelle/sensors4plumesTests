########################################################
# test            SpatialDataFrame                     #
########################################################
data(SIndexDF)
data(SPointsDF)
data(SPixelsDF)
data(SPolygridDF)
data(SPolygonsDF)
data(SLinesDF)

# coordinates
expect_warning(coordinates(SIndexDF))
coordinates(SPointsDF)
coordinates(SPixelsDF)
coordinates(SPolygridDF)
coordinates(SPolygonsDF)
coordinates(SLinesDF)
# proj4string
expect_warning(proj4string(SIndexDF))
proj4string(SPointsDF)
proj4string(SPixelsDF)
proj4string(SPolygridDF)
proj4string(SPolygonsDF)
proj4string(SLinesDF)
# bbox
expect_warning(bbox(SIndexDF))
bbox(SPointsDF)
bbox(SPixelsDF)
bbox(SPolygridDF)
bbox(SPolygonsDF)
bbox(SLinesDF)
# spplot
expect_warning(spplot(SIndexDF, zcol = 1))
spplot(SPointsDF, zcol = 1)
spplot(SPixelsDF, zcol = 1)
spplot(SPolygridDF, zcol = 1)
spplot(SPolygonsDF, zcol = 1)
spplot(SLinesDF, zcol = 1)

# coordCentral (not exported)
expect_warning(spplot(SIndexDF, 1, sp.layout = list("sp.points", sensors4plumes:::coordCentral(SIndexDF))))
spplot(SPointsDF, 1, sp.layout = list("sp.points", sensors4plumes:::coordCentral(SPointsDF)))
spplot(SPixelsDF, 1, sp.layout = list("sp.points", sensors4plumes:::coordCentral(SPixelsDF)))
spplot(SPolygridDF, 1, sp.layout = list("sp.points", sensors4plumes:::coordCentral(SPolygridDF)))
spplot(SPolygonsDF, 1, sp.layout = list("sp.points", sensors4plumes:::coordCentral(SPolygonsDF)))
spplot(SLinesDF, 1, sp.layout = list("sp.points", sensors4plumes:::coordCentral(SLinesDF)))






