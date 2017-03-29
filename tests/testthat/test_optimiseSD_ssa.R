################################################################
# test optimiseSD_ssa                                          #
################################################################

data(SimulationsSmall)
locDel1 = sample.int(nLocations(SimulationsSmall), 2)
locKeep1 = sample(setdiff(1:nLocations(SimulationsSmall), locDel1), 2)
locAll1 = c(sample(setdiff(1:nLocations(SimulationsSmall), c(locDel1, locKeep1)), 4), locDel1)
data(radioactivePlumes_local)
locDel2 = sample.int(nLocations(radioactivePlumes_local), 5)
locKeep2 = sample(setdiff(1:nLocations(radioactivePlumes_local), locDel2), 100)
locAll2 = c(sample(setdiff(1:nLocations(radioactivePlumes_local), c(locDel2, locKeep2)), 10), locDel2)
data(radioactivePlumes_area)
locDel3 = sample.int(nLocations(radioactivePlumes_area), 5)
locKeep3 = sample(setdiff(1:nLocations(radioactivePlumes_area), locDel3), 100)
locAll3 = c(sample(setdiff(1:nLocations(radioactivePlumes_area), c(locDel3, locKeep3)), 10), locDel3)

threshold = 1e-7
radioactivePlumes_area@values$detectable = calc(
  radioactivePlumes_area@values$maxdose,
  fun = function(x){x >= threshold})
radioactivePlumes_area@plumes$totalDose =
  summaryPlumes(radioactivePlumes_area, fun = sum, kinds = "finaldose")[[2]]
radioactivePlumes_area@plumes$nDetectable =
  summaryPlumes(radioactivePlumes_area, fun = sum, kinds = "detectable")[[2]]
radioactivePlumes_area@plumes$earliestDetection =
  summaryPlumes(radioactivePlumes_area, fun = min, kinds = "time", na.rm = TRUE)[[2]]

meanFun = function(x){mean(x, na.rm = TRUE)}
minDist = replaceDefault(
  spatialSpread, newDefaults = list(
    fun = minimalDistance,
    fun_R = meanFun
  ), type = "costFun.optimiseSD"
)[["fun"]]
costInitial_minDist_local = minDist(simulations = radioactivePlumes_local,
                                                      locations = c(locKeep2, locDel2))
costInitial_multipleDetection_area = multipleDetection(simulations = radioactivePlumes_area,
                                                      locations = c(locKeep3, locDel3))

# ----------- test subfunctions --------------------------------#

# - - - - -- - -  adjustShift - - -- - - - - -- - -- - - - - - #
# input: all locations
## random
allLocations = subsetSDF(radioactivePlumes_local@locations,
                         locations = setdiff(1:nLocations(radioactivePlumes_local), locDel2))
## small in x-direction, big in y-direction
locationsAll3 = c(1:441, 5850:5949)
allLocations2 = subsetSDF(radioactivePlumes_local@locations,
                         locations = locationsAll3)

# run function
adjShift1 = sensors4plumes:::adjustShift(simulations = radioactivePlumes_local,
                        locationsAll = setdiff(1:nLocations(radioactivePlumes_local), locDel2),
                        maxShiftFactor = 0.2,
                        minShiftFactor = 0.01)

expect_warning(
  adjShift2 <- sensors4plumes:::adjustShift(simulations = radioactivePlumes_local,
                           locationsAll = setdiff(1:nLocations(radioactivePlumes_local), locDel2),
                           maxShiftFactor = 0.2,
                           minShiftFactor = 0.001) # too small, changed -> warning
)

expect_warning(
  adjShift3 <- sensors4plumes:::adjustShift(simulations = radioactivePlumes_local,
                           locationsAll = locationsAll3, # smaller total area -> minShift too small -> warning
                           maxShiftFactor = 0.2,
                           minShiftFactor = 0.01)
)

# 1) no correction needed
# expect shift to be shiftFactor * extent
maxShift = apply(FUN = diff, X = bbox(allLocations), MARGIN = 1)
minShift = radioactivePlumes_local@locations@grid@cellsize/maxShift
expect_equivalent(adjShift1[["shift"]][,"max"],
                  maxShift * 0.2)
expect_equivalent(adjShift1[["shift"]][,"min"],
                  maxShift * 0.01)

expect_true(0.01 >= max(minShift)) # min shift factor not too small
expect_gt(min(adjShift1[["shift"]]), min(radioactivePlumes_local@locations@grid@cellsize))
expect_equal(adjShift1[["allLocations"]], allLocations)

# 2) minShiftFactor too small, minimal shift set to cellsize
expect_true(0.001 <= max(minShift)) # min shift factor too small
expect_equal(min(adjShift2[["shift"]]), # change to meet smallest change = cellsize
             min(radioactivePlumes_local@locations@grid@cellsize))

# 3) smaller extent -> minShift too small; minimal shift set to cellsize
maxShift2 = apply(FUN = diff, X = bbox(allLocations2), MARGIN = 1)
minShift2 = radioactivePlumes_local@locations@grid@cellsize/maxShift2
expect_lt(0.01, max(minShift2)) # min shift factor too small
expect_equal(min(adjShift3[["shift"]]), # change to meet smallest change = cellsize
             min(radioactivePlumes_local@locations@grid@cellsize))
## only shift in y direction is set to cellsize; in x direction minShift is big enough:
expect_equivalent(adjShift3[["shift"]]["y","min"],
             radioactivePlumes_local@locations@grid@cellsize["y_sim"])# set to cellsize
expect_equal(adjShift3[["shift"]]["x",], adjShift1[["shift"]]["x",])# not affected, area has same extent

# 4) minShift too small, non-square cells
expect_warning(
  adjShift1S <- sensors4plumes:::adjustShift(simulations = SimulationsSmall,
                            maxShiftFactor = 0.2,
                            minShiftFactor = 0.1)
)
## shift in each direction adjusted to respective cell size
expect_equivalent(adjShift1S[["shift"]],
                  matrix(c(1,1.5,1,1.5), nrow = 2))

## other
expect_warning(
  adjShift4 <- sensors4plumes:::adjustShift(simulations = radioactivePlumes_area,
                                            locationsAll = setdiff(1:nLocations(radioactivePlumes_local), locDel3), # smaller total area -> minShift too small -> warning
                                            maxShiftFactor = 1,
                                            minShiftFactor = 0.01)
)



# - - - - -- - -  - - estimateCostDiff  - - -- - - - - -- - #
# average cost difference in first (random move) and final (neighbour move) iterations
#  [always starting from initalLocations, moving 1 sensor]
# is not taking into account shift! nor irregular order!
set.seed(18121965)
costDiff1 = sensors4plumes:::estimateCostDiff(
  simulations = radioactivePlumes_local,
  locationsInitial = locDel2,
  locationsFix = locKeep2,
  costFun = minDist,
  costInitial = costInitial_minDist_local[[1]])

set.seed(18121965)
del = integer(10)
repl = list()
costDiff = numeric(10)
# delete random locations
for (i in 1:10){
  del[i] = sample.int(length(locDel2), 1)
}
# replace by random (other, not fix) location [first sample half]
for (i in 1:5){
  repl[[i]] = c(setdiff(locDel2, locDel2[del[i]]),
                sample.int(nLocations(radioactivePlumes_local), 1),
                locKeep2)
}
# replace by neighbour (in terms of index) [second sample half]
for (i in 6:10){
  if (locDel2[del[i]] <= nLocations(radioactivePlumes_local)/2){
    repl[[i]] = c(setdiff(locDel2, locDel2[del[i]]),
                  locDel2[del[i]] + 1,
                  locKeep2)
  } else {
    repl[[i]] = c(setdiff(locDel2, locDel2[del[i]]),
                  locDel2[del[i]] - 1,
                  locKeep2)
  }
}
# compute cost difference
for (i in 1:10){
  costDiff[i] = abs(costInitial_minDist_local[[1]] - minDist(simulations = radioactivePlumes_local,
                                               locations = repl[[i]])[[1]])
}
# expect result to be median of first and second sample half respectively
expect_equal(costDiff1[1],
             median(costDiff[1:5]))
expect_equal(costDiff1[2],
             median(costDiff[6:10]))

## test with different cost function and simulations
costDiff2 = sensors4plumes:::estimateCostDiff(
  simulations = radioactivePlumes_area,
  locationsInitial = locDel3,
  locationsFix = locKeep3,
  costFun = multipleDetection,
  costInitial = multipleDetection(simulations = radioactivePlumes_area,
                                locations = c(locDel3, locKeep3))[[1]])

# - - - - -- - - - test moveSensors  -- - - - - - -- -- - - -- - -- - #
# coordinates(subsetSDF(simulations@locations, locations = locationsPrev[which_shifts[i]]))[1]) <=
## as in optimiseSD_ssa: add original indices to selected coordinates
allLocations = adjShift1[["allLocations"]]
allLocations@data$indexOrig = setdiff(1:nLocations(radioactivePlumes_local), locDel2)

maxShiftNumber_i = 6
maxIterations_i = 5000
k_i = 2000
startMoveProb_i = 1
set.seed(04121942)
moveSensors1 = sensors4plumes:::moveSensors(
  allLocations = allLocations,
  locationsPrev = locAll2[1:10],
  maxShiftNumber = maxShiftNumber_i,
  maxIterations = maxIterations_i,
  k = k_i,
  startMoveProb = startMoveProb_i,
  shift = adjShift1[["shift"]]
)

# shift decreases linearly
shift = adjShift1[["shift"]][,"max"] - apply(FUN = diff, X = adjShift1[["shift"]], MARGIN = 1) * k_i/maxIterations_i

# potential neighbours are all within coordinate limits;
#                 that are part of allLocations;
#                 and not part of the current locations (locAll2)
coordAll = coordinates(radioactivePlumes_local@locations)
neighbours = list()
for (i in 1:10){
  coordOld = coordAll[locAll2[i],]
  neighbours[[i]] = setdiff(
    intersect(
      which(abs(coordAll[,1] - coordOld[1]) < shift[1] & abs(coordAll[,2] - coordOld[2]) < shift[2]),
      allLocations@data$indexOrig),
    locAll2)
}

set.seed(04121942)
# randomly replace maxShiftNumber_i locations
change = data.frame(del = sample.int(10, maxShiftNumber_i),
                    repl = integer(maxShiftNumber_i))
# replace them only with limited probability (the first one is replaced in any case)
ifChange = runif(maxShiftNumber_i) < startMoveProb_i * (1 - k_i/maxIterations_i)
ifChange[1] = TRUE
change = change[ifChange,]
# replace by value from neighbourhood
for (i in seq(along = change$repl)){
  change$repl[i] = sample(neighbours[[change$del[i]]], 1)
}
newLoc = c(locAll2[1:10][-change$del], change$repl)
expect_equal(
  sort(newLoc),
  sort(moveSensors1[[1]])
)

# visualisation showing neighbourhood etc. (neighbourhoods are squares, replacement is in neighbourhood)
plot(coordAll, pch = ".")
points(coordAll[unlist(neighbours),], pch = ".", col = 3)
points(coordAll[locAll2[change$del],], pch = as.character(1:sum(ifChange)), col = 2)
points(coordAll[change$repl,], pch =  as.character(1:sum(ifChange)), col = 4)

## other simulations
allLocations3 = adjShift4[["allLocations"]]
allLocations3@data$indexOrig = setdiff(1:nLocations(radioactivePlumes_area), locDel3)

moveSensors2 = sensors4plumes:::moveSensors(
  allLocations = allLocations3,
  locationsPrev = locAll3[1:10],
  maxShiftNumber = 10,
  maxIterations = 10000,
  k = 1000,
  startMoveProb = 1,
  shift = adjShift4[["shift"]]
)


# - - - - - --  - - - test formulas of cooling parameters --- - - #
#    - - - - - -  van Groenigen
## cooling and startAcceptance are calibrated correctly
start_acc_vG = runif(1)
end_acc_vG = runif(1) * 0.001
maxIterations = sample.int(5000, 1) + 2500
cooling = ((log(start_acc_vG) * costDiff1[2]) /
             (log(end_acc_vG) * costDiff1[1])) ^ {1 / (maxIterations - 1)}
startAcceptance = - costDiff1[1]/
  (cooling * log(start_acc_vG))
chiBegin = exp(-costDiff1[1] / (startAcceptance * cooling ^ 1)) # formula from l 299; only used if costDiff < 0
chiEnd = exp(-costDiff1[2] / (startAcceptance * cooling ^ maxIterations))

expect_equal(chiBegin, start_acc_vG)
expect_equal(chiEnd, end_acc_vG)

## effect of all parameters in chi
chiExamples = data.frame(costDiff = rep(c(1, 0.01), times = 8),
                         k = rep(c(1,1,1000,1000), times = 4),
                         startAcceptance = rep(c(0.8, 0.2, 0.8, 0.2), each = 4),
                         cooling = rep(c(0.999, 0.9999), each = 8),
                         chi = numeric(16))
for (i in 1:16){
  chiExamples$chi[i] = exp(-chiExamples$costDiff[i] /
                          (chiExamples$startAcceptance[i] * chiExamples$cooling[i] ^ chiExamples$k[i]))
}
expect_true(# bigger cost difference -> smaller chi
  all(chiExamples$chi[c(1,3,5,7,9,11,13,15)] <= chiExamples$chi[c(2,4,6,8,10,12,14,16)])
)
expect_true(# earlier iteration -> bigger chi
  all(chiExamples$chi[c(1,2,5,6,9,10,13,14)] >= chiExamples$chi[c(3,4,7,8,11,12,15,16)])
)
expect_true(# bigger startAcceptance -> bigger chi
  all(chiExamples$chi[c(1:4, 9:12)] >= chiExamples$chi[c(5:8, 13:16)])
)
expect_true(# moore cooling (cooling value farther from 1) -> chi decreases faster
  all(chiExamples$chi[c(1,2,5,6)] / chiExamples$chi[c(3,4,7,8)] >=
      chiExamples$chi[c(9,10,13,14)] / chiExamples$chi[c(11,12,15,16)])
)

## chi always in [0,1]?
# chi = exp(costDiff / (startAcceptance * cooling ^ k))
# costDiff < 0; startAcceptance > 0, cooling > 0; => chi in (0,1]

#   -  - - -- -- - - intamapInteractive
## is it the same? yes:
# accept worse if [intamapInteractive: ssaMap.R l153] runif(1) <= (start_p*exp(-k/(nr_iterations/10)))

## effect of parameters
chiExamples_iI = data.frame(start_p = rep(c(0.1, 0.9), times = 2),
                            k = rep(c(1, 1000), each = 2),
                            chi = numeric(4))
for (i in 1:4){
  chiExamples_iI$chi[i] = chiExamples_iI$start_p[i] * exp(- chiExamples_iI$k[i] / (1000/10) )
}
expect_true(# higher start acceptance -> higher chi
  all(chiExamples_iI$chi[c(1,3)] <= chiExamples_iI$chi[c(2,4)])
)
expect_true(# later iteration -> lower chi
  all(chiExamples_iI$chi[1:2] >= chiExamples_iI$chi[3:4])
)

## chi always in [0,1]
# start_p in [0,1]; k > 0; nr_iterations > 0 => chi in [0,1]



# ---------- test complete function: fix, to test effect of certain parameters ----- #
# - - - - - - iteration - - - - - - - -
# is a single iteration processed correctly? (looks correct from plot)
maxIt = 20
set.seed(10112004)
expect_warning( # minShiftFactor set to cellsize
  optSSA1 = optimiseSD_ssa(
    simulations = radioactivePlumes_area,
    costFun = multipleDetection,
    locationsAll = setdiff(1:nLocations(radioactivePlumes_area), locKeep3),
    locationsFix = locKeep3,
    locationsInitial = locDel3,
    aimCost = 0.05 * costInitial_multipleDetection_area[[1]],
    verbatim = TRUE,
    maxIterations = maxIt
  )
)

#save(optSSA1, file = "optSSA1.Rdata")

# iteration types (1st iteration excluded as type unknown)
costDecreases_optSSA1 = which((c(optSSA1$report$iterations$costTest, 0) <=
                                 c(0, optSSA1$report$iterations$cost))[1:maxIt]) # 1st value FALSE: in fact: unknown
costIncreasesButAccepted_optSSA1 = setdiff(which(optSSA1$report$iterations$accepted), c(costDecreases_optSSA1, 1))
notAccepted_optSSA1 = setdiff(which(!optSSA1$report$iterations$accepted), 1)
costMin_optSSA1 = which(optSSA1$report$iterations$costTest <= cummin(optSSA1$report$iterations$costTest))

expect_equal(# cost is only minimal when it decreases
  setdiff(costMin_optSSA1, costDecreases_optSSA1), 1
)
expect_equal(
  sort(c(costDecreases_optSSA1,
         costIncreasesButAccepted_optSSA1,
         notAccepted_optSSA1)),
  2:maxIt
)
expect_true(# cost increases -> chi computed
  all(optSSA1$report$iterations$chi[- c(1, costDecreases_optSSA1)] > 0)
)

# cost decreases
expect_true(# -> accepted
  all(optSSA1$report$iterations$accepted[costDecreases_optSSA1])
)
expect_equal(# -> costTest becomes cost
  optSSA1$report$iterations$cost[costDecreases_optSSA1],
  optSSA1$report$iterations$costTest[costDecreases_optSSA1]
)
expect_equal(# -> SDtest becomes SD
  optSSA1$report$SDs[costDecreases_optSSA1,],
  optSSA1$report$SDs_test[costDecreases_optSSA1,]
)
expect_equal(# -> costTest becomes costBest
  optSSA1$report$iterations$costBest[costMin_optSSA1],
  optSSA1$report$iterations$costTest[costMin_optSSA1]
)
expect_equal(# -> SD becomes SDbest
  optSSA1$report$SDs_best[costMin_optSSA1,],
  optSSA1$report$SDs_test[costMin_optSSA1,]
)

# cost increases but accepted
expect_equal(# -> costTest becomes cost
  optSSA1$report$iterations$cost[costIncreasesButAccepted_optSSA1],
  optSSA1$report$iterations$costTest[costIncreasesButAccepted_optSSA1]
)
expect_equal(# -> SDtest becomes SD
  optSSA1$report$SDs[costIncreasesButAccepted_optSSA1,],
  optSSA1$report$SDs_test[costIncreasesButAccepted_optSSA1,]
)
expect_equal(# -> costBest is kept
  optSSA1$report$iterations$costBest[costIncreasesButAccepted_optSSA1],
  optSSA1$report$iterations$costBest[costIncreasesButAccepted_optSSA1 - 1]
)
expect_equal(# -> SDbest is kept
  optSSA1$report$SDs_best[costIncreasesButAccepted_optSSA1,],
  optSSA1$report$SDs_best[costIncreasesButAccepted_optSSA1 - 1,]
)

# not accepted
expect_equal(# -> cost is kept
  optSSA1$report$iterations$cost[notAccepted_optSSA1],
  optSSA1$report$iterations$cost[notAccepted_optSSA1 - 1]
)
expect_equal(# -> SD is kept
  optSSA1$report$SDs[notAccepted_optSSA1,],
  optSSA1$report$SDs[notAccepted_optSSA1 - 1,]
)
expect_equal(# -> costBest is kept
  optSSA1$report$iterations$costBest[notAccepted_optSSA1],
  optSSA1$report$iterations$costBest[notAccepted_optSSA1 - 1]
)
expect_equal(# -> SDbest is kept
  optSSA1$report$SDs_best[notAccepted_optSSA1,],
  optSSA1$report$SDs_best[notAccepted_optSSA1 - 1,]
)

# cost fits SDs
costSDs_test = numeric(maxIt)
for (i in 1:maxIt){
  costSDs_test[i] = multipleDetection(simulations = radioactivePlumes_area,
                                      locations = optSSA1$report$SDs_test[i,])[[1]]
}
expect_equal(
  optSSA1$report$iterations$costTest,
  costSDs_test
)
costSDs = numeric(maxIt)
for (i in 1:maxIt){
  costSDs[i] = multipleDetection(simulations = radioactivePlumes_area,
                                      locations = optSSA1$report$SDs[i,])[[1]]
}
expect_equal(
  optSSA1$report$iterations$cost,
  costSDs
)
costSDs_best = numeric(maxIt)
for (i in 1:maxIt){
  costSDs_best[i] = multipleDetection(simulations = radioactivePlumes_area,
                                 locations = optSSA1$report$SDs_best[i,])[[1]]
}
expect_equal(
  optSSA1$report$iterations$costBest,
  costSDs_best
)

# - - - - - - total  - - - - - - - -
## global best
expect_equal(# global best is lowest costTest
  min(optSSA1$report$iterations$costTest),
  optSSA1$report$cost_best
)
expect_true(# global best decreases monotonously
  all((c(1, optSSA1$report$iterations$costBest) - c(optSSA1$report$iterations$costBest, 0))[2:maxIt] >= 0)
)
expect_equal(# global best SD is SD of global best)
  optSSA1$report$SD_best,
  optSSA1$report$SDs_test[which.min(optSSA1$report$iterations$costTest),]
)
expect_equal(
  multipleDetection(simulations = radioactivePlumes_area,
                    locations = optSSA1$report$SD_best)[[1]],
  optSSA1$report$cost_best
)
## final
expect_equal(# final SD is last accepted SD
  optSSA1$report$SD,
  optSSA1$SDs[max(which(optSSA1$report$iterations$accepted)),]
)
expect_equal(
  multipleDetection(simulations = radioactivePlumes_area,
                    locations = optSSA1$SD[[1]][1,])[[1]],
  optSSA1$evaluation$cost[1]
)

# -  -- - -- -stop and jump back
# stopping when aim reached
set.seed(20122008)
optSSA2 = optimiseSD_ssa(
  simulations = radioactivePlumes_area,
  costFun = multipleDetection,
  locationsAll = setdiff(1:nLocations(radioactivePlumes_area), locKeep3),
  locationsFix = locKeep3,
  locationsInitial = locDel3,
  start_acc_vG = 0.01,
  aimCost = 0.99 * costInitial_multipleDetection_area[[1]],
  verbatim = TRUE,
  maxIterations = 5000
)

expect_lt(# last cost is below aim
  optSSA2$evaluation$cost,
#  optSSA2$report$iterations$cost[nrow(optSSA2$report$iterations)],
  0.99 * costInitial_multipleDetection_area[[1]]
)
expect_gt(# earlier cost is above aim
  min(optSSA2$report$iterations$cost[-nrow(optSSA2$report$iterations)]),
  0.99 * costInitial_multipleDetection_area[[1]]
)

# stopping when k = maxIterations
maxIt = 1000
set.seed(29061928)
optSSA3 = optimiseSD_ssa(
  simulations = radioactivePlumes_area,
  costFun = multipleDetection,
  locationsAll = setdiff(1:nLocations(radioactivePlumes_area), locKeep3),
  locationsFix = locKeep3,
  locationsInitial = locDel3,
  start_acc_vG = 0.01,
  aimCost = 0.05 * costInitial_multipleDetection_area[[1]],
  verbatim = TRUE,
  maxIterations = maxIt
)
expect_equal(# number of (reported) iterations equals maxIterations
  nrow(optSSA3$report$iterations),
  maxIt
)
plot(optSSA3$report$iterations$costTest, pch = 20, cex = 0.2)
lines(optSSA3$report$iterations$costBest, col = 4)
lines(optSSA3$report$iterations$cost, col = 3, pch = 20, cex = 0.2)
points(optSSA3$report$iterations$costTest, col = optSSA3$report$iterations$accepted + 2, pch = 20, cex = 0.2)


# stopping if no improvement too long
maxStabIt = 200
set.seed(21121993)
optSSA4 = optimiseSD_ssa(
  simulations = radioactivePlumes_area,
  costFun = multipleDetection,
  locationsAll = setdiff(1:nLocations(radioactivePlumes_area), locKeep3),
  locationsFix = locKeep3,
  locationsInitial = locDel3,
  start_acc_vG = 0.01,
  aimCost = 0.05 * costInitial_multipleDetection_area[[1]],
  verbatim = TRUE,
  maxIterations = 5000,
  maxStableIterations = maxStabIt,
  maxIterationsJumpBack = 200
)
costRepoptSSA4 = rle(optSSA4$report$iterations$cost)
expect_equal(# last cost repeated maxStableIterations times
  rev(costRepoptSSA4$lengths)[1],
  maxStabIt + 1
)
expect_lt(# earlier, cost repeated less often
  max(rev(costRepoptSSA4$lengths)[-1]),
  maxStabIt + 1
)

# jumping back if no improvement too long (if jump back desired)
maxItJumpBack = 200
set.seed(15012005)
optSSA5 = optimiseSD_ssa(
  simulations = radioactivePlumes_area,
  costFun = multipleDetection,
  locationsAll = setdiff(1:nLocations(radioactivePlumes_area), c(locKeep3, locAll3)),
  locationsFix = locKeep3,
  locationsInitial = locDel3,
  start_acc_vG = 0.1,
  aimCost = 0.05 * costInitial_multipleDetection_area[[1]],
  verbatim = TRUE,
  maxIterations = 3000,
  maxStableIterations = 500,
  maxIterationsJumpBack = maxItJumpBack
)
bestRepoptSSA5 = rle(optSSA5$report$iterations$costBest)

# jump back when best cost has not improved for maxIterationsJumpBack
whichJBi = which(bestRepoptSSA5$lengths >= maxItJumpBack)
whichJB = integer(0)
for (i in seq(along = whichJBi)){
  k = 1
  while (k * maxItJumpBack <= bestRepoptSSA5$lengths[whichJBi[i]]){
    whichJB = c(whichJB, sum(bestRepoptSSA5$lengths[1:whichJBi[i]-1]) + k * maxItJumpBack + 1)
    k = k + 1
  }
}
expect_equal(
  whichJB,
  optSSA5$report$jumpBack
)
# when jumping back, SD is replaced by last best SD (and cost accordingly)
for (i in seq(along = optSSA5$report$jumpBack)){
  expect_equal(
    optSSA5$report$SDs[optSSA5$report$jumpBack[i], 1:length(locDel3)],
    optSSA5$report$SDs_best[optSSA5$report$jumpBack[i], 1:length(locDel3)]
  )
}

# ---- sensor selection
# fix, all... used correctly (derive from report)
## fix is always part of SD
i = sample.int(nrow(optSSA5$report$SDs), 1)
expect_true(
  all(is.element(locKeep3, optSSA5$report$SDs[i,]))
)
## not in all is never part of SD
expect_false(
  any(is.element(setdiff(locAll3, locDel3), unique(as.integer(optSSA5$report$SDs))))
)

#- saving
# file.remove("data/optSSA6.Rdata")
# filesBefore = list.files("data")
# optSSA6 = optimiseSD_ssa(
#   simulations = radioactivePlumes_area,
#   costFun = multipleDetection,
#   locationsAll = setdiff(1:nLocations(radioactivePlumes_area), locAll3),
#   locationsFix = locKeep3,
#   locationsInitial = locDel3,
#   start_acc_vG = 0.1,
#   aimCost = 0.05 * costInitial_multipleDetection_area[[1]],
#   verbatim = TRUE,
#   maxIterations = 1000,
#   maxStableIterations = 500,
#   maxIterationsJumpBack = 250,
#   nameSave = "data/optSSA6"
# )
# filesAfter = list.files("data")
# expect_true(
#   setdiff(filesAfter, filesBefore) == "optSSA6.Rdata"
# )
# load("data/optSSA6.Rdata")
# expect_equal(
#    SDs,
#    optSSA6$SDs
# )
# file.remove("data/optSSA6.Rdata")

# - plotting (disabled)

# ---------- call via optimiseSD -----------------
optimSD_ssa = replaceDefault(
  optimiseSD_ssa,  newDefaults = list(
    start_acc_vG = 0.1,
    aimCost = 0,
    verbatim = TRUE,
    maxIterations = 3000,
    maxStableIterations = 500,
    maxIterationsJumpBack = 200
  ),
  type = "optimisationFun.optimiseSD")[[1]]

set.seed(15012005)
optSSA7 = optimiseSD(
  simulations = radioactivePlumes_area,
  costFun = multipleDetection,
  locationsAll = setdiff(1:nLocations(radioactivePlumes_area), c(locKeep3, locAll3)),
  locationsFix = locKeep3,
  locationsInitial = locDel3,
  aimCost = 0.05 * costInitial_multipleDetection_area[[1]],
  optimisationFun = optimSD_ssa
)
# equal except verbatim
expect_equal(optSSA5[1:2], optSSA7[1:2])
expect_equal(optSSA5[[3]][1:4], optSSA7[[3]][1:4])
expect_equal(optSSA5[[3]][[5]][,c(2,5)], optSSA7[[3]][[5]])


# other simulations and cost function
# cost functions from test_costFunction
optSSA8 = optimiseSD_ssa(
  simulations = radioactivePlumes_local,
  costFun = delineationError_kl1,
  locationsAll = setdiff(1:nLocations(radioactivePlumes_local), locAll2),
  locationsFix = locKeep2,
  locationsInitial = locDel2,
  start_acc_vG = 0.1,
  aimCost = 0,
  verbatim = TRUE,
  maxIterations = 5,
  maxStableIterations = 5,
  nameSave = "data/optSSA8"
)
set.seed(07021916)
set.seed(377)
optSSA9 = optimiseSD_ssa(
  simulations = SimulationsSmall,
  costFun = absError_i1,
  locationsAll = setdiff(1:nLocations(SimulationsSmall), locAll1),
  locationsFix = locKeep1,
  locationsInitial = locDel1[1:2],
  start_acc_vG = 0.1,
  aimCost = 0,
  verbatim = FALSE,
  maxIterations = 50,
  maxStableIterations = 30,
  nameSave = paste0(path1, "/test1/optSSA9")
)
# plot(optSSA9$report$costTest, pch = ".", col = optSSA9$report$accepted + 2, main = "optSSA9")
# lines(optSSA9$report$costBest, col = 4)
# lines(optSSA9$report$cost, col = 3)
# points(optSSA9$report$costTest, pch = ".", col = optSSA9$report$accepted + 2)

# start from random locations
expect_error(# aimNumber below number of fix sensors
  optimiseSD_ssa(
    simulations = radioactivePlumes_area,
    costFun = multipleDetection,
    locationsAll = setdiff(1:nLocations(radioactivePlumes_area), c(locKeep3, locAll3)),
    locationsFix = locKeep3,
    start_acc_vG = 0.1,
    aimCost = 0.05 * costInitial_multipleDetection_area[[1]],
    aimNumber = 8,
    verbatim = TRUE,
    maxIterations = 30,
    maxStableIterations = 5,
    maxIterationsJumpBack = 6
  )
)

radioactivePlumes_area@values$detectable = calc(
    radioactivePlumes_area@values$maxdose,
   fun = function(x){x >= threshold})
radioactivePlumes_area@plumes$nDetectable =
     summaryPlumes(radioactivePlumes_area, fun = sum, kinds = "detectable")[["summaryPlumes"]]
set.seed(15011915) # Elisabeth & Margarethe Muehlbayer *
optSSA10 = optimiseSD_ssa(
  simulations = radioactivePlumes_area,
  costFun = multipleDetection,
  locationsAll = setdiff(1:nLocations(radioactivePlumes_area), c(locKeep3, locAll3)),
  locationsFix = locKeep3,
  start_acc_vG = 0.1,
  aimCost = 0.05 * costInitial_multipleDetection_area[[1]],
  aimNumber = length(locKeep3) + 8,
  verbatim = TRUE,
  maxIterations = 30,
  maxStableIterations = 5,
  maxIterationsJumpBack = 6
)
expect_equal(
  length(setdiff(optSSA10$SD[[1]][1,], locKeep3)),
  8
)

set.seed(15011915)
optSSA11 = optimiseSD_ssa(
  simulations = radioactivePlumes_area,
  costFun = multipleDetection,
  locationsAll = setdiff(1:nLocations(radioactivePlumes_area), c(locKeep3, locAll3)),
  locationsFix = locKeep3,
  locationsInitial = sample(setdiff(setdiff(1:nLocations(radioactivePlumes_area), c(locKeep3, locAll3)), locKeep3), 8),
  start_acc_vG = 0.1,
  aimCost = 0.05 * costInitial_multipleDetection_area[[1]],
  aimNumber = length(locKeep3) + 7,
  verbatim = TRUE,
  maxIterations = 30,
  maxStableIterations = 5,
  maxIterationsJumpBack = 6
)
expect_equal(
  optSSA10,
  optSSA11
)

