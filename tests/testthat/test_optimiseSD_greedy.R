################################################################
# test optimiseSD_greedy                                       #
################################################################
# deleteSensor
# addSensor
# determineNextStep
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

meanFun = function(x){mean(x, na.rm = TRUE)}
minDist = replaceDefault(
  spatialSpread, newDefaults = list(
    fun = minimalDistance,
    fun_R = meanFun
  ), type = "costFun.optimiseSD"
)[["fun"]]
mean1 = function(x, nout = 1, weight_l = 1, weight_p = 1){
  mean(x[,1], na.rm = TRUE)
}

delineationError_kl1 = replaceDefault(
  interpolationError,
  newDefaults = list(
  kinds = 2,
  fun_interpolation = idw0z,
  fun_error = delineationError,
  fun_Rpl = mean1,
  fun_Rpl_cellStats = "mean",
  fun_l = delineationErrorMap),
  type = "costFun.optimiseSD"
)[["fun"]]
# ---------------------------------- deleteSensor ------------------------ #
del1 = sensors4plumes:::deleteSensor(simulations = SimulationsSmall,
             costFun = minDist,
             locationsDeletable = locDel1,
             locationsKeep = locKeep1)

costWithout1 = numeric(length(locDel1))
for (i in seq(along = locDel1)){
  costWithout1[i] = minDist(SimulationsSmall, c(locDel1[-i], locKeep1))[[1]]
}
expect_equal(
  del1[[3]],
  costWithout1)
expect_equal(
  del1[[2]],
  which(costWithout1 == min(costWithout1))
  )
expect_true(
  is.element(setdiff(locDel1, del1[[1]]), locDel1[which(costWithout1 == min(costWithout1))])
  )

#
del2 = sensors4plumes:::deleteSensor(simulations = radioactivePlumes_local,
                    costFun = minDist,
                    locationsDeletable = locDel2,
                    locationsKeep = locKeep2)

costWithout2 = numeric(length(locDel2))
for (i in seq(along = locDel2)){
  costWithout2[i] = minDist(radioactivePlumes_local, c(locDel2[-i], locKeep2))[[1]]
}
expect_equal(
  del2[[3]],
  costWithout2)
expect_equal(
  del2[[2]],
  which(costWithout2 == min(costWithout2))
)
expect_true(
  is.element(setdiff(locDel2, del2[[1]]), locDel2[which(costWithout2 == min(costWithout2))])
)


# ---------------- addSensor ----------------------------- #
add_delinE_kl1 = sensors4plumes:::addSensor(
  costFun = delineationError_kl1,
  simulations = radioactivePlumes_local,
  locationsAddable = locDel2,        # currently empty potential sensor locations; one of them shall be added
  locationsCurrent = locKeep2         # current sampling design
)

costWith1 = numeric(length(locDel2))
for (i in seq(along = locDel2)){
  costWith1[i] = delineationError_kl1(radioactivePlumes_local, c(locDel2[i], locKeep2))[[1]]
}
expect_equal(
  add_delinE_kl1[[3]],
  costWith1)
expect_equal(
  add_delinE_kl1[[2]],
  which(costWith1 == min(costWith1))
)
expect_true(
  is.element(add_delinE_kl1[[1]], locDel2[which(costWith1 == min(costWith1))])
)

# ------------------------ determineNextStep ------------------------ #
# - - - - -  kindAim = "number" - - - - -
# stop because aim reached (aimNumber = length(locCurrent))
next_stopIt = sensors4plumes:::determineNextStep(
  locationsCurrent = locDel1,
  locationsFix = locKeep1,
  locationsAll = locAll1,
  kindAim = "number",
  aim = 2,
  maxIterations = 10,
  l = 10,
  costFun = minDist,
  simulations = SimulationsSmall
)
expect_equal(next_stopIt[[1]], "finish")

# del, because aim not reached (aimNumber < length(locCurrent))
next_stopIt = sensors4plumes:::determineNextStep(
  locationsCurrent = locDel1,                      # indices of current sensors
  locationsFix = locKeep1,                      # indices of current sensors that are fix
  locationsAll = locAll1,            # indices of all possible sensor locations, including fix ones
  kindAim = "number",                                # if aim is a certain number of sensors ("number") or a certain cost ("cost")
  aim = 1,                                        # desired number of sensors or desired cost
  maxIterations = 10,                       # maximal number of iterations, usually optimisation should stop earlier by other criteria
  l = 10,                                            # iteration step
  costFun = minDist,                                # to be forwarded to computeCost
  simulations = SimulationsSmall
)
expect_equal(next_stopIt[[1]], "del")

# add, because aim not reached (aimNumber > length(locCurrent))
next_stopIt = sensors4plumes:::determineNextStep(
  locationsCurrent = locDel1,                      # indices of current sensors
  locationsFix = locKeep1,                      # indices of current sensors that are fix
  locationsAll = locAll1,            # indices of all possible sensor locations, including fix ones
  kindAim = "number",                                # if aim is a certain number of sensors ("number") or a certain cost ("cost")
  aim = 5,                                        # desired number of sensors or desired cost
  maxIterations = 10,                       # maximal number of iterations, usually optimisation should stop earlier by other criteria
  l = 10,                                            # iteration step
  costFun = minDist,                                # to be forwarded to computeCost
  simulations = SimulationsSmall
)
expect_equal(next_stopIt[[1]], "add")

# - - - - -  kindAim = "cost" - - - - -
# stop because aim reached (aimCost = cost(locCurrent))
next_stopIt = sensors4plumes:::determineNextStep(
  locationsCurrent = locDel1,                      # indices of current sensors
  locationsFix = locKeep1,                      # indices of current sensors that are fix
  locationsAll = locAll1,            # indices of all possible sensor locations, including fix ones
  kindAim = "cost",                                # if aim is a certain number of sensors ("number") or a certain cost ("cost")
  aim = minDist(simulations = SimulationsSmall,
                locations = locDel1)[[1]],                                        # desired number of sensors or desired cost
  maxIterations = 10,                       # maximal number of iterations, usually optimisation should stop earlier by other criteria
  l = 10,                                            # iteration step
  costFun = minDist,                                # to be forwarded to computeCost
  simulations = SimulationsSmall
)
expect_equal(next_stopIt[[1]], "finish")

# del because aim not reached (aimCost > cost(locCurrent))
next_stopIt = sensors4plumes:::determineNextStep(
  locationsCurrent = locDel1,
  locationsFix = locKeep1,
  locationsAll = locAll1,
  kindAim = "cost",
  aim = minDist(simulations = SimulationsSmall,
                locations = locDel1)[[1]] * 2,
  maxIterations = 10,
  l = 10,                                            # iteration step
  costFun = minDist,                                # to be forwarded to computeCost
  simulations = SimulationsSmall
)
expect_equal(next_stopIt[[1]], "del")

# add, because aim not reached (aimCost < cost(locCurrent))
next_stopIt = sensors4plumes:::determineNextStep(
  locationsCurrent = locDel1,                      # indices of current sensors
  locationsFix = locKeep1,                      # indices of current sensors that are fix
  locationsAll = locAll1,            # indices of all possible sensor locations, including fix ones
  kindAim = "cost",                                # if aim is a certain number of sensors ("number") or a certain cost ("cost")
  aim = minDist(simulations = SimulationsSmall,
                locations = locDel1)[[1]] * 1/2,                                        # desired number of sensors or desired cost
  maxIterations = 10,                       # maximal number of iterations, usually optimisation should stop earlier by other criteria
  l = 10,                                            # iteration step
  costFun = minDist,                                # to be forwarded to computeCost
  simulations = SimulationsSmall
)
expect_equal(next_stopIt[[1]], "add")

# stop
## stop because iterations full (l > maxIterations)
expect_warning(
  next_stopIt <- sensors4plumes:::determineNextStep(
    locationsCurrent = locDel1,                      # indices of current sensors
    locationsFix = locKeep1,                      # indices of current sensors that are fix
    locationsAll = locAll1,            # indices of all possible sensor locations, including fix ones
    kindAim = "number",                                # if aim is a certain number of sensors ("number") or a certain cost ("cost")
    aim = 4,                                        # desired number of sensors or desired cost
    maxIterations = 10,                       # maximal number of iterations, usually optimisation should stop earlier by other criteria
    l = 12,                                            # iteration step
    costFun = minDist,                                # to be forwarded to computeCost
    simulations = SimulationsSmall
  )
)
expect_equal(next_stopIt[[1]], "stop")

## stop because all sensors used and cost not reached
expect_warning(
  next_stopIt <- sensors4plumes:::determineNextStep(
    locationsCurrent = locAll1,                      # indices of current sensors
    locationsFix = locKeep1,                      # indices of current sensors that are fix
    locationsAll = locAll1,            # indices of all possible sensor locations, including fix ones
    kindAim = "cost",                                # if aim is a certain number of sensors ("number") or a certain cost ("cost")
    aim = 0,                                        # desired number of sensors or desired cost
    maxIterations = 10,                       # maximal number of iterations, usually optimisation should stop earlier by other criteria
    l = 10,                                            # iteration step
    costFun = minDist,                                # to be forwarded to computeCost
    simulations = SimulationsSmall
  )
)
expect_equal(next_stopIt[[1]], "stop")

## stop because all sensors included and number not reached
expect_warning(
  next_stopIt <- sensors4plumes:::determineNextStep(
    locationsCurrent = 1:9,
    locationsFix = locKeep1,
    locationsAll = locAll1,
    kindAim = "number",
    aim = 10,
    maxIterations = 10,
    l = 10,
    costFun = minDist,
    simulations = SimulationsSmall
  )
)
expect_equal(next_stopIt[[1]], "stop")

## stop because only fix sensors used and cost still below aimCost
expect_warning(
  next_stopIt <- sensors4plumes:::determineNextStep(
    locationsCurrent = locDel1,                      # indices of current sensors
    locationsFix = locDel1,                      # indices of current sensors that are fix
    locationsAll = locAll1,            # indices of all possible sensor locations, including fix ones
    kindAim = "cost",                                # if aim is a certain number of sensors ("number") or a certain cost ("cost")
    aim = minDist(
      simulations = SimulationsSmall,
      locations = locDel1[-1])[[1]],                                        # desired number of sensors or desired cost
    maxIterations = 10,                       # maximal number of iterations, usually optimisation should stop earlier by other criteria
    l = 10,                                            # iteration step
    costFun = minDist,                                # to be forwarded to computeCost
    simulations = SimulationsSmall
  )
)
expect_equal(next_stopIt[[1]], "stop")


## stop because only fix sensors used and number still above aimNumber (case does not occur if fixed sensors are not taken into account)
# case should not occur as locationsFix is a subset of locationsCurrent
# expect_warning(
 #   next_stopIt <- sensors4plumes:::determineNextStep(
 #     locationsCurrent = integer(0),                      # indices of current sensors
 #     locationsFix = locDel1,                      # indices of current sensors that are fix
 #     locationsAll = locAll1,            # indices of all possible sensor locations, including fix ones
 #     kindAim = "number",                                # if aim is a certain number of sensors ("number") or a certain cost ("cost")
 #     aim = 0,                                        # desired number of sensors or desired cost
 #     maxIterations = 10,                       # maximal number of iterations, usually optimisation should stop earlier by other criteria
 #     l = 10,                                            # iteration step
 #     costFun = minDist,                                # to be forwarded to computeCost
 #     simulations = SimulationsSmall
 #   )
 # )
 # expect_equal(next_stopIt[[1]], "stop")


# ---------------------- optimiseSD_greedy ---------------------------- #
# aimNumber (exact), add, no swap
optGreedy1 = optimiseSD_greedy(
  simulations = SimulationsSmall,
  costFun = minDist,
  locationsAll = locAll1,
  locationsFix = locKeep1,
  locationsInitial = locDel1,
  aimNumber = 8,
  maxIterations = 20,
  swap = FALSE
)
i = sample.int(5, 1)
expect_equal(
  optGreedy1$SD[[i]][1,],
  sort(optGreedy1$report$SDs[[i]])
)
expect_equal(# holds only here (only addition)
  optGreedy1$evaluation,
  optGreedy1$report$evalSDs
)
expect_equal(
  optGreedy1$report$evalSDs$number,
  4:8
)
expect_equal(
  optGreedy1$evaluation$cost[i],
  minDist(simulations = SimulationsSmall, locations = optGreedy1$SD[[i]][1,])[[1]]
)
# aimNumber (exact), delete, no swap
optGreedy2 = optimiseSD_greedy(
  simulations = SimulationsSmall,
  costFun = minDist,
  locationsAll = locAll1,
  locationsFix = locKeep1,
  locationsInitial = locDel1,
  aimNumber = 3,
  maxIterations = 20,
  swap = FALSE
)
expect_equal(
  optGreedy2$report$evalSDs$number,
  4:3
)
expect_equal(
  sapply(X = optGreedy2$SD, FUN = dim)[2,],
  3:4
)
expect_equivalent(
  optGreedy2$report$evalSDs[2:1,],
  optGreedy2$evaluation
)


# aimCost (to undershoot), add, no swap
optGreedy3 = optimiseSD_greedy(
  simulations = SimulationsSmall,
  costFun = minDist,
  locationsAll = locAll1,
  locationsFix = locKeep1,
  locationsInitial = locDel1,
  aimCost = 0.17,
  maxIterations = 20,
  swap = FALSE
)
expect_equal(optGreedy1, optGreedy3)

# aimCost (to undershoot), del, no swap
optGreedy4 = optimiseSD_greedy(
  simulations = SimulationsSmall,
  costFun = minDist,
  locationsAll = locAll1,
  locationsFix = locKeep1,
  locationsInitial = locDel1,
  aimCost = 1.3,
  maxIterations = 20,
  swap = FALSE
)
expect_equal(
  optGreedy2,
  optGreedy4
)

# aimNumber(exact), add, swap
optGreedy5 = optimiseSD_greedy(
  simulations = SimulationsSmall,
  costFun = minDist,
  locationsAll = locAll1,
  locationsFix = locKeep1,
  locationsInitial = locDel1,
  aimNumber = 8,
  maxIterations = 20,
  swap = TRUE
)
expect_equal(
  optGreedy5$report$evalSDs$number,
  c(4:8, 7)
)
expect_equal(
  sapply(X = optGreedy5$SD, FUN = dim)[2,],
  4:8
)
expect_equivalent(
  optGreedy5$report$evalSDs[c(1:3,6,5),],
  optGreedy5$evaluation
)

# aimNumber(exact), delete, swap
optGreedy6 = optimiseSD_greedy(
  simulations = SimulationsSmall,
  costFun = minDist,
  locationsAll = locAll1,
  locationsFix = locKeep1,
  locationsInitial = locDel1,
  aimNumber = 4,
  maxIterations = 20,
  swap = TRUE
)
expect_equal(
  optGreedy6$report$evalSDs$number,
  c(4,3,4,3,4,3)
)
expect_equivalent(
  optGreedy6$report$evalSDs[6:5,],
  optGreedy6$evaluation
)
expect_equivalent(
  sort(optGreedy6$report$SDs[[6]]),
  optGreedy6$SD[[1]][1,]
)
expect_equivalent(
  sort(optGreedy6$report$SDs[[5]]),
  optGreedy6$SD[[2]][1,]
)
# aimCost (to undershoot), add, swap
optGreedy7 = optimiseSD_greedy(
  simulations = SimulationsSmall,
  costFun = minDist,
  locationsAll = locAll1,
  locationsFix = locKeep1,
  locationsInitial = locDel1,
  aimCost = 0.17,
  maxIterations = 20,
  swap = TRUE
)
expect_equal(
  optGreedy7$report$evalSDs$number,
  c(4:8,7)
)
# aimCost (to undershoot), del, swap
optGreedy8 = optimiseSD_greedy(
  simulations = SimulationsSmall,
  costFun = minDist,
  locationsAll = locAll1,
  locationsFix = locKeep1,
  locationsInitial = locDel1,
  aimCost = 1.3,
  maxIterations = 20,
  swap = TRUE
)
expect_equal(
  optGreedy8$report$evalSDs$number,
  c(4,3,4,3,4,3,4)
)
expect_equal(# SDs of relevant cost: all different except the last ones
  nrow(unique(matrix(unlist(optGreedy8$report$SDs[optGreedy8$report$evalSDs$cost <= 1.3]), byrow = TRUE, ncol = 4))),
  3
)
# various reasons to stop (with/without swap)
# cost, add, maxIterations, no swap
optGreedy9 = optimiseSD_greedy(
  simulations = SimulationsSmall,
  costFun = minDist,
  locationsAll = locAll1,
  locationsFix = locKeep1,
  locationsInitial = locDel1,
  aimCost = 0.17,
  maxIterations = 2,
  swap = FALSE
)
expect_equal(
  optGreedy9$report$evalSDs$number,
  c(4:6)
)
# expect_equal(
#   optGreedy9$finalSDwhich,
#   integer(0)
# )
# cost, del, not reachable, no swap
optGreedy10 = optimiseSD_greedy(
  simulations = SimulationsSmall,
  costFun = minDist,
  locationsAll = locAll1,
  locationsFix = locKeep1,
  locationsInitial = locDel1,
  aimCost = 0,
  maxIterations = 20,
  swap = FALSE
)
expect_equal(
  optGreedy10$report$evalSDs$number,
  c(4:8)
)
# expect_equal(
#   optGreedy10$finalSDwhich,
#   integer(0)
# )
# number, add, not reachable, no swap
optGreedy11 = optimiseSD_greedy(
  simulations = SimulationsSmall,
  costFun = minDist,
  locationsAll = locAll1,
  locationsFix = locKeep1,
  locationsInitial = locDel1,
  aimNumber = 9,
  maxIterations = 20,
  swap = FALSE
)
expect_equal(
  optGreedy11$report$evalSDs$number,
  c(4:8)
)
# number, del, maxIterations, swap
expect_error( # number of fix sensors above aimNumber
  optGreedy12 = optimiseSD_greedy(
  simulations = SimulationsSmall,
  costFun = minDist,
  locationsAll = locAll1,
  locationsFix = locKeep1,
  locationsInitial = locDel1,
  aimNumber = 2,
  maxIterations = 1,
  swap = TRUE
)
)

# cost, del, not reachable, no swap
optGreedy13 = optimiseSD_greedy(
  simulations = SimulationsSmall,
  costFun = minDist,
  locationsAll = locAll1,
  locationsFix = locKeep1,
  locationsInitial = locDel1,
  aimCost = 4,
  maxIterations = 20,
  swap = FALSE
)
expect_equal(
  optGreedy13$report$evalSDs$number,
  c(4:2)
)
# number, add, not reachable, swap
optGreedy14 = optimiseSD_greedy(
  simulations = SimulationsSmall,
  costFun = minDist,
  locationsAll = locAll1,
  locationsFix = locKeep1,
  locationsInitial = locDel1,
  aimNumber = 9,
  maxIterations = 20,
  swap = TRUE
)
expect_equal(
  optGreedy14$report$evalSDs$number,
  c(4:8)
)

# number (ignored) and cost, swap
optGreedy15 = optimiseSD_greedy(
  simulations = SimulationsSmall,
  costFun = minDist,
  locationsAll = locAll1,
  locationsFix = locKeep1,
  locationsInitial = locDel1,
  aimCost = 1.3,
  aimNumber = 9,
  maxIterations = 20,
  swap = TRUE
)
expect_equal(
  optGreedy15,
  optGreedy8
)
# big dataset
optGreedy_l1 = optimiseSD_greedy(
  simulations = radioactivePlumes_local,
  costFun = minDist,
  locationsAll = locAll2,
  locationsFix = locKeep2,
  locationsInitial = locDel2,
  aimNumber = 101,
  maxIterations = 20,
  swap = TRUE
)
optGreedy_l2 = optimiseSD_greedy(
  simulations = radioactivePlumes_local,
  costFun = minDist,
  locationsAll = locAll2,
  locationsFix = locKeep2,
  locationsInitial = locDel2,
  aimNumber = 109,
  maxIterations = 20,
  swap = TRUE
)
# saving
# file.remove("data/optGr8.Rdata")
# filesBefore = list.files("data")
# optGreedy8_ = optimiseSD_greedy(
#   simulations = SimulationsSmall,
#   costFun = minDist,
#   locationsAll = locAll1,
#   locationsFix = locKeep1,
#   locationsInitial = locDel1,
#   aimCost = 1.3,
#   maxIterations = 20,
#   swap = TRUE,
#   nameSave = "data/optGr8"
# )
# filesAfter = list.files("data")
# setdiff(filesAfter, filesBefore) == "optGr8.Rdata"
# load("data/optGr8.Rdata")
# expect_equal(
#   SDs,
#   optGreedy8_$SDs
# )
# file.remove("data/optGr8.Rdata")

# plotting: was parameter to be forwarded to cost function each time it is called; else no functionality

# optimiseSD_greedy called by optimiseSD
optGreedy8 = optimiseSD_greedy(
  simulations = SimulationsSmall,
  costFun = minDist,
  locationsAll = locAll1,
  locationsFix = locKeep1,
  locationsInitial = locDel1,
  aimCost = 1.3,
  maxIterations = 20,
  swap = TRUE
)

optimSD_greedy = replaceDefault(
  optimiseSD_greedy,  newDefaults = list(
    maxIterations = 20,
    swap = TRUE
    ),
  type = "optimisationFun.optimiseSD")[[1]]

opt8 = optimiseSD(
  simulations = SimulationsSmall,
  costFun = minDist,
  locationsAll = locAll1,
  locationsFix = locKeep1,
  locationsInitial = locDel1,
  aimCost = 1.3,
  optimisationFun = optimSD_greedy,
  nameSave = paste0(path1,"/test1/optimiseSD8")
)
expect_equal(opt8, optGreedy8)

optimSD_greedy2 = replaceDefault(
  optimiseSD_greedy,  newDefaults = list(
    maxIterations = 20,
    swap = FALSE
  ),
  type = "optimisationFun.optimiseSD")[[1]]

opt4 = optimiseSD(
  simulations = SimulationsSmall,
  costFun = minDist,
  locationsAll = locAll1,
  locationsFix = locKeep1,
  locationsInitial = locDel1,
  aimCost = 1.3,
  optimisationFun = optimSD_greedy2,
  nameSave = paste0(path1,"/test1/optimiseSD8")
)
expect_equal(opt4, optGreedy4)
