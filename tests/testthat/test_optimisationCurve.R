###############################################################
# test postprocessing                                         #
###############################################################
# use results form optimiseSD of all different algorithms

# greedy
data(SDgreedy)
curve_greedy1 = optimisationCurve(
  optSD = SDgreedy,
  type = "greedy")
curve_greedy2 = optimisationCurve(
  optSD = SDgreedy,
  type = "greedy",
  nameSave = "curve_greedy1.png",
  width = 300, height = 400, bg = "yellow")

# ssa
data(SDssa)
curve_ssa1 = optimisationCurve(
  optSD = SDssa,
  type = "ssa")
curve_ssa2 = optimisationCurve(
  optSD = SDssa,
  type = "ssa",
  nameSave = "curve_ssa1.png")

# report = optSD_gen1$report
data(SDgenetic)
curve_genetic1 = optimisationCurve(
  optSD = SDgenetic,
  type = "genetic")
curve_genetic2 = optimisationCurve(
  optSD = SDgenetic,
  type = "genetic",
  nameSave = "curve_genetic1.png")

# report = optSD_global2
data(SDglobal)
curve_global1 = optimisationCurve(
  optSD = SDglobal,
  type = "global")
curve_global2 = optimisationCurve(
  optSD = SDglobal,
  type = "global",
  nameSave = "curve_global1.png")

# report = optSD_manual2
data(SDmanual)
curve_manual1 = optimisationCurve(
  optSD = SDmanual,
  type = "manual")
curve_manual2 = optimisationCurve(
  optSD = SDmanual,
  type = "manual",
  nameSave = "curve_manual1.png")
