# End-To-End Distance Analyzer
add_library(AnalyzerEndToEndDistance MODULE AnalyzerEndToEndDistance.h)
# utility for calculating property statists
add_library(StatisticMoment MODULE StatisticMoment.h)
# Analyzer for calculating properties during constant force simulations
add_library(AnalyzerConstantForce MODULE AnalyzerConstantForce.h)
# utility for accumulating histogram statistics
add_library(HistogramGeneralStatistik1D MODULE HistogramGeneralStatistik1D.h)
# Analyzer for calculating bond bond angular distribution during simulations
add_library(AnalyzerBondBondDistribution MODULE AnalyzerBondBondDistribution.h)
