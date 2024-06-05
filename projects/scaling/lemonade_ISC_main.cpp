#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureAttributes.h>

#include <LeMonADE/updater/UpdaterAddLinearChains.h>
#include <LeMonADE/updater/UpdaterSimpleSimulator.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
#include <LeMonADE/analyzer/AnalyzerRadiusOfGyration.h>

#include <LeMonADE/utility/RandomNumberGenerators.h>
#include <LeMonADE/utility/TaskManager.h>

// modules not found in LeMonADE Library
#include <AnalyzerEndToEndDistance.h>

int main(int argc, char* argv[])
{

  //  first argument :: integer for chain length
  int chainLength = atoi(argv[1]);
  // second argument :: integer for number of monte carlo steps
  int nMCS = atoi(argv[2]);
  // third argument :: integer for number of runs
  int nRuns = atoi(argv[3]);
  // fourth argument :: double for equilibriation time
  int t_equil = atof(argv[4]);

  int nChains(1),type1(1);
  
  typedef LOKI_TYPELIST_2(
    FeatureMoleculesIO, 
    FeatureAttributes<>) Features;
  const uint max_bonds=4;
  typedef ConfigureSystem<VectorInt3,Features,max_bonds> Config;
  typedef Ingredients<Config> IngredientsType;
  IngredientsType ingredients;
    
  RandomNumberGenerators rng;
  rng.seedAll();
  
  ingredients.setBoxX(chainLength);
  ingredients.setBoxY(chainLength);
  ingredients.setBoxZ(chainLength);
  ingredients.setPeriodicX(true);
  ingredients.setPeriodicY(true);
  ingredients.setPeriodicZ(true);
  ingredients.modifyBondset().addBFMclassicBondset();
  ingredients.synchronize();

  TaskManager taskManager;
  taskManager.addUpdater(new UpdaterAddLinearChains<IngredientsType>(ingredients, nChains,chainLength,type1,type1),0);
  taskManager.addUpdater(new UpdaterSimpleSimulator<IngredientsType,MoveLocalSc>(ingredients,nMCS));
  taskManager.addAnalyzer(new AnalyzerWriteBfmFile<IngredientsType>("config.bfm",ingredients,AnalyzerWriteBfmFile<IngredientsType>::APPEND));
  taskManager.addAnalyzer(new AnalyzerRadiusOfGyration<IngredientsType>(ingredients, "ROG.dat"));
  taskManager.addAnalyzer(new AnalyzerEndToEndDistance<IngredientsType>(ingredients, "ROG2.dat"));
  // TODO :: add RouseTimeScale Property Calculations
  // TODO :: add BondBondCorrelation Property Calculations
  
  taskManager.initialize();
  taskManager.run(nRuns);
  taskManager.cleanup();
  
  return 0;
}
