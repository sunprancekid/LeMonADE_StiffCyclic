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

  //  read in the chain length from the executable command line
  int n_monomers = atoi(argv[1]);

  int nChains(1),chainLength(n_monomers),type1(1),nMCS(100),nRuns(10);
  
  typedef LOKI_TYPELIST_2(
    FeatureMoleculesIO, 
    FeatureAttributes<>) Features;
  const uint max_bonds=4;
  typedef ConfigureSystem<VectorInt3,Features,max_bonds> Config;
  typedef Ingredients<Config> IngredientsType;
  IngredientsType ingredients;
    
  RandomNumberGenerators rng;
  rng.seedAll();
  
  ingredients.setBoxX(n_monomers);
  ingredients.setBoxY(n_monomers);
  ingredients.setBoxZ(n_monomers);
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
