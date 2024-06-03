#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureAttributes.h>

#include <LeMonADE/updater/UpdaterAddLinearChains.h>
#include <LeMonADE/updater/UpdaterSimpleSimulator.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>

#include <LeMonADE/utility/RandomNumberGenerators.h>
#include <LeMonADE/utility/TaskManager.h>

int main(int argc, char* argv[])
{
  int nChains(1),chainLength(64),type1(1),nMCS(100),nRuns(10);
  
  typedef LOKI_TYPELIST_2(
    FeatureMoleculesIO, 
    FeatureAttributes<>) Features;
  const uint max_bonds=4;
  typedef ConfigureSystem<VectorInt3,Features,max_bonds> Config;
  typedef Ingredients<Config> IngredientsType;
  IngredientsType ingredients;
    
  RandomNumberGenerators rng;
  rng.seedAll();
  
  ingredients.setBoxX(64);
  ingredients.setBoxY(64);
  ingredients.setBoxZ(64);
  ingredients.setPeriodicX(true);
  ingredients.setPeriodicY(true);
  ingredients.setPeriodicZ(true);
  ingredients.modifyBondset().addBFMclassicBondset();
  ingredients.synchronize();

  TaskManager taskManager;
  taskManager.addUpdater(new UpdaterAddLinearChains<IngredientsType>(ingredients, nChains,chainLength,type1,type1),0);
  taskManager.addUpdater(new UpdaterSimpleSimulator<IngredientsType,MoveLocalSc>(ingredients,nMCS));
  taskManager.addAnalyzer(new AnalyzerWriteBfmFile<IngredientsType>("config.bfm",ingredients,AnalyzerWriteBfmFile<IngredientsType>::APPEND));
  
  taskManager.initialize();
  taskManager.run(nRuns);
  taskManager.cleanup();
  
  return 0;
}
