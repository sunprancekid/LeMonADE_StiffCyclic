#include <iostream>
#include <string>
using namespace std;

#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
#include <LeMonADE/feature/FeatureLinearForce.h>

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
    int t_equil = atoi(argv[4]);
    // fifth argument :: double for constant force parameter
    double conForce = stod(argv[5]);

    int nChains(1),type1(1);

    typedef LOKI_TYPELIST_4(
        FeatureMoleculesIO,
        FeatureAttributes<>,
        FeatureExcludedVolumeSc< FeatureLatticePowerOfTwo < > >,
        FeatureLinearForce) Features;
    const uint max_bonds=4;
    typedef ConfigureSystem<VectorInt3,Features,max_bonds> Config;
    typedef Ingredients<Config> IngredientsType;
    IngredientsType ingredients;

    RandomNumberGenerators rng;
    rng.seedAll();

    // set box parameters
    ingredients.setBoxX(256);
    ingredients.setBoxY(256);
    ingredients.setBoxZ(256);
    // turn on periodic boundaries
    ingredients.setPeriodicX(true);
    ingredients.setPeriodicY(true);
    ingredients.setPeriodicZ(true);
    // add monomer bonding
    ingredients.modifyBondset().addBFMclassicBondset();
    // add constant linear force
    ingredients.setForceOn(true);
    ingredients.setAmplitudeForce(conForce);
    // ingredients.modifyMolecules().setAge(0);
    // synchronize
    ingredients.synchronize();

    TaskManager taskManager;
    taskManager.addUpdater(new UpdaterAddLinearChains<IngredientsType>(ingredients, nChains,chainLength,type1,type1),0);
    taskManager.addUpdater(new UpdaterSimpleSimulator<IngredientsType,MoveLocalSc>(ingredients,nMCS));
    taskManager.addAnalyzer(new AnalyzerWriteBfmFile<IngredientsType>("config_ev.bfm",ingredients,AnalyzerWriteBfmFile<IngredientsType>::APPEND));
    taskManager.addAnalyzer(new AnalyzerRadiusOfGyration<IngredientsType>(ingredients, "ROG.dat"));
    taskManager.addAnalyzer(new AnalyzerEndToEndDistance<IngredientsType>(ingredients, "RE2E.dat", t_equil)); // TODO :: equilibriation time
    // TODO :: add RouseTimeScale Property Calculations
    // TODO :: add BondBondCorrelation Property Calculations
    // TODO :: add radial distribution accumulation

    taskManager.initialize();
    taskManager.run(nRuns);
    taskManager.cleanup();

    return 0;
}