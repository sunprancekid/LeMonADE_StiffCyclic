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
#include <AnalyzerConstantForce.h>
#include <AnalyzerEndToEndDistance.h>
#include <AnalyzerBondBondDistribution.h>

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

    int nChains(1),type1(1), bb_bins(30);

    typedef LOKI_TYPELIST_3(
        FeatureMoleculesIO,
        FeatureAttributes<>,
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
    // taskManager.addAnalyzer(new AnalyzerRadiusOfGyration<IngredientsType>(ingredients, "ROG.dat"));
    taskManager.addAnalyzer(new AnalyzerEndToEndDistance<IngredientsType>(ingredients, "RE2E.dat", t_equil));
    taskManager.addAnalyzer(new AnalyzerBondBondDistribution<IngredientsType>(ingredients, "BBD.dat", t_equil, bb_bins));
    // TODO :: add RouseTimeScale Property Calculations
    // TODO :: add BondBondCorrelation Property Calculations

    taskManager.initialize();

    // assign the molecules that experience forces
    for (uint32_t i=0; i < ingredients.getMolecules().size();i++){
        if (i == 0){
            // the first molecule in the chain experiences a positive force
            ingredients.modifyMolecules()[i].setAttributeTag(4);
        } else if (i == (ingredients.getMolecules().size() - 1)) {
            // the last molecule in the chain experiences a negative force
            ingredients.modifyMolecules()[i].setAttributeTag(5);
        } else {
            // otherwise, set the attribute tag to 0
            ingredients.modifyMolecules()[i].setAttributeTag(0);
        }
        // cout << ingredients.getMolecules()[i].getAttributeTag() << "\n";
    }
    // run simulation
    taskManager.run(nRuns);
    taskManager.cleanup();

    return 0;
}
