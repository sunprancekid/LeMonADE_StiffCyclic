#include <iostream>
#include <string>
#include <unistd.h> //for getopt
#include <cstring>
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
// #include <clara>
#include <UpdaterCreateRingMelt.h>

int main(int argc, char* argv[])
{
    try{
        // constants
        const uint max_bonds=4;
        int type1(1);
        // establish default parameters
        int chainLength = 100; // number of monomers in single polymers
        int numChains = 1; // number of unique polymers (if rings, interlocking)
        int boxSize = 256; // box length in one dimensions
        bool ring = false; // determines whether a ring should be generated
        bool force = false; // determines whether molecules will experience a force
        double conForce = 0.; // base strength of force
        std::string outfile = "config.bfm"; // output file that contains the configurations
        bool bendingPot = false; // determines whether a bending potential should be added
        double k_theta = 0.; // parameterized bending potential strength

        // determine if any options were passed to the executable
        // read in options by getopt
        int option_char(0);
        while ((option_char = getopt (argc, argv, "n:m:o:rb:k:f:h"))  != EOF){
            switch (option_char)
            {
                // TODO add force oscilation and amplitude
                case 'n':
                    chainLength = atoi(optarg);
                    break;
                case 'm':
                    numChains = atoi(optarg);
                    break;
                case 'o':
                    outfile=optarg;
                    break;
                case 'r':
                    ring = true;
                    break;
                case 'b':
                    boxSize = atoi(optarg);
                    break;
                case 'k':
                    bendingPot = true;
                    k_theta = stod(optarg);
                    break;
                case 'f':
                    force = true;
                    conForce = stod(optarg);
                    break;
                case 'h':
                default:
                    std::cerr << "Usage: " << argv[0] << " \ngen_config << OPTIONS >> \n[-o filenameOutput] \n[-n monomer in ring / chain] \n[-m number of rings / chains] \n[-r generate ring (otherwise generate chain)] \n[-k bending potential strenth (otherwise no bending potential)] \n[-f constant force that molecules experience (otherwise no force is applied)] \n[-b box size]\n";
                    return 0;
            }
        }

        // establish ingredients based on input parameters
        if (bendingPot) {
            // TODO add bending potential package
            std::cout << "TODO add bending potential package to configurations generation routine." << std::endl;
            exit(0);
            // ingredients have four features
            // typedef LOKI_TYPELIST_4(FeatureMoleculesIO, FeatureAttributes< >,FeatureExcludedVolumeSc<>, FeaturePotentialBending) Features;
            // typedef ConfigureSystem<VectorInt3,Features, max_bonds> Config;
            // typedef Ingredients<Config> IngredientsType;
            // IngredientsType ingredients;
        } else {
            // do nothing
        }
        // generate configuration
        // ingredients have three features
        typedef LOKI_TYPELIST_3(
            FeatureMoleculesIO,
            FeatureAttributes<>,
            FeatureExcludedVolumeSc< FeatureLatticePowerOfTwo < > >) Features;
        typedef ConfigureSystem<VectorInt3,Features,max_bonds> Config;
        typedef Ingredients<Config> IngredientsType;
        IngredientsType ingredients;

        // seed randomness
        RandomNumberGenerators rng;
        rng.seedAll();

        // task manager contains updaters and synchronizers
        // ingredients.synchronize(ingredients);
        TaskManager taskManager;
        // build molecule(s) based on input arguments
        if (ring) {
            // build ring molecule
            UpdaterCreateRingMelt<IngredientsType> UCRM(ingredients, numChains, chainLength, boxSize, boxSize, boxSize);
            UCRM.initialize();
            UCRM.execute();
            UCRM.cleanup();
        } else {
            // build chain molecule
            // set box parameters
            ingredients.setBoxX(boxSize);
            ingredients.setBoxY(boxSize);
            ingredients.setBoxZ(boxSize);
            // turn on periodic boundaries
            ingredients.setPeriodicX(true);
            ingredients.setPeriodicY(true);
            ingredients.setPeriodicZ(true);
            // add monomer bonding
            ingredients.modifyBondset().addBFMclassicBondset();
            ingredients.synchronize(ingredients);
            taskManager.addUpdater(new UpdaterAddLinearChains<IngredientsType>(ingredients, numChains,chainLength,type1,type1),0);
            taskManager.initialize();
            taskManager.run();
            taskManager.cleanup();
        }

        // add constant linear force
        // ingredients.setForceOn(true);
        // ingredients.setAmplitudeForce(conForce);

        // create character array for output file name
        // remove the config file if it already exists
        char* outfile_char_array = new char[outfile.length() + 1];
        strcpy(outfile_char_array, outfile.c_str());
        remove(outfile_char_array);

        ingredients.synchronize(ingredients);
        AnalyzerWriteBfmFile<IngredientsType> AWBFM(outfile,ingredients);
        AWBFM.initialize();
        AWBFM.execute();
        AWBFM.cleanup();
    }
    catch(std::exception& err){std::cerr<<err.what();}


    return 0;
}
