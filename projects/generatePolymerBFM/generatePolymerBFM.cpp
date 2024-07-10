#include <iostream>
#include <string>
#include <unistd.h> //for getopt
#include <cstring>
using namespace std;

#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>

#include <LeMonADE/updater/UpdaterAddLinearChains.h>
#include <LeMonADE/updater/UpdaterSimpleSimulator.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
#include <LeMonADE/analyzer/AnalyzerRadiusOfGyration.h>

#include <LeMonADE/utility/RandomNumberGenerators.h>
#include <LeMonADE/utility/TaskManager.h>
#include <LeMonADE/utility/Vector3D.h>

// modules not found in LeMonADE Library
// #include <clara>
#include <UpdaterCreateRingMelt.h>
#include <FeaturePotentialBending.h>
#include <FeatureOscillatoryForce.h>

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
        std::string outfile = "config_init.bfm"; // output file that contains the configurations
        bool bendingPot = false; // determines whether a bending potential should be added
        double k_theta = 0.; // parameterized bending potential strength
        bool bendingPot_CA = false; // determines if a CA potential should be used for bending potential interactions (default is CSA)
        std::string forceVecString = "111";

        // determine if any options were passed to the executable
        // read in options by getopt
        int option_char(0);
        while ((option_char = getopt (argc, argv, "n:m:o:rb:k:f:v:ch"))  != EOF){
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
                case 'v':
                    forceVecString = optarg;
                    break;
                case 'c':
                    bendingPot_CA = true;
                    break;
                case 'h':
                default:
                    std::cerr << "\n\nUsage: ./generatePolymerBFM << OPTIONS >> \n[-o filenameOutput] \n[-n number of monomers in ring / chain] \n[-m number of rings / chains] \n[-r generate ring (otherwise generate chain)] \n[-k bending potential strength (otherwise no bending potential)] \n[-f constant force that molecules experience (otherwise no force is applied)] \n[-v string with three integers xyz denoting the force orientation in each dimension (default is " << forceVecString << ")] \n[-b box size]\n[-c use cosine angle potential for bending potential (default is cosine square angle potential)]\n\n";
                    return 0;
            }
        }

        // generate ingredients
        typedef LOKI_TYPELIST_5(FeatureMoleculesIO,
                                FeatureAttributes< >,
                                FeatureExcludedVolumeSc<>,
                                FeaturePotentialBending,
                                FeatureOscillatoryForce) Features;
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

        // assign bending potential to polymers, if specified
        if (bendingPot) {
            // TODO assumes only one macromolecule, adjust algoritm to consider multiple
            // add the bending information for ring chain (middle monomers affect 3 angles)
            for(uint32_t i=2;i<chainLength-2;i++){
                ingredients.modifyMolecules()[i].setBendingBondInformation(std::make_pair(i-2,i-1) , std::make_pair(i-1,i  ));
                ingredients.modifyMolecules()[i].setBendingBondInformation(std::make_pair(i-1,i  ) , std::make_pair(i  ,i+1));
                ingredients.modifyMolecules()[i].setBendingBondInformation(std::make_pair(i  ,i+1) , std::make_pair(i+1,i+2));
            }
            if (ring) {
                // if the macromolecule is a ring, the first and last monomers are bonded
                // add the bending information for ring chain (first and last monomers affect 3 angles)
                size_t monomerIdxEnd = chainLength-1;
                ingredients.modifyMolecules()[0].setBendingBondInformation(std::make_pair(0,1) , std::make_pair(1,2));
                ingredients.modifyMolecules()[0].setBendingBondInformation(std::make_pair(monomerIdxEnd,0) , std::make_pair(0,1));
                ingredients.modifyMolecules()[0].setBendingBondInformation(std::make_pair(monomerIdxEnd-1,monomerIdxEnd) , std::make_pair(monomerIdxEnd,0));

                ingredients.modifyMolecules()[monomerIdxEnd].setBendingBondInformation(std::make_pair(monomerIdxEnd-2,monomerIdxEnd-1) , std::make_pair(monomerIdxEnd-1,monomerIdxEnd));
                ingredients.modifyMolecules()[monomerIdxEnd].setBendingBondInformation(std::make_pair(monomerIdxEnd-1,monomerIdxEnd) , std::make_pair(monomerIdxEnd,0));
                ingredients.modifyMolecules()[monomerIdxEnd].setBendingBondInformation(std::make_pair(monomerIdxEnd,0) , std::make_pair(0,1));

                // add the bending information for ring chain (second and before last monomers affect 2 angles)
                size_t monomerIdx = 1;
                ingredients.modifyMolecules()[monomerIdx].setBendingBondInformation(std::make_pair(monomerIdxEnd,0)  , std::make_pair(0, monomerIdx));
                ingredients.modifyMolecules()[monomerIdx].setBendingBondInformation(std::make_pair(monomerIdx-1,monomerIdx )  , std::make_pair(monomerIdx  ,monomerIdx+1));
                ingredients.modifyMolecules()[monomerIdx].setBendingBondInformation(std::make_pair(monomerIdx  ,monomerIdx+1) , std::make_pair(monomerIdx+1,monomerIdx+2));

                monomerIdx = chainLength-2;
                ingredients.modifyMolecules()[monomerIdx].setBendingBondInformation(std::make_pair(monomerIdx-2,monomerIdx-1) , std::make_pair(monomerIdx-1,monomerIdx));
                ingredients.modifyMolecules()[monomerIdx].setBendingBondInformation(std::make_pair(monomerIdx-1,monomerIdx )  , std::make_pair(monomerIdx  ,monomerIdx+1));
                ingredients.modifyMolecules()[monomerIdx].setBendingBondInformation(std::make_pair(monomerIdx  ,monomerIdx+1)  , std::make_pair(monomerIdx+1 ,0));
            } else {
                // if the macro molecule is a chain, the first and second monomers in the chain are partially bonded to the rest of the chain
                // the beginning of the chain
                ingredients.modifyMolecules()[0].setBendingBondInformation(std::make_pair(0,1), std::make_pair(1,2));
                ingredients.modifyMolecules()[1].setBendingBondInformation(std::make_pair(0,1), std::make_pair(1,2));
                ingredients.modifyMolecules()[1].setBendingBondInformation(std::make_pair(1,2), std::make_pair(2,3));
                // the end of the chain
                size_t monomerIdx = chainLength - 1;
                ingredients.modifyMolecules()[monomerIdx].setBendingBondInformation(std::make_pair(monomerIdx-2, monomerIdx-1), std::make_pair(monomerIdx-1, monomerIdx));
                monomerIdx = chainLength - 2;
                ingredients.modifyMolecules()[monomerIdx].setBendingBondInformation(std::make_pair(monomerIdx-2, monomerIdx-1), std::make_pair(monomerIdx-1, monomerIdx));
                ingredients.modifyMolecules()[monomerIdx].setBendingBondInformation(std::make_pair(monomerIdx-1, monomerIdx), std::make_pair(monomerIdx, monomerIdx + 1));
            }
            // set bending constant
            ingredients.setBending_Potential_Constant(k_theta);
            if (bendingPot_CA) {
                ingredients.setBending_Potential_Type_CA();
            } else {
                ingredients.setBending_Potential_Type_CSA();
            }
            // synchronize
            ingredients.synchronize(ingredients);
        }

        // add force
        if (force) {
            // add constant linear force
            ingredients.setForceOn(true);
            ingredients.setBaseForce(conForce);
            // parse the force vector
            VectorDouble3 fv;
            for(std::string::size_type i = 0; i < forceVecString.size(); ++i) {
                std::string s(1, forceVecString[i]);
                if (i == 0) {
                    // x component
                    fv.setX(stod(s));
                } else if (i == 1) {
                    // y component
                    fv.setY(stod(s));
                } else if (i == 2) {
                    // z component
                    fv.setZ(stod(s));

                }
            }
            ingredients.setForceVector(fv);
            // synchronize
            ingredients.synchronize(ingredients);
        }

        // assign the molecules that experience forces
        if (ring) {
            for (uint32_t i=0; i < ingredients.getMolecules().size();i++){
                if ((i % chainLength) == 0){
                    // the first molecule in the ring experiences a positive force
                    ingredients.modifyMolecules()[i].setAttributeTag(4);
                } else if ((i % chainLength) == floor(ingredients.getMolecules().size() / 2.)) {
                    // the middle molecule in the ring experiences a negative force
                    ingredients.modifyMolecules()[i].setAttributeTag(5);
                } else {
                    // otherwise, set the attribute tag to 0
                    ingredients.modifyMolecules()[i].setAttributeTag(0);
                }
            }
        } else {
            for (uint32_t i=0; i < ingredients.getMolecules().size();i++){
                if ((i % chainLength) == 0){
                    // the first molecule in the chain experiences a positive force
                    ingredients.modifyMolecules()[i].setAttributeTag(4);
                } else if ((i % chainLength) == (ingredients.getMolecules().size() - 1)) {
                    // the last molecule in the chain experiences a negative force
                    ingredients.modifyMolecules()[i].setAttributeTag(5);
                } else {
                    // otherwise, set the attribute tag to 0
                    ingredients.modifyMolecules()[i].setAttributeTag(0);
                }
            }
        }

        // create character array for output file name
        // remove the output file if it already exists
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
