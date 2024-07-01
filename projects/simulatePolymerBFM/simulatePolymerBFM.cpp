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
#include <LeMonADE/updater/UpdaterReadBfmFile.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
#include <LeMonADE/analyzer/AnalyzerRadiusOfGyration.h>

#include <LeMonADE/utility/RandomNumberGenerators.h>
#include <LeMonADE/utility/TaskManager.h>

// modules not found in LeMonADE Library
// #include <clara>
#include <FeaturePotentialBending.h>
#include <AnalyzerEndToEndDistance.h>
#include <AnalyzerBondBondCorrelation.h>
#include <AnalyzerBondBondDistribution.h>

int main(int argc, char* argv[])
{
    try{
        // constants
        const uint max_bonds=4;
        int type1(1);
        // establish default parameters
        double save_interval = 1000000; // frequency of property calculations, position save
        double max_MCs = 100000000; // total number of Monte Carlo steps
        std::string infile = "config_init.bfm"; // file that contains initial configuraiton for bfm simulation
        std::string outfile = "config.bfm"; // file that contains system configuratun at each save interval
        double t_equil = 10000000; // simulation time before which properties are not calculated
        // double t_equil = 0;

        // determine if any options were passed to the executable
        // read in options by getopt
        int option_char(0);
        while ((option_char = getopt (argc, argv, "i:o:s:n:e:h"))  != EOF){
            switch (option_char)
            {
                case 'i':
                    infile = optarg;
                    break;
                case 'o':
                    outfile = atoi(optarg);
                    break;
                case 's':
                    save_interval=atoi(optarg);
                    break;
                case 'n':
                    max_MCs = atoi(optarg);
                    break;
                case 'e':
                    t_equil = stod(optarg);
                    break;
                case 'h':
                default:
                    std::cerr << "\n\nUsage: ./simulatePolymerBFM << OPTIONS >> \n[-i load file] \n[-o output file] \n[-n number of total Monte Carlo steps] \n[-s save frequency (in Monte Carlo steps)]\n[-e equilibriation time]\n\n";
                    return 0;
            }
        }

        // generate ingredients
        typedef LOKI_TYPELIST_5(FeatureMoleculesIO,
                                FeatureAttributes< >,
                                FeatureExcludedVolumeSc<>,
                                FeaturePotentialBending,
                                FeatureLinearForce) Features;
        typedef ConfigureSystem<VectorInt3,Features,max_bonds> Config;
        typedef Ingredients<Config> IngredientsType;
        IngredientsType ingredients;

        // seed randomness
        RandomNumberGenerators rng;
        rng.seedAll();

        // add updaters and analyzers to task manager
        TaskManager taskmanager;
        taskmanager.addUpdater(new UpdaterReadBfmFile<IngredientsType>(infile,ingredients,UpdaterReadBfmFile<IngredientsType>::READ_LAST_CONFIG_SAVE),0);
        taskmanager.addUpdater(new UpdaterSimpleSimulator<IngredientsType,MoveLocalSc>(ingredients,save_interval));
        taskmanager.addAnalyzer(new AnalyzerWriteBfmFile<IngredientsType>(outfile,ingredients,AnalyzerWriteBfmFile<IngredientsType>::APPEND));
        taskmanager.addAnalyzer(new AnalyzerEndToEndDistance<IngredientsType>(ingredients, "RE2E.dat", t_equil));
        taskManager.addAnalyzer(new AnalyzerRadiusOfGyration<IngredientsType>(ingredients, "ROG.dat"));
        // taskmanager.addAnalyzer(new AnalyzerBondBondDistribution<IngredientsType>(ingredients, "BBD.dat", t_equil));
        taskmanager.addAnalyzer(new AnalyzerBondBondCorrelation<IngredientsType>(ingredients, "BBC.dat", t_equil));

        // if the outfile exists, delete it
        char* outfile_char_array = new char[outfile.length() + 1];
        strcpy(outfile_char_array, outfile.c_str());
        remove(outfile_char_array);

        // run
        taskmanager.initialize();
        taskmanager.run(ceil(max_MCs/save_interval));
        taskmanager.cleanup();
    }
    catch(std::exception& err){std::cerr<<err.what();}


    return 0;
}
