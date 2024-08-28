#include <iostream>
#include <string>
#include <unistd.h> //for getopt
#include <cstring>
using namespace std;

#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
// #include <LeMonADE/feature/FeatureLinearForce.h>

#include <LeMonADE/updater/UpdaterAddLinearChains.h>
#include <LeMonADE/updater/UpdaterSimpleSimulator.h>
#include <LeMonADE/updater/UpdaterReadBfmFile.h>
#include <LeMonADE/analyzer/AnalyzerWriteBfmFile.h>
#include <LeMonADE/analyzer/AnalyzerRadiusOfGyration.h>

#include <LeMonADE/utility/RandomNumberGenerators.h>
#include <LeMonADE/utility/TaskManager.h>
#include <LeMonADE/utility/Vector3D.h>

// modules not found in LeMonADE Library
// #include <clara>
#include <FeaturePotentialBending.h>
#include <FeatureOscillatoryForce.h>
#include <AnalyzerEndToEndDistance.h>
#include <AnalyzerBondBondCorrelation.h>
#include <AnalyzerBondVectorDistribution.h>
#include <AnalyzerHysteresis.h>
#include <AnalyzerScatteringSingleObject.h>

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
        double t_equil = 0; // simulation time before which properties are not calculated
        bool add_end2end_analyzer = false;
        bool add_hysteresis_analyzer = false;
        bool add_bondbondcorr_analyzer = false;
        bool add_bondvecdist_analyzer = false;
        bool add_radiusgyr_analyzer = false;
        bool add_scatter_analyzer = false;

        // determine if any options were passed to the executable
        // read in options by getopt
        int option_char(0);
        while ((option_char = getopt (argc, argv, "abcdgqi:o:s:n:e:h"))  != EOF){
            switch (option_char)
            {
                case 'a':
                    add_end2end_analyzer = true;
                    break;
                case 'b':
                    add_hysteresis_analyzer = true;
                    break;
                case 'c':
                    add_bondbondcorr_analyzer = true;
                    break;
                case 'd':
                    add_bondvecdist_analyzer = true;
                    break;
                case 'g':
                    add_radiusgyr_analyzer = true;
                    break;
                case 'q':
                    add_scatter_analyzer = true;
                    break;
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
                    std::cerr << "\n\nUsage: ./simulatePolymerBFM << OPTIONS >> \n[-i load file] \n[-o output file] \n[-n number of total Monte Carlo steps] \n[-s save frequency (in Monte Carlo steps)]\n[-e equilibriation time]\n[-a add end-to-end analyzer]\n[-b add hysteresis analyzer]\n[-c add bond-bond correlation analyzer]\n[-d add bond vector distribution analyzer]\n[-g add radius of gyration analyzer]\n[-q add static scattering factor analyzer]\n\n";
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

        // add updaters and analyzers to task manager
        TaskManager taskmanager;
        taskmanager.addUpdater(new UpdaterReadBfmFile<IngredientsType>(infile,ingredients,UpdaterReadBfmFile<IngredientsType>::READ_LAST_CONFIG_SAVE),0);
        taskmanager.addUpdater(new UpdaterSimpleSimulator<IngredientsType,MoveLocalSc>(ingredients,save_interval));
        taskmanager.addAnalyzer(new AnalyzerWriteBfmFile<IngredientsType>(outfile,ingredients,AnalyzerWriteBfmFile<IngredientsType>::APPEND));
        if (add_end2end_analyzer) {
            taskmanager.addAnalyzer(new AnalyzerEndToEndDistance<IngredientsType>(ingredients, "RE2E.dat", t_equil));
        }
        if (add_radiusgyr_analyzer) {
            taskmanager.addAnalyzer(new AnalyzerRadiusOfGyration<IngredientsType>(ingredients, "ROG.dat"));
        }
        if (add_bondbondcorr_analyzer) {
            taskmanager.addAnalyzer(new AnalyzerBondVectorDistribution<IngredientsType>(ingredients, "BVD.dat", t_equil));
        }
        if (add_bondvecdist_analyzer) {
            taskmanager.addAnalyzer(new AnalyzerBondBondCorrelation<IngredientsType>(ingredients, "BBC.dat", t_equil));
        }
        if (add_hysteresis_analyzer) {
            // TODO :: here is the problem. the save interval is fixed with the UpdaterSimpleSimulator. The save_interval is only passed here
            // to the AnalyzerHysteresis in order to calculate hysteresis (as A). Futher, the period is unspecified in ingredients until the
            // initialize function is called. This means that right now, in order for the save_interval and the period to synchronize to a divisable value,
            // the program calling generate and simulate PolymerBFM need to be speaking to each other.
            // At the very least, it would be nice if the simulatePolymerBFM could add a hysteresis analyzer and then work without collapsing, even if it
            // isnt optimal. Can I initialize the taskmanager / ingredients, and then re-analyze? Can I provide suggestions for save interval
            // if it does not work?
            taskmanager.addAnalyzer(new AnalyzerHysteresis<IngredientsType>(ingredients, "HYS.dat", t_equil, save_interval));
        }
        if (add_scatter_analyzer) {
            taskmanager.addAnalyzer(new AnalyzerScatteringSingleObject<IngredientsType>(ingredients, "SKQ.dat", t_equil), (1));
        }

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
