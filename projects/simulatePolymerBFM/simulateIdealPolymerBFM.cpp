#include <iostream>
#include <string>
#include <unistd.h> //for getopt
#include <cstring>
using namespace std;

#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureAttributes.h>
// #include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
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
        int type1(1); // type assigned to all "monomers"
        // establish default parameters
        double save_interval = 1000000; // frequency of property calculations, position save
        long int max_MCs = 100000000; // total number of Monte Carlo steps
        std::string infile = "config_init.bfm"; // file that contains initial configuraiton for bfm simulation
        std::string outfile = "config.bfm"; // file that contains system configuratun at each save interval
        double t_equil = 0; // simulation time before which properties are not calculated
        bool add_end2end_analyzer = false;
        bool add_hysteresis_analyzer = false;
        int num_hysanal_points = 5000;
        bool add_bondbondcorr_analyzer = false;
        bool add_bondvecdist_analyzer = false;
        bool add_radiusgyr_analyzer = false;
        bool add_scatter_analyzer = false;
        bool add_movie = false;

        // determine if any options were passed to the executable
        // read in options by getopt
        int option_char(0);
        while ((option_char = getopt (argc, argv, "abcdgmqf:o:s:n:e:h"))  != EOF){
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
                case 'm':
                    add_movie = true;
                    break;
                case 'q':
                    add_scatter_analyzer = true;
                    break;
                case 'f':
                    infile = optarg;
                    break;
                case 'o':
                    outfile = optarg;
                    break;
                case 's':
                    save_interval=atoi(optarg);
                    break;
                case 'n':
                    max_MCs = atol(optarg);
                    break;
                case 'e':
                    t_equil = stod(optarg);
                    break;
                case 'h':
                default:
                    std::cerr << "\n\nUsage: ./simulatePolymerBFM << OPTIONS >> \n[-f load file containing initial coordinates] \n[-o output file] \n[-n number of total Monte Carlo steps] \n[-s save frequency (in Monte Carlo steps)]\n[-e equilibriation time]\n[-a add end-to-end analyzer]\n[-b add hysteresis analyzer]\n[-c add bond-bond correlation analyzer]\n[-d add bond vector distribution analyzer]\n[-g add radius of gyration analyzer]\n[-m create movie]\n[-q add static scattering factor analyzer]\n\n";
                    return 0;
            }
        }

        // generate ingredients
        typedef LOKI_TYPELIST_4(FeatureMoleculesIO,
                                FeatureAttributes< >,
                                FeaturePotentialBending,
                                FeatureOscillatoryForce) Features;
        typedef ConfigureSystem<VectorInt3,Features,max_bonds> Config;
        typedef Ingredients<Config> IngredientsType;
        IngredientsType ingredients;

        // seed randomness
        RandomNumberGenerators rng;
        rng.seedAll();

        // initialize task manager, load initial config file
        TaskManager taskmanager;
        taskmanager.addUpdater(new UpdaterReadBfmFile<IngredientsType>(infile,ingredients,UpdaterReadBfmFile<IngredientsType>::READ_LAST_CONFIG_SAVE),0);
        taskmanager.initialize();

        // if hysteresis analyzer is being called, check for force oscillation and tinker with the save interval
        if (add_hysteresis_analyzer) {
            // first, check that the force is on, and that it is oscillating
            if (ingredients.isForceOscillationOn()) {
                // if force is on, adjust save interval to a value relative to oscillation period
                int p = ingredients.getForceOscillationPeriod(); // period assigned to oscillatory force
                // check if the period and the number of points are divisable by one another
                if (p % num_hysanal_points == 0) {
                    // if they are evenly divisable by one another, assign the save interval as the period over the number of save points
                    save_interval = p / num_hysanal_points;
                }
                // TODO :: Add routine for if the two numbers are not evenly divisable by one another
                // (if the two numbers are not evenly divisable, the hysteresis analyzer will throw an error)
            } else {
                // force oscillation is not on, but hysteresis analyzer is being called
                // report error and quit
                std::cout << "Hysteresis Analyzer cannot be used when the constant force is not on, or is not oscillating." << std::endl;
                std::cout << "Reinitialize config file with oscillating force or turn off Hysteresis Analyzer." << std::endl;
                exit(0);
            }
        }

        // if the outfile exists, delete it
        char* outfile_char_array = new char[outfile.length() + 1];
        strcpy(outfile_char_array, outfile.c_str());
        remove(outfile_char_array);

        // add updaters and analyzers to task manager
        taskmanager.addUpdater(new UpdaterSimpleSimulator<IngredientsType,MoveLocalSc>(ingredients,save_interval));
        if (add_movie) {
            taskmanager.addAnalyzer(new AnalyzerWriteBfmFile<IngredientsType>(outfile,ingredients,AnalyzerWriteBfmFile<IngredientsType>::APPEND));
        }
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
            // the save_interval is not added to the analyzer to control the frequency that the hysteresis analyzer samples
            // but rather just used when calculating hysteresis (as A)
            taskmanager.addAnalyzer(new AnalyzerHysteresis<IngredientsType>(ingredients, "HYS.dat", t_equil, save_interval));
        }
        if (add_scatter_analyzer) {
            taskmanager.addAnalyzer(new AnalyzerScatteringSingleObject<IngredientsType>(ingredients, "SKQ.dat", t_equil), (1));
        }

        // run
        taskmanager.initialize();
        taskmanager.run(ceil(max_MCs/save_interval));
        taskmanager.cleanup();
    }
    catch(std::exception& err){std::cerr<<err.what();}


    return 0;
}
