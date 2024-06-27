/*--------------------------------------------------------------------------------
    ooo      L   attice-based  | LeMonADE: An Open Source Implementation of the
  o\.|./o    e   xtensible     |           Bond-Fluctuation-Model for Polymers
 o\.\|/./o   Mon te-Carlo      |           
oo---0---oo  A   lgorithm and  | StiffCyclic-simulations
 o/./|\.\o   D   evelopment    | Copyright (C) 2024 by
  o/.|.\o    E   nvironment    | LeMonADE Principal Developers
    ooo      StiffCyclic       | Ron Dockhorn
----------------------------------------------------------------------------------

This file is part of LeMonADE and LeMonADE_StiffCyclic extension.

LeMonADE and LeMonADE_StiffCyclic extension is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

LeMonADE  and LeMonADE_StiffCyclic extension is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with LeMonADE and LeMonADE_StiffCyclic extension. If not, see <http://www.gnu.org/licenses/>.

--------------------------------------------------------------------------------*/

#include <cstring>
#include <stdlib.h> //for atoi
#include <unistd.h> //for getopt

#include <LeMonADE/utility/RandomNumberGenerators.h>
#include <LeMonADE/core/ConfigureSystem.h>
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/feature/FeatureAttributes.h>
#include <LeMonADE/feature/FeatureExcludedVolumeSc.h>
//#include <LeMonADE/feature/FeatureFixedMonomers.h>
#include <LeMonADE/utility/TaskManager.h>
#include <LeMonADE/updater/UpdaterReadBfmFile.h>
#include <LeMonADE/updater/UpdaterSimpleSimulator.h>

#include "UpdaterCreateTwoConcatenatedRings.h"


int main(int argc, char* argv[])
{
  try{
	//std::string infile;
	std::string outfile;
	uint32_t max_mcs=0;
	uint32_t save_interval=0;
	
	outfile="outfile.bfm";
	
	int option_char(0);

	int numMonomerInRingOne = 8;
	int numMonomerInRingTwo = 8;

	int box_size = 128;

	int crosslinker = 0;


		//read in options by getopt
		while ((option_char = getopt (argc, argv, "o:n:m:b:h"))  != EOF){
			switch (option_char)
			{
			//case 'f':
			//	infile=optarg;
			//	break;
			case 'o':
							outfile=optarg;
							break;
            case 'm':
                numMonomerInRingOne = atoi(optarg);
                break;
			case 'n':
				numMonomerInRingTwo = atoi(optarg);
				break;
			case 'b':
				box_size = atoi(optarg);
				break;
			case 'h':
			default:
				std::cerr << "Usage: " << argv[0] << "  [-o filenameOutput] [-m monomer in first ring] [-n monomer in second ring] [-b box size]\n";
				return 0;
			}
		}
		std::cout << "write out options:  filenameOutput=" << outfile <<  " Monos in first ring:" << numMonomerInRingOne<< " Monos in second ring:" << numMonomerInRingTwo << " box size" << box_size << std::endl;


	

	//seed the globally available random number generators
	RandomNumberGenerators rng;
	rng.seedAll();
	
	// FeatureExcludedVolumeSc<> is equivalent to FeatureExcludedVolumeSc<FeatureLattice<bool> >
	typedef LOKI_TYPELIST_3(FeatureMoleculesIO, FeatureAttributes<>,FeatureExcludedVolumeSc<>) Features;
	//typedef LOKI_TYPELIST_2(FeatureMoleculesIO, FeatureAttributes) Features;
	
	typedef ConfigureSystem<VectorInt3,Features, 4> Config;
	typedef Ingredients<Config> Ing;
	Ing myIngredients;
	//myIngredients.setName(outfile);

	UpdaterCreateTwoConcatenatedRings<Ing> UCRM(myIngredients, numMonomerInRingOne, numMonomerInRingTwo, box_size, box_size, box_size);
	UCRM.initialize();
	UCRM.execute();
	UCRM.cleanup();

	//taskmanager.addAnalyzer(new AnalyzerWriteBfmFile<Ing>( outfile,myIngredients));
	AnalyzerWriteBfmFile<Ing> AWBFM( outfile,myIngredients);
	AWBFM.initialize();
	AWBFM.execute();
	AWBFM.cleanup();

	}
	catch(std::exception& err){std::cerr<<err.what();}
	return 0;
  
}

