/******************************************************************************
 *   Copyright (C) 2014 - 2023 Jan Fostier (jan.fostier@ugent.be)             *
 *   This file is part of HaploDetox                                          *
 *                                                                            *
 *   This program is free software; you can redistribute it and/or modify     *
 *   it under the terms of the GNU Affero General Public License as published *
 *   by the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                      *
 *                                                                            *
 *   This program is distributed in the hope that it will be useful,          *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *   GNU General Public License for more details.                             *
 *                                                                            *
 *   You should have received a copy of the GNU Affero General Public License *
 *   along with this program; if not, see <https://www.gnu.org/licenses/>.    *
 ******************************************************************************/

#include "settings.h"
#include "util.h"

#include <thread>
#include <iostream>
#include <fstream>
#include <string.h>
#include <numeric>

using namespace std;

void Settings::printProgramVersion() const
{
        cout << "HaploDetox, determines strain-aware de Bruijn graphs via node/edge multiplicity estimation \n\t-- version "
             << HAPLODETOX_MAJOR_VERSION << "." << HAPLODETOX_MINOR_VERSION
             << "." << HAPLODETOX_PATCH_LEVEL << " (commit " << GIT_SHA1 << ")" << "\n";

        cout << "Copyright (C) 2019-2023 Jan Fostier, Aranka Steyaert and Marie Van Hecke\n"
                "Report bugs to Jan Fostier <jan.fostier@ugent.be>\n"
                "This is free software; see the source for copying conditions. "
                "There is NO\nwarranty; not even for MERCHANTABILITY or "
                "FITNESS FOR A PARTICULAR PURPOSE." << endl;
}

void Settings::printUsage() const
{
        cout << "HaploDetox, determines strain-aware de Bruijn graphs via node/edge multiplicity estimation \n\t-- version "
        << HAPLODETOX_MAJOR_VERSION << "." << HAPLODETOX_MINOR_VERSION
        << "." << HAPLODETOX_PATCH_LEVEL << " (commit " << GIT_SHA1 << ")" << "\n";
        
        cout <<

        "Usage: haplodetox [options] unitigs.fa readfile\n\n  File "
        "\"unitigs.fa\" is the input de Bruijn graph produced by BCALM 2\n"
        "  File \"readfile\" contains a list of input FASTQ files\n\n"

        " [options without arg]\n"
        "  -no-correct\t\tSkip the removal of sequencing errors based on initial model in stage 2\n"
        "  -use-qual\t\tUse Phred quality scores to weigh k-mer coverage. When using this\n"
        "            \t\toption, it is advised to also pass the -abundance-min parameter\n"
        "  -no-approx-inf\tUse exact inference on subgraph-based CRFs for all multiplicity computations\n"
        "  -map-assignment\tUse in combination with approx-inf; compute a MAP assignment\n"
        "  -em-no-eqweight\tDo not equalise certain (strain-specific) weights in model training\n"
        "                \tEqualising weights usually improves model fit,"
        "                \tbut can be switched of with this option\n"
        "  -em-fix-zero\t\tFix the error distribution as fitted during initialisation\n"
        "               \t\tinstead of refitting after graph cleaning\n\n"
//        "  -mm-subgenus\t\tUse an underlying histogram model that assumes two strains differ above species level\n"
//        "               \t\t(not advised to use haplodetox for differences above genus level)[default=false]\n"
//        "              \t\tWARNING: option under development\n" 
        "  -help\t\t\tdisplay help page\n" 
        "  -version\t\tdisplay version\n\n"

        " [options with 1 arg]\n"
        "  -num-threads\t\tNumber of threads [default = #cores]\n\n"

        "  -abundance-min\tMin abundance threshold as passed to BCALM 2 [default = auto]\n\n"
        
        "  -max-cut-frac\t\tMaximum fraction of estimated avgerage coverage below\n"
        "               \t\twich we check for erroneous nodes/arcs during graph cleaning [default=0.05]\n\n"

        "  -crf-nb-size\t\tCRF neighborhood size [default = 3]\n"
        "  -crf-max-mult\t\tCRF maximum multiplicity [default = 2]\n"
        "  -crf-flow\t\tCRF flow conservation strength [default = 1.0e7]\n"
        "  -crf-max-fact\t\tCRF maximum factor size [default = 1.0e6]\n\n"

        "  -mm-coverage\t\tInitial coverage est. in mixture model [default = auto]\n"
        "  -mm-err-cov\t\tInitial error coverage est. in mixture model [default = 1.0]\n"
        "  -mm-odf\t\tInitial overdispersion factor est. in mixture model [default = 1.5]\n"
        "  -mm-frac\t\tInitial est. fractions at which multiple strains present might occur\n"
        "          \t\tformat: f1;f2;...;fn [default = equal distribution of fractions]\n"
        "  -mm-initstrains\tInitial number of strains estimate\n"
        "               \t\t(nstrains in [initstrains, initstrains+2] will be tested). [default = 2]\n"
        "  -mm-nstrains\t\tKnown number of strains, if this parameter is set, no test of\n"
        "              \t\tnumber of strains is performed, but a model with this number will be trained\n"
        "              \t\tIf mm-frac is set, this parameter will be overwritten with\n"
        "              \t\tthe length of the fractions vector [default = -1]\n\n"

        "  -em-max-iter\t\tMaximum number of EM iterations [default = 25]\n"
        "  -em-conv-eps\t\tRelative EM convergence epsilon [default = 1.0e-6]\n"
        "  -em-max-change\tEM convergence number of assignment changes epsilon [default  =1.0e-3]\n"
        "  -em-train-size\tNumber of nodes/arcs to use for EM training\n"
        "                \t(-1 for training on all nodes/arcs) [default = 1e4]\n\n"

        "  -phred-base\t\tASCII value corresponding to Phred score Q = 0 [default = 33]\n\n"
        "  -vis-subgraph\t\tCentral nodeID around which a subgraph will be visualised [default = 0]\n"
        "               \t\t if value = 0, no visualisation and normal pipeline will be run\n\n"
        

        "Report bugs to Jan Fostier <jan.fostier@ugent.be>" << endl;
}


Settings::Settings(int argc, char ** argv) : correctFlag(true),
        useQualFlag(false), approxInfFlag(true), mapAssignment(false), singleCRF(false), numThreads(thread::hardware_concurrency()),
        abundanceMin(-1), maxCutoffFrac(0.05), crfDepth(3), crfMaxMult(2), crfFlow(1e7), crfMaxFact(1e6), mmCov(-1.0),
        mmErrCov(1.0), mmODF(1.5), mmErrODF(1.0), mmStrains(-1), initStrains(2), mmGenus(false), emMaxIter(25),
        emConvEps(1e-6), emMaxChange(1e-3), emTrainSize(1e4), emEqWeights(true), emFixZero(false), phredBase(33), visGraphNode(0)
{
        const int reqArguments = 2;     // not counting argument 0

        //mmFrac = {0.3, 0.7};

        // no arguments provided
        if (argc <= 1) {
                printUsage();
                exit(EXIT_FAILURE);
        }

        // only one argument provided
        if (argc == 2) {
                string arg(argv[1]);

                if (arg == "-help") {
                        printUsage();
                        exit(EXIT_SUCCESS);
                } else if (arg == "-version") {
                        printProgramVersion();
                        exit(EXIT_SUCCESS);
                } else {
                        printUsage();
                        exit(EXIT_FAILURE);
                }
        }

        // process optional arguments (we know there are at least two arguments)
        for (int i = 1; i < argc-reqArguments; i++) {
                string arg(argv[i]);

                // process options without arguments
                if (arg == "-no-correct") {
                        correctFlag = false;
                        continue;
                } else if (arg == "-use-qual") {
                        useQualFlag = true;
                        continue;
                } else if (arg == "-em-no-eqweight") {
                        emEqWeights = false;
                        continue;
                } else if (arg == "-em-fix-zero") {
                        emFixZero = true;
                        continue;
                } else if (arg == "-no-approx-inf") {
                        approxInfFlag = false;
                        continue;
                } else if (arg == "-map-assignment") {
                        mapAssignment = true;
                        continue;
                } else if (arg == "-single-crf"){
                        singleCRF = true;
                } else if (arg == "-mm-subgenus") {
                        mmGenus = true;
                        continue;
                } else if (arg == "-help") {
                        printUsage();
                        exit(EXIT_SUCCESS);
                } else if (arg == "-version") {
                        printProgramVersion();
                        exit(EXIT_SUCCESS);
                } else if (i == argc-reqArguments-1) {
                        cerr << "Unknown option or missing argument: "
                             << argv[i] << endl;
                        printUsage();
                        exit(EXIT_FAILURE);
                }

                // process options with arguments
                if (arg == "-num-threads") {
                        numThreads = atoi(argv[i+1]);
                } else if (arg == "-abundance-min") {
                        abundanceMin = atoi(argv[i+1]);
                } else if (arg == "-max-cut-frac") {
                        maxCutoffFrac = atof(argv[i+1]);
                } else if (arg == "-crf-nb-size") {
                        crfDepth = atoi(argv[i+1]);
		} else if (arg == "-crf-max-mult") {
			crfMaxMult = atoi(argv[i+1]);
                } else if (arg == "-crf-flow") {
                        crfFlow = atof(argv[i+1]);
                } else if (arg == "-crf-max-fact") {
                        crfMaxFact = atof(argv[i+1]);
                } else if (arg == "-mm-coverage") {
                        mmCov = atof(argv[i+1]);
                } else if (arg == "-mm-err-cov") {
                        mmErrCov = atof(argv[i+1]);
                } else if (arg == "-mm-odf") {
                        mmODF = atof(argv[i+1]);
                } else if (arg == "-mm-frac"){
                        mmFrac.clear();
                        char* token;
                        token = strtok(argv[i+1], ";");
                        while(token != NULL){
                                mmFrac.push_back(atof(token));
                                token = strtok(NULL, ";");
                        }
                } else if (arg == "-mm-nstrains") {
                        mmStrains = atoi(argv[i+1]);
                } else if (arg == "mm-initstrains") {
                        initStrains = atoi(argv[i+1]);
                } else if (arg == "-em-max-iter") {
                        emMaxIter = atoi(argv[i+1]);
                } else if (arg == "-em-conv-eps") {
                        emConvEps = atof(argv[i+1]);
                } else if (arg == "-em-max-change") {
                        emMaxChange = atof(argv[i+1]);
                } else if (arg == "-em-train-size") {
                        emTrainSize = atof(argv[i+1]);
                } else if (arg == "-phred-base") {
                        phredBase = atoi(argv[i+1]);
                } else if (arg == "-vis-subgraph") {
                        visGraphNode = atoi(argv[i+1]);
                } else {
                        cerr << "Unknown argument: " << argv[i] << endl;
                        printUsage();
                        exit(EXIT_FAILURE);
                }

                i++;    // if we reach this point, an argument was processed
        }

        if ( ! mmFrac.empty() ) { // If fractions were passed set nstrains to length of fractions vector
                mmStrains = mmFrac.size();
        } else if (mmStrains > 0) { // If only number of strains was passed use evenly spaced fraction estimates
                mmFrac.resize(mmStrains); // TODO: check what is really the best way to initialise fractions
                int sum = 0;
                for (int i = 0; i < mmStrains; i++) {
                        mmFrac[i] = i+1;
                        sum += (i+1);
                }
                for (int i = 0; i < mmStrains; i++)
                        mmFrac[i] /= sum;
        }

        graphFilename = argv[argc-2];
        readFilename = argv[argc-1];

        // check graph file
        if (graphFilename.empty()) {
                cerr << "Specify input de Bruijn graph file" << endl;
                printUsage();
                exit(EXIT_FAILURE);
        }

        if (!Util::fileExists(graphFilename))
                throw runtime_error("cannot open file " + graphFilename);

        // check read file
        if (readFilename.empty()) {
                cerr << "Specify input read file" << endl;
                printUsage();
                exit(EXIT_FAILURE);
        }

        if (!Util::fileExists(readFilename))
                throw runtime_error("cannot open file " + readFilename);
        
        // perform sanity checks on input paramters
        if (numThreads < 1)
                throw runtime_error("Number of threads must be >= 1");
        if (crfDepth < 0 || crfDepth > 10)
                if (! approxInfFlag)
                        throw runtime_error("CRF subgraph depth must be "
                                            "in range [0..10]");
        if (mapAssignment && ! approxInfFlag)
                throw runtime_error("MAP assignments not supported with exact inference on neightbourhoods");
	if (crfMaxMult < 2 || crfMaxMult > 10)
		throw runtime_error("CRF maximum multiplicity must be "
				    "in range [2..10]");
        if (crfFlow < 1.0)
                throw runtime_error("CRF flow conservation strength "
                                    "must be >= 1.0");
        if (crfMaxFact < 1e3)
                throw runtime_error("CRF maximum factor size must "
                                    "be >= 1e3");
        // Note: negative avgCov means "auto detect" and is hence allowed
        if (mmCov >= 0.0 && mmCov < 1.0)
                throw runtime_error("Mixture model initial coverage "
                                    "must be >= 1.0");
        if (mmCov >= 0.0 && mmCov <= mmErrCov)
                throw runtime_error("Mixture model initial coverage must be "
                                    "higher than the initial error coverage");
        if (mmErrCov <= 0.0)
                throw runtime_error("Mean sequencing error coverage "
                                    "must be > 0.0");
        if (mmODF < 1.0)
                throw runtime_error("Mixture model initial overdispersion "
                                    "factor must be >= 1.0");
        if (emMaxIter < 1)
                throw runtime_error("Maximum number of EM iterations "
                                    "must be >= 1");
        if (maxCutoffFrac > 0.5 || maxCutoffFrac < 1e-8)
                throw runtime_error("Maximum fraction of estimated average coverage to determine error coverage upper bound"
                                        "should be in range [1e-8..0.5]");
        if (emConvEps >= 1.0 || emConvEps < 1e-8)
                throw runtime_error("EM relative converage tolerance "
                                    "should be in range [1e-8..1[");
        
        if (emMaxChange >= 1.0 || emMaxChange < 1e-8)
                throw runtime_error("EM relative converage tolerance "
                                    "should be in range [1e-8..1[");
        if (emTrainSize < 1e2 && emTrainSize!=-1)
                throw runtime_error("EM traing size should be >= 1e2");
        if (phredBase < 0 || phredBase > 255)
                throw runtime_error("ASCII value for Phred score Q = 0 must "
                                    "be in range [0..255]");
        if (useQualFlag && abundanceMin == -1) {
                cerr << "WARNING: no -abundance-min argument passed while using q-mers.\n"
                     << "Detox cannot automatically infer the correct value.\n"
                     << "Please consider providing the appropriate -abudance-min argument to Detox.\n";
        }

}
