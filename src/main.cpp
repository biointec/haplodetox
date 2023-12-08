/******************************************************************************
 *   Copyright (C) 2014 - 2023 Jan Fostier (jan.fostier@ugent.be)             *
 *   This file is part of HaploHaploDetox                                          *
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

#include <cstdlib>
#include <numeric>
#include <iomanip>

#include "dbgraph.h"
#include "settings.h"
#include "refcomp.h"
#include "crfmult.h"
#include "readaln.h"
#include "util.h"

using namespace std;

void populateNodeMult(NodeMap<Multiplicity>& nodeMult,
                      const vector<NodeRep>& nodes, size_t numStrains)
{
        nodeMult = NodeMap<Multiplicity>(nodes.size());
        for (size_t i = 0; i < nodes.size(); i++)
                nodeMult[nodes[i]] = Multiplicity(vector<int>(numStrains, 0));
}

void populateEdgeMult(EdgeMap<Multiplicity>& edgeMult,
                      const vector<EdgeRep>& edges, size_t numStrains)
{
        edgeMult = EdgeMap<Multiplicity>(edges.size());
        for (size_t i = 0; i < edges.size(); i++)
                edgeMult[edges[i]] = Multiplicity(vector<int>(numStrains, 0));
}

void writeTrueStrainMultiplicities(const DBGraph& dBG, 
                                  const Settings& settings, 
                                  const KmerNPPTable& table,
                                  NodeMap<vector<int>>& trueNodeMult,
                                  EdgeMap<vector<int>>& trueEdgeMult,
                                  bool correct)
{
        // Get all strain reference filenames
        int numstrains = 0;
        vector<string> strainRefs;
        
        ifstream ifs("strains.mf");
        if (!ifs){
                cout << "cannot open file strains.mf, no true multiplicity files will be created." << endl ;
                return;
        }
        
        string line;
        while (getline(ifs, line)) {
                if (line.empty())
                        continue;
                strainRefs.push_back(line);
                numstrains++;
        }
        
        // Get all reference comparisons for all strains refs
        trueNodeMult.clear();
        for (NodeRep n: dBG.getNodeReps(dBG.getNumValidNodes())){
                trueNodeMult[n] = vector<int>(numstrains, 0);
        }
        trueEdgeMult.clear();
        for (EdgeRep e: dBG.getEdgeReps(dBG.getNumValidArcs())){
                trueEdgeMult[e] = vector<int>(numstrains, 0);
        }
        
        for (int i=0; i<strainRefs.size(); i++){
                RefComp refComp(dBG, settings, table, strainRefs[i]);
                refComp.ignoreLowerCase();
                refComp.getTrueStrainMultiplicity(trueNodeMult, trueEdgeMult, i);
        }
        
        ofstream nodeOFS(settings.getTrueNodeStrainFilename(correct));
        ofstream edgeOFS(settings.getTrueEdgeStrainFilename(correct));
        
        
        // write the true node multiplicities
        for (const auto& it : trueNodeMult){
                nodeOFS << it.first; 
                for (int mult: it.second)
                        nodeOFS << "\t" << mult;
                nodeOFS << "\n";
        }
        nodeOFS.close();
        
        // write the true edge multiplicities
        for (const auto& it : trueEdgeMult){
                edgeOFS << it.first.getSrcID() << "\t"
                << it.first.getDstID();
                for (int mult: it.second)
                        edgeOFS << "\t" << mult;
                edgeOFS << "\n";
        }
        edgeOFS.close();
        
}

void readTrueStrainMultiplicities(const Settings& settings, 
                                  const vector<NodeRep>& nodeReps,
                                  const vector<EdgeRep>& edgeReps, 
                                  NodeMap<vector<int>>& trueNodeMult,
                                  EdgeMap<vector<int>>& trueEdgeMult,
                                  bool correct)
{
        ifstream ifs_("strains.mf");
        int trueNstrains = 0;
        if (!ifs_){
                cerr << "cannot open file strains.mf, assuming you passed the true number of strains via settings\n";
                trueNstrains = settings.getInitNStrains();
        } else {
                string line;
                while (getline(ifs_, line)) {
                        if (line.empty())
                                continue;
                        trueNstrains++;
                }
        }
        
        for (const auto& nr : nodeReps)
                trueNodeMult[nr] = vector<int>(trueNstrains, -1);
        
        for (const auto& er : edgeReps)
                trueEdgeMult[er] = vector<int>(trueNstrains, -1);
        
        ifstream ifs(settings.getTrueNodeStrainFilename(correct));
        if (!ifs)
                cerr << "Could not find true nodes multiplicities file...\n"
                "True node strain multiplicities will be set to -1 in "
                "resulting Cytoscape graph\n";
        
        while (ifs) {
                NodeID nodeID; int m; vector<int> strainM;
                ifs >> nodeID;
                for (int i=0; i < trueNstrains; i++){
                        ifs >> m;
                        strainM.push_back(m);
                }
                if (!ifs)
                        break;
                auto it = trueNodeMult.find(nodeID);
                if (it != trueNodeMult.end())
                        it->second = strainM;
        }
        ifs.close();
        
        ifs.open(settings.getTrueEdgeStrainFilename(correct));
        if (!ifs)
                cerr << "Could not find true edges multiplicities file...\n"
                "True edge multiplicities will be set to -1 in "
                "resulting Cytoscape graph\n";
        
        while (ifs) {
                NodeID srcID, dstID; int m; vector<int> strainM;
                ifs >> srcID >> dstID;
                for (int i=0; i < trueNstrains; i++){
                        ifs >> m;
                        strainM.push_back(m);
                }
                if (!ifs)
                        break;
                auto it = trueEdgeMult.find(EdgeRep(srcID, dstID));
                if (it != trueEdgeMult.end())
                        it->second = strainM;
        }
        ifs.close();
}

void stageOne(Settings& settings, LibraryContainer& libraries)
{
        cout << "\nEntering stage 1\n";
        cout << "================\n" << endl;

        if (!settings.stageOneNecessary()) {
                cout << "File " << settings.getStage1GraphFilename()
                     << " exists. Skipping stage 1...\n";
                return;
        }

        Util::startChrono();

        DBGraph dBG(settings);

        // load the graph from BCALM2 file
        dBG.loadBCalm(settings.getGraphFilename());

        // create a <kmer, nodePosPair> table
        KmerNPPTable table(dBG, true);
        dBG.createKmerNPPTable(table);

        // stream the reads to get the coverage
        dBG.getCovFromReads(libraries, table);
        
        // write stage 1 binary output file
        dBG.writeBinary(settings.getStage1GraphFilename());

        if (Util::fileExists("strains.mf")) {
                NodeMap<vector<int>> trueNodeMult(dBG.getNumValidNodes());
                EdgeMap<vector<int>> trueEdgeMult(dBG.getNumValidArcs()/2);
                writeTrueStrainMultiplicities(dBG, settings, table,
                                             trueNodeMult,
                                             trueEdgeMult, false);
                
                vector<NodeID> nodes;
                vector<pair<NodeID, NodeID>> edges;
                nodes.clear();
                edges.clear();
                dBG.getFullDirGraph(nodes, edges);

                dBG.writeStage1CytoscapeGraph("Cytograph_strains", nodes, edges, trueNodeMult, trueEdgeMult);

                dBG.writeStage1GFA("stage1", nodes, edges, trueNodeMult, trueEdgeMult);

        }

        cout << "Stage 1 finished in " << Util::stopChronoStr() << endl;
}

void stageTwoCorrectGraph(Settings& settings, const CovModel& nodeModel,
                          const CovModel& edgeModel, DBGraph& dBG,
                          int& abundanceMin, int& numNodesRemoved, bool write=true)
{
        cout << "\nCorrecting graph\n";
        cout << "================\n" << endl;
        cout << "\tGraph has " << dBG.getNumValidNodes() << " nodes and "
             << dBG.getNumValidArcs() << " arcs" << endl;

        Util::startChrono();

	numNodesRemoved = 0;

        CRFMult myCRFMult(dBG, settings.getCRFDepth(),
                          settings.getCRFMaxFactSize(),
                          settings.getCRFFlowStrength(),
                          settings.approxInf(),
                          settings.getNumThreads(),
                          settings.getThreadGraphWorkSize());

        // remove arcs with zero coverage
        cout << "Removing nodes/arcs with zero coverage..." << endl;
        int oldNumValidNodes = dBG.getNumValidNodes();
	dBG.removeCoverage(0.0, 0);
	numNodesRemoved += oldNumValidNodes - dBG.getNumValidNodes();
        dBG.concatenateNodes();

        cout << "\tGraph has " << dBG.getNumValidNodes() << " nodes and "
             << dBG.getNumValidArcs() << " arcs" << endl;

        // If they exist, also load the true multiplicities from disk
        // Note: for efficiency, story only true multiplicities that are needed
        
        vector<NodeRep> nodeReps = dBG.getNodeReps(dBG.getNumNodes());
        vector<EdgeRep> edgeReps = dBG.getEdgeReps(dBG.getNumArcs());

        NodeMap<vector<int>> trueNodeMult;
        EdgeMap<vector<int>> trueEdgeMult;
        readTrueStrainMultiplicities(settings, 
                                     nodeReps, edgeReps,
                                     trueNodeMult, trueEdgeMult,
                                     false);

        const int numRounds = 4;
	double nodeCO = min(nodeModel.getCovCutOff(1), nodeModel.getSosLambda()*settings.maxcutF());
        for (int i = 1; i <= numRounds; i++)
        {
                double fraction = double(i) / double(numRounds);
                double nodeCutoff = fraction * nodeCO;
                abundanceMin = nodeCutoff+1;

                // a) TIPS
                {
                        vector<NodeRep> nodeReps = dBG.getLowCovSmallTips(nodeCutoff, 2*Kmer::getK());
                        cout << "Selected " << nodeReps.size()
                        << " tips with coverage <= " << nodeCutoff << endl;

                        vector<bool> flowOK(nodeReps.size());
                        myCRFMult.checkFlow(nodeReps, flowOK, nodeModel, edgeModel);

                        size_t numRemove = count(flowOK.begin(), flowOK.end(), true);

                        vector<NodeRep> toRemove;
                        toRemove.reserve(numRemove);

                        for (size_t i = 0; i < nodeReps.size(); i++)
                                if (flowOK[i]){
                                        toRemove.push_back(nodeReps[i]);
#ifdef DEBUG
                                        if (accumulate(trueNodeMult[nodeReps[i]].begin(), trueNodeMult[nodeReps[i]].end(), 0) > 0)
                                                cout << "ATTENTION !!! REMOVING TRUE NODE "  << nodeReps[i].getNodeID() << ", (marg length, cov): " <<  dBG.getDSNode(nodeReps[i].getNodeID()).getMarginalLength()<< ", " << dBG.getDSNode(nodeReps[i].getNodeID()).getCov() << endl;
#endif
                                }

                cout << "\tRemoving " << numRemove << " nodes" << endl;
                numNodesRemoved += numRemove;
		dBG.removeNodes(toRemove);
                dBG.concatenateNodes();

                cout << "\tGraph has " << dBG.getNumValidNodes() << " nodes and "
                <<  dBG.getNumValidArcs() << " arcs" << endl;
                }

                // b) BUBBLES
                {
                        vector<NodeRep> nodeReps = dBG.getLowCovBubbles(nodeCutoff);
                        cout << "Selected " << nodeReps.size()
                             << " bubbles with coverage <= " << nodeCutoff << endl;

                        vector<bool> flowOK(nodeReps.size());
                        myCRFMult.checkFlow(nodeReps, flowOK, nodeModel, edgeModel);

                        size_t numRemove = count(flowOK.begin(), flowOK.end(), true);
                        vector<NodeRep> toRemove;
                        toRemove.clear();
                        toRemove.reserve(numRemove);

                        cout << "\tRemoving " << numRemove << " nodes" << endl;
                        for (size_t i = 0; i < nodeReps.size(); i++)
                                if (flowOK[i]){
                                        toRemove.push_back(nodeReps[i]);
#ifdef DEBUG
                                        if (accumulate(trueNodeMult[nodeReps[i]].begin(), trueNodeMult[nodeReps[i]].end(), 0) > 0)
                                                cout << "ATTENTION !!! REMOVING TRUE NODE "  << nodeReps[i].getNodeID() << ", (marg length, cov): " <<  dBG.getDSNode(nodeReps[i].getNodeID()).getMarginalLength()<< ", " << dBG.getDSNode(nodeReps[i].getNodeID()).getCov() << endl;
#endif
                                }
                        
                        dBG.removeNodes(toRemove);
                        dBG.concatenateNodes();
                        cout << "\tGraph has " << dBG.getNumValidNodes() << " nodes and "
                             <<  dBG.getNumValidArcs() << " arcs" << endl;
                }
                
                // c) ALL
                {
                        vector<NodeRep> nodeReps = dBG.getLowCovNodes(nodeCutoff);
                        cout << "Selected " << nodeReps.size()
                        << " nodes with coverage <= " << nodeCutoff << endl;
                        
                        vector<bool> flowOK(nodeReps.size());
                        myCRFMult.checkFlow(nodeReps, flowOK, nodeModel, edgeModel);
                        
                        size_t numRemove = count(flowOK.begin(), flowOK.end(), true);
                        vector<NodeRep> toRemove;
                        toRemove.clear();
                        toRemove.reserve(numRemove);
                        
                        cout << "\tRemoving " << numRemove << " nodes" << endl;
                        for (size_t i = 0; i < nodeReps.size(); i++)
                                if (flowOK[i]){
                                        toRemove.push_back(nodeReps[i]);
#ifdef DEBUG
                                        if (accumulate(trueNodeMult[nodeReps[i]].begin(), trueNodeMult[nodeReps[i]].end(), 0) > 0)
                                                cout << "ATTENTION !!! REMOVING TRUE NODE "  << nodeReps[i].getNodeID() << ", (marg length, cov): " <<  dBG.getDSNode(nodeReps[i].getNodeID()).getMarginalLength()<< ", " << dBG.getDSNode(nodeReps[i].getNodeID()).getCov() << endl;
#endif
                                }
                                
                        dBG.removeNodes(toRemove);
                        dBG.concatenateNodes();
                        
                        cout << "\tGraph has " << dBG.getNumValidNodes()
                        << " nodes and " <<  dBG.getNumValidArcs() << " arcs\n";
                }
        }

        
        for (int i = 1; i <= numRounds; i++)
        {
                double fraction = double(i) / double(numRounds);
                double covCutoff = fraction * nodeCO;
                
                // a) nodes
                vector<NodeRep> nodeReps = dBG.getLowCovNodes(covCutoff);
                cout << "Selected " << nodeReps.size()
                << " nodes with coverage <= " << covCutoff << endl;
                
                NodeMap<Multiplicity> nodeMult(nodeReps.size());
                for (size_t i = 0; i < nodeReps.size(); i++)
                        nodeMult[nodeReps[i]] = Multiplicity();
                EdgeMap<Multiplicity> edgeMult;
                myCRFMult.computeMult(nodeMult, nodeModel, edgeMult, edgeModel);
                
                vector<NodeRep> toRemove;
                for (const auto& it : nodeMult){
                        vector<int> expM = it.second.getExpMult();
                        if (accumulate(expM.begin(), expM.end(), 0) == 0)
                                toRemove.push_back(it.first);
                        if (accumulate(trueNodeMult[it.first].begin(), trueNodeMult[it.first].end(), 0) > 0){
#ifdef DEBUG
                                cout << "ATTENTION !!! REMOVING TRUE NODE "  << it.first.getNodeID() 
                                        << ", (marg length, cov): " <<  dBG.getDSNode(it.first.getNodeID()).getMarginalLength()
                                        << ", " << dBG.getDSNode(it.first.getNodeID()).getCov() << endl;
#endif
                        }
                }
                                
                        
                cout << "\tRemoving " << toRemove.size() << " nodes" << endl;
                dBG.removeNodes(toRemove);
                dBG.concatenateNodes();
                
                // b) edges
                vector<EdgeRep> edgeReps = dBG.getLowCovEdges(covCutoff);
                cout << "Selected " << edgeReps.size()
                << " edges with coverage <= " << covCutoff << endl;
                
                nodeMult.clear();
                edgeMult = EdgeMap<Multiplicity>(edgeReps.size());
                for (size_t i = 0; i < edgeReps.size(); i++)
                        edgeMult[edgeReps[i]] = Multiplicity();
                
                myCRFMult.computeMult(nodeMult, nodeModel, edgeMult, edgeModel);
                
                vector<EdgeRep> edgesToRemove;
                for (const auto& it : edgeMult){
                        vector<int> expM = it.second.getExpMult();
                        if (accumulate(expM.begin(), expM.end(), 0) == 0)
                                edgesToRemove.push_back(it.first);
                }
                        
                cout << "\tRemoving " << edgesToRemove.size() << " arcs" << endl;
                dBG.removeEdges(edgesToRemove);
                dBG.concatenateNodes();
                
                cout << "\tGraph has " << dBG.getNumValidNodes()
                << " nodes and " <<  dBG.getNumValidArcs() << " arcs\n";
        }

        
        {
                vector<NodeRep> nodeReps = dBG.getLowCovSmallTips(1e100,
                                                             2*Kmer::getK());
                cout << "Selected " << nodeReps.size() << " tips " << endl;
                
                NodeMap<Multiplicity> nodeMult(nodeReps.size());
                for (size_t i = 0; i < nodeReps.size(); i++)
                        nodeMult[nodeReps[i]] = Multiplicity();
                EdgeMap<Multiplicity> edgeMult;
                myCRFMult.computeMult(nodeMult, nodeModel,
                                      edgeMult, edgeModel);
                
                vector<NodeRep> toRemove;
                for (const auto& it : nodeMult){
                        vector<int> expM = it.second.getExpMult();
                        if (accumulate(expM.begin(), expM.end(), 0) == 0)
                                toRemove.push_back(it.first);
#ifdef DEBUG
                        if (accumulate(trueNodeMult[it.first].begin(), trueNodeMult[it.first].end(), 0) > 0){
                                cout << "ATTENTION !!! REMOVING TRUE NODE "  << it.first.getNodeID() 
                                << ", (marg length, cov): " <<  dBG.getDSNode(it.first.getNodeID()).getMarginalLength()
                                << ", " << dBG.getDSNode(it.first.getNodeID()).getCov() << endl;
                        }
#endif
                }
                
                cout << "\tRemoving " << toRemove.size() << " nodes" << endl;
                dBG.removeNodes(toRemove);
                dBG.concatenateNodes();
                
                cout << "\tGraph has " << dBG.getNumValidNodes()
                << " nodes and " <<  dBG.getNumValidArcs() << " arcs\n";
        }

        /*dBG.defrag();
        dBG.sanityCheck();

        // create a <kmer, nodePosPair> table
        KmerNPPTable table(dBG, true);
        dBG.createKmerNPPTable(table);
        
        writeTrueStrainMultiplicities(dBG, settings, table,
                                      trueNodeMult, trueEdgeMult,
                                      true);*/
        
        // write stage 2 binary output file
        if(write)
                dBG.writeBinary(to_string(nodeModel.getNumStrains()) + "strains." + settings.getStage2GraphFilename());


        double elapsed = Util::stopChrono();
        cout << "Stage 2 graph cleaning finished in " << elapsed << endl;
}

int stageTwoTestStrainNum(Settings& settings, DBGraph& dBG, 
                           int numStrains, 
                           CovModel& nodeModel, CovModel& edgeModel,
                           NodeMap<Multiplicity>& nodeMult,
                           EdgeMap<Multiplicity>& edgeMult,
                           bool checkGoodness = true)
{
        // get an initial (rough) k-mer model fit
        CovModel kmerModel = dBG.getInitCovModel(settings, numStrains);
        
        cout << "Initial (rough) coverage estimate [error: " << fixed
             << std::setprecision(2) << kmerModel.getErrLambda()
             << " (" << kmerModel.getErrorODF() << "); "
             << "sum-of-strains: " << kmerModel.getSosLambda()
             << " (" << kmerModel.getODF() << ")]\n";
//        cout << kmerModel << endl;

        // get the kmer histogram
        map<unsigned int, double> nodeCovHist;
        dBG.getNodeCovHist(nodeCovHist);
        uint xMin = max<int>(0, settings.getAbundanceMin());
        uint xMax = (kmerModel.getMaxMult() + 0.5) * kmerModel.getSosLambda();
        
        // write node histogram to file
        ofstream histOFS("nodeKmerHist.tsv");
        for (const auto it : nodeCovHist)
                histOFS << it.first << "\t" << it.second << "\n";
        histOFS.close();
        
        // fit the model the kmer histogram
        cout << "Fit model to k-mer histogram using EM..." << endl;
        
        int numIt; double relErr;
        double convergence = max(settings.getEMConvEps(), settings.getEMMaxChange());
        tie(numIt, relErr) = kmerModel.fitEM(nodeCovHist, xMin, xMax,
                                             convergence,
                                             settings.getEMMaxIter(), settings.useGenusLevelModel());


        cout << "\tEM algorithm " << (relErr > convergence ?
        "did NOT converge:" : "converged:") << "\n\tno. of iterations: "
        << numIt << "; rel. error: " << scientific << relErr << endl;
        
        cout << "Model after EM fitting: " << kmerModel << endl;
                
        
        kmerModel.writeGnuplot(to_string(numStrains) + "strains." + "haplotypeNBMix", nodeCovHist);
        
        // MODEL SCORES
        cout << "Init Model RMSE: " << sqrt(kmerModel.getMSE(nodeCovHist)) << endl;
        vector<double> nCov;
        nCov.reserve(dBG.getNumValidNodes());
        for (auto node: dBG.getNodeReps(dBG.getNumNodes()))
                nCov.push_back(dBG.getSSNode(node.getNodeID()).getAvgCov());
        double LL = kmerModel.getLL(nCov);
        cout << "Init Model LL: " << LL << "\n"
        << "Init Model fit entropy: " << kmerModel.getEntropy(nCov) << endl;
        
        cout    << "\t AIC: " << kmerModel.getAIC(LL) << "\n"
        << "\t BIC: " << kmerModel.getBIC(LL, nCov.size()) << "\n"
        << "\t CLC: " << kmerModel.getCLC(LL, kmerModel.getEntropy(nCov)) << "\n"
        << "\t ICL: " << kmerModel.getICL(LL, kmerModel.getEntropy(nCov), nCov.size()) << "\n"
        << "\t Entropy IC: " << setprecision(7) << kmerModel.entropyCriterion(kmerModel.getEntropy(nCov), nCov.size()) << endl;
        
        if( checkGoodness && kmerModel.getStrainLambda(0) < 1.0){
                cout << "Smallest strain coverage < 1, assuming less strains are present and aborting estimation." << endl;
                return -1;
        }
                
        // EM on mixture of NB for CRFMaxMult mults
        int newAbMin = settings.getAbundanceMin();
        int removed = 0;
        
        if (settings.correct()){
                stageTwoCorrectGraph(settings, kmerModel, kmerModel,
                                     dBG, newAbMin, removed);
                
                cout << "Number of nodes removed: " << removed << endl;
                cout << "New minimum abundance: " << newAbMin << endl;
                
                settings.updateAbundanceMin(newAbMin);
                //wN[0] *= (1.0 - ( (double) removed / (double) dBG.getNumNodes() ));
                //wE[0] *= (1.0 - ( (double) removed / (double) dBG.getNumNodes() ));
                //if (! settings.fixEMZeroMult())
                //        kmerModel.setWeight(vector<int>(numStrains, 0), 2.0);
                
                //wN[0] = 2.0;
                //wE[0] = 2.0;
        }
        
        
        CRFMult myCRFMult(dBG, settings.getCRFDepth(),
                          settings.getCRFMaxFactSize(),
                          settings.getCRFFlowStrength(),
                          settings.approxInf(),
                          settings.getNumThreads(),
                          settings.getThreadGraphWorkSize());
        
        vector<NodeRep> nodes;
        vector<EdgeRep> edges;
        
        if (settings.getEMTrainingSize() != -1) {
                nodes = dBG.getNodeReps(settings.getEMTrainingSize());
                edges = dBG.getEdgeReps(settings.getEMTrainingSize());
        } else { // train on full data
                nodes = dBG.getNodeReps(dBG.getNumNodes());
                edges = dBG.getEdgeReps(dBG.getNumArcs());
        }
        
        dBG.rescaleKmerToNodeCts(nodes,kmerModel);
        
        cout << kmerModel << endl;
        
        if (! settings.fixEMZeroMult() && settings.correct()) { 
                cout << "zero-multiplicity model will be re-estimated after cleaning, setting initial weight to 2.0" << endl;
                kmerModel.setZeroWeight(2.0); 
        } 
        
        edgeModel = kmerModel;
        nodeModel = kmerModel;
        
        Multiplicity initMult(vector<int>(numStrains, 0));
        populateNodeMult(nodeMult, nodes, numStrains);
        populateEdgeMult(edgeMult, edges, numStrains);
        
        cout.precision(2);
        cout << "Fitting mixture model using EM for " << nodes.size()
        << " nodes and " << edges.size() << " edges." << endl;
        
        /*cout << "Initial node model: " << endl;
         c *out << nodeModel << endl;
         cout << "Initial edge model" << endl;
         cout << edgeModel << endl;*/
        
        int nIters = myCRFMult.computeMultEM(nodeMult, nodeModel,
                                             edgeMult, edgeModel,
                                             settings.getEMConvEps(),
                                             settings.getEMMaxChange(),
                                             settings.getEMMaxIter(),
                                             settings.useEMEqualWeights(),
                                             settings.approxInf(),
                                             settings.computeMAP(), settings.fixEMZeroMult());
        
        // write fitted coverage models
        nodeModel.write(to_string(numStrains) + "strains." + settings.getNodeModelFilename());
        edgeModel.write(to_string(numStrains) + "strains." + settings.getEdgeModelFilename());
        
        if (nIters <= settings.getEMMaxIter())
                cout << "EM algorithm converged after "
                << nIters << " iterations\n";
        else
                cout << "WARNING: maximum number of iterations reached. "
                << "Convergence is not guaranteed." << endl;
        return 0;
}


void stageTwo(Settings& settings)
{
        cout << "\nEntering stage 2\n";
        cout << "================\n" << endl;

        if (!settings.stageTwoNecessary()) {
                cout << "Files " << settings.getNodeModelFilename() << " and "
                     << settings.getEdgeModelFilename()
                     << " exist. Skipping stage 2...\n";
		return;
        }
        
        Util::startChrono();
        
        int ns = settings.getNStrains();
        // If number of strains is not passed: test for 2, 3 or 4 strains present
        if (ns < 0){
        
                int minStrainN = settings.getInitNStrains();
                int numStrainTests = 3;
                
                double initAbmin = settings.getAbundanceMin();
                
                vector<pair<CovModel,CovModel>> covmodels;
                covmodels.reserve(3);
                DBGraph dBGInit(settings);
                dBGInit.loadBinary(settings.getStage1GraphFilename());
                vector<NodeRep> remainingNodes(dBGInit.getNodeReps(dBGInit.getNumValidNodes()));
                vector<EdgeRep> remainingEdges(dBGInit.getEdgeReps(dBGInit.getNumValidArcs()));
                
                vector<NodeMap<Multiplicity>> estNodeMults;
                vector<EdgeMap<Multiplicity>> estEdgeMults;
                
                vector<double> BICscores(numStrainTests,0.0);
                
                // Train all models
                for (int nstrains = minStrainN; nstrains < (minStrainN+numStrainTests); nstrains ++){
                        // Load stage 1 de Bruijn graph
                        DBGraph dBG(settings);
                        dBG.loadBinary(settings.getStage1GraphFilename());
                        cout << "Loaded graph with " << dBG.getNumValidNodes() << " nodes and "
                        << dBG.getNumValidArcs() << " edges" << endl;

                        CovModel nodeModel, edgeModel;
                        vector<NodeRep> remaining;
                        vector<EdgeRep> remainingE;
                        NodeMap<Multiplicity> nodeMult;
                        EdgeMap<Multiplicity> edgeMult;
                        
                        cout << "Cleaning graph and building model under the assumption of " << nstrains << " strains..."  << endl;
                        Util::startChrono();
                        
                        int goodModel = stageTwoTestStrainNum(settings, dBG, 
                                        nstrains, 
                                        nodeModel, edgeModel,
                                        nodeMult, edgeMult);
                        
                        remaining = dBG.getNodeReps(dBG.getNumValidNodes());
                        remainingE = dBG.getEdgeReps(dBG.getNumValidArcs());
                        
                        // If smallest lambda is < 1 during model initialisation, we assume less strains are present than we are testing now
                        if (goodModel < 0){
                                for (int i=nstrains-minStrainN; i < numStrainTests; i++)
                                        BICscores[i] = goodModel;
                                break;
                        }
                        
                        // Keep track of all nodes/arcs in the current cleaned graph, for fair model score comparison 
                        vector<NodeRep> temp(remainingNodes);
                        remainingNodes.clear();
                        std::sort(temp.begin(), temp.end());
                        std::sort(remaining.begin(), remaining.end());
                        std::set_intersection(temp.begin(), temp.end(),
                                        remaining.begin(), remaining.end(),
                                        std::back_inserter(remainingNodes));
                        vector<EdgeRep> tempE(remainingEdges);
                        remainingEdges.clear();
                        std::sort(tempE.begin(), tempE.end());
                        std::sort(remainingE.begin(), remainingE.end());
                        std::set_intersection(tempE.begin(), tempE.end(),
                                        remainingE.begin(), remainingE.end(),
                                        std::back_inserter(remainingEdges));
                        covmodels.push_back(pair<CovModel,CovModel>(nodeModel, edgeModel));
                        estNodeMults.push_back(nodeMult);
                        estEdgeMults.push_back(edgeMult);
                        settings.updateAbundanceMin(initAbmin);
                        
                        cout << "... finished " << nstrains << " strains model (" << Util::stopChronoStr() << ")\n------------------------------------------\n" << endl;
                }
                
                vector<NodeRep> nodes(remainingNodes);
                vector<EdgeRep> edges(remainingEdges);
#ifdef DEBUG               
                ofstream remNodeOut("modelscorenodes.tsv");
                ofstream remEdgeOut("modelscoreedges.tsv");
#endif
                
                // Compute model scores based on all nodes and arcs present in all 3 corrected graphs
                // And write out models for all strain assumptions
                for (int nstrains = minStrainN; nstrains < (minStrainN+numStrainTests); nstrains ++){
                        if(BICscores[nstrains-minStrainN] < 0){
                                BICscores[nstrains-minStrainN] = std::numeric_limits<double>::max();
                                break;
                        }
                        DBGraph dBG(settings);
                        dBG.loadBinary(settings.correct() ? (to_string(nstrains) + "strains." + settings.getStage2GraphFilename()) : settings.getStage1GraphFilename());
                        CovModel nodeModel = covmodels[nstrains-2].first;
                        CovModel edgeModel = covmodels[nstrains-2].second;
                        
                        vector<double> nodeCov, edgeCov;
                        nodeCov.reserve(nodes.size());
                        edgeCov.reserve(edges.size());
                        for (auto node: nodes){
                                nodeCov.push_back(dBG.getSSNode(node.getNodeID()).getAvgCov());
#ifdef DEBUG                                
                                remNodeOut << node.getNodeID() << "\t" << nstrains << "\t" << dBG.getSSNode((node.getNodeID())).getMarginalLength() << "\t" << dBG.getSSNode(node.getNodeID()).getAvgCov() << "\n";
#endif
                        }
                        for (auto edge: edges){
                                edgeCov.push_back(dBG.getArc(dBG.getArcID(edge)).getCov());
#ifdef DEBUG                                
                                remEdgeOut << edge.getSrcID() << "\t" << edge.getDstID() << "\t" << nstrains << "\t" << dBG.getArc(dBG.getArcID(edge)).getCov() << "\n";
#endif
                        }
                
                        
                        map<unsigned int, double> nodeHist;
                        for (const auto& it : estNodeMults[nstrains-2]) {
                                SSNode node = dBG.getSSNode(it.first);
                                double f = node.getAvgCov() - floor(node.getAvgCov());
                                nodeHist[node.getAvgCov()] += (1.0-f);// * node.getMarginalLength();
                                nodeHist[node.getAvgCov() + 1] += f;// * node.getMarginalLength();
                        }
                        
                        cout << "\n" << nstrains << " strains model scores:" << endl;
                        
                        nodeModel.writeGnuplot(to_string(nstrains) + "strains."+ "nodes", nodeHist);
                        cout << "Node Model RMSE: " << sqrt(nodeModel.getMSE(nodeHist)) << endl;
                        double nodeLL = nodeModel.getLL(nodeCov);
                        double crfEntropy = nodeModel.getCRFassignmentEntropy(estNodeMults[nstrains-2]);
                        double singletonEntropy = nodeModel.getEntropy(nodeCov);
                        cout << "Node Model LL: " << nodeLL  << "\n"
                        << "CRF assignment entropy: " << crfEntropy
                        << "\nModel fit entropy: " << singletonEntropy
                        << "\n\t AIC: " << nodeModel.getAIC(nodeLL) <<"\n"
                        << "\t BIC: " << nodeModel.getBIC(nodeLL, nodeCov.size()) << "\n"
                        << "\t CLC: " << nodeModel.getCLC(nodeLL, nodeModel.getEntropy(nodeCov)) << "\n" 
                        << "\t CLC (crf entropy): " << "\t" << nodeModel.getCLC(nodeLL, nodeModel.getCRFassignmentEntropy(estNodeMults[nstrains-2])) << "\n"
                        << "\t ICL: " << nodeModel.getICL(nodeLL, nodeModel.getEntropy(nodeCov), nodeCov.size()) << "\n"
                        << "\t ICL (crf entropy): " <<  nodeModel.getICL(nodeLL, nodeModel.getCRFassignmentEntropy(estNodeMults[nstrains-2]), estNodeMults[nstrains-2].size()) << "\n"
                        << "\t Entropy IC: " << setprecision(7) << nodeModel.entropyCriterion(nodeModel.getEntropy(nodeCov), nodeCov.size())
                        << "\t CRF Entropy IC: " << nodeModel.entropyCriterion(nodeModel.getCRFassignmentEntropy(estNodeMults[nstrains-2]), estNodeMults[nstrains-2].size()) << endl;
                        
                        map<unsigned int, double> edgeHist;
                        for (const auto& it : estEdgeMults[nstrains-2]) {
                                Arc& arc = dBG.getArc(dBG.getArcID(it.first));
                                double f = arc.getCov() - floor(arc.getCov());
                                edgeHist[arc.getCov()] += (1.0-f);
                                edgeHist[arc.getCov() + 1] += f;
                        }
                        
                        edgeModel.writeGnuplot(to_string(nstrains) + "strains."+ "edges", edgeHist);
                        cout << "Edge Model RMSE: " << sqrt(edgeModel.getMSE(edgeHist)) << endl;
                        double edgeLL = edgeModel.getLL(edgeCov);
                        crfEntropy = edgeModel.getCRFassignmentEntropy(estEdgeMults[nstrains-2]);
                        singletonEntropy = edgeModel.getEntropy(edgeCov);
                        cout << "Edge Model LL: " << edgeLL << "\n"
                        << "CRF assignment entropy: " << crfEntropy
                        << "\nModel fit entropy: " << singletonEntropy
                        << "\n\t AIC: " << edgeModel.getAIC(edgeLL) << "\n"
                        << "\t BIC: " << edgeModel.getBIC(edgeLL, edgeCov.size()) << "\n"
                        << "\t CLC: " << edgeModel.getCLC(edgeLL, edgeModel.getEntropy(edgeCov)) << "\n"
                        << "\t CLC (crf entropy): " <<  edgeModel.getCLC(edgeLL, edgeModel.getCRFassignmentEntropy(estEdgeMults[nstrains-2])) << "\n"
                        << "\t ICL: " << edgeModel.getICL(edgeLL, edgeModel.getEntropy(edgeCov), edgeCov.size()) << "\n" 
                        << "\t ICL (crf entropy): " <<  edgeModel.getICL(edgeLL, edgeModel.getCRFassignmentEntropy(estEdgeMults[nstrains-2]), estEdgeMults[nstrains-2].size()) << "\n"
                        << "\t Entropy IC: " << setprecision(7) << edgeModel.entropyCriterion(edgeModel.getEntropy(edgeCov), edgeCov.size())
                        << "\t CRF Entropy IC: " << edgeModel.entropyCriterion(edgeModel.getCRFassignmentEntropy(estEdgeMults[nstrains-2]), estEdgeMults[nstrains-2].size()) << endl;
                        
                        // We will be using the BIC score for the node model to select the number of strains
                        BICscores[nstrains-minStrainN] = nodeModel.getBIC(nodeLL, nodeCov.size());
                }
                
                int optNStrains = minStrainN +  std::distance(std::begin(BICscores), std::min_element(std::begin(BICscores), std::end(BICscores)));
                
                cout << "Choosing model with " << optNStrains << " strains." << endl;
                
                settings.setNStrains(optNStrains);
                
                // write fitted coverage models
                covmodels[optNStrains-minStrainN].first.write(settings.getNodeModelFilename());
                covmodels[optNStrains-minStrainN].second.write(settings.getEdgeModelFilename());
                
                DBGraph dBGFinal(settings);
                if(settings.correct()){
                        dBGFinal.loadBinary(to_string(optNStrains) + "strains." + settings.getStage2GraphFilename());
                
                        dBGFinal.defrag();
                        dBGFinal.sanityCheck();
                        
                        // create a <kmer, nodePosPair> table
                        KmerNPPTable table(dBGFinal, true);
                        dBGFinal.createKmerNPPTable(table);
                        
                        NodeMap<vector<int>> trueNodeMult;
                        EdgeMap<vector<int>> trueEdgeMult;
                        writeTrueStrainMultiplicities(dBGFinal, settings, table,
                                                trueNodeMult, trueEdgeMult, true);
                        
                        // write stage 2 binary output file
                        dBGFinal.writeBinary(settings.getStage2GraphFilename());
                }
        } else {
        // The number of strains was passed
                
                // load stage 1 de Bruijn graph
                DBGraph dBG(settings);
                dBG.loadBinary(settings.getStage1GraphFilename());
                cout << "Loaded graph with " << dBG.getNumValidNodes() << " nodes and "
                << dBG.getNumValidArcs() << " edges" << endl;
                
                CovModel nodeModel, edgeModel;
                NodeMap<Multiplicity> nodeMult;
                EdgeMap<Multiplicity> edgeMult;
                
                stageTwoTestStrainNum(settings, dBG, ns, 
                                      nodeModel, edgeModel, 
                                      nodeMult, edgeMult, false);
                
                
                map<unsigned int, double> nodeHist;
                for (const auto& it : nodeMult) {
                        SSNode node = dBG.getSSNode(it.first);
                        double f = node.getAvgCov() - floor(node.getAvgCov());
                        nodeHist[node.getAvgCov()] += (1.0-f);// * node.getMarginalLength();
                        nodeHist[node.getAvgCov() + 1] += f;// * node.getMarginalLength();
                }
                
                nodeModel.writeGnuplot("nodes", nodeHist);
                
                map<unsigned int, double> edgeHist;
                for (const auto& it : edgeMult) {
                        Arc& arc = dBG.getArc(dBG.getArcID(it.first));
                        double f = arc.getCov() - floor(arc.getCov());
                        edgeHist[arc.getCov()] += (1.0-f);
                        edgeHist[arc.getCov() + 1] += f;
                }
                
                edgeModel.writeGnuplot("edges", edgeHist);
                
                dBG.defrag();
                dBG.sanityCheck();
                
                // create a <kmer, nodePosPair> table
                KmerNPPTable table(dBG, true);
                dBG.createKmerNPPTable(table);
                
                NodeMap<vector<int>> trueNodeMult;
                EdgeMap<vector<int>> trueEdgeMult;
                writeTrueStrainMultiplicities(dBG, settings, table,
                                              trueNodeMult, trueEdgeMult, true);
                
               
                // write fitted coverage models
                nodeModel.write(settings.getNodeModelFilename());
                edgeModel.write(settings.getEdgeModelFilename());
                
                // write stage 2 binary output file
                dBG.writeBinary(settings.getStage2GraphFilename());
        }
        
        cout << "Stage 2 finished in " << Util::stopChronoStr() << endl;
}


void stageThreeCompMult(Settings& settings)
{
        cout << "\nEntering stage 3\n";
        cout << "================\n" << endl;
        
        Util::startChrono();
        
        DBGraph dBG(settings);
        if (settings.correct()){
                if (! Util::fileExists(settings.getStage2GraphFilename()))
                        throw runtime_error("File " + settings.getStage2GraphFilename() + " cannot be found, aborting." );
                dBG.loadBinary(settings.getStage2GraphFilename());
        } else {
                if (! Util::fileExists(settings.getStage1GraphFilename()))
                        throw runtime_error("File " + settings.getStage1GraphFilename() + " cannot be found, aborting." );
                dBG.loadBinary(settings.getStage1GraphFilename());
        }
        // create a <kmer, nodePosPair> table
        KmerNPPTable table(dBG, true);
        dBG.createKmerNPPTable(table);
                
        NodeMap<vector<int>> trueNodeMult;
        EdgeMap<vector<int>> trueEdgeMult;
        writeTrueStrainMultiplicities(dBG, settings, table,
                                      trueNodeMult, trueEdgeMult, true);
 
        cout << "Loaded graph with " << dBG.getNumValidNodes() << " nodes and "
        << dBG.getNumValidArcs() << " edges." << endl;
        
        if (! Util::fileExists(settings.getNodeModelFilename()) || ! Util::fileExists(settings.getEdgeModelFilename()))
                throw runtime_error("Coverage model files not found, aborting" );
        
        CovModel nodeModel(settings.getNodeModelFilename());
        CovModel edgeModel(settings.getEdgeModelFilename());
        
        cout << "Node model: " << nodeModel << endl;
        cout << "Edge model: " << edgeModel << endl;
        
        CRFMult myCRFMult(dBG, settings.getCRFDepth(),
                          settings.getCRFMaxFactSize(),
                          settings.getCRFFlowStrength(),
                          settings.approxInf(),
                          settings.getNumThreads(),
                          settings.getThreadGraphWorkSize());
        
        NodeMap<Multiplicity> nodeMult(dBG.getNumValidNodes());
        EdgeMap<Multiplicity> edgeMult(dBG.getNumValidArcs());
        vector<NodeRep> nodes = dBG.getNodeReps(dBG.getNumValidNodes());
        vector<EdgeRep> edges = dBG.getEdgeReps(dBG.getNumValidArcs());
        populateNodeMult(nodeMult, nodes, nodeModel.getNumStrains());
        populateEdgeMult(edgeMult, edges, edgeModel.getNumStrains());
        
        double logZ=-1.0;
        
        if(settings.approxInf())
                myCRFMult.approxMultAll(nodeMult, nodeModel, edgeMult, edgeModel, logZ,
                                        settings.computeMAP(), settings.useSingleCRF());
        else
                myCRFMult.computeMult(nodeMult, nodeModel, edgeMult, edgeModel);
        
        
        ofstream nodeOFS(settings.getEstNodeMultFilename());
        ofstream edgeOFS(settings.getEstEdgeMultFilename());
        
        for (auto& it: nodeMult) {
                SSNode node = dBG.getSSNode(it.first);
                nodeOFS << it.first.getNodeID() << "\t";
                //int s = 0;
                for ( int sMult : it.second.getExpMult()){
                        nodeOFS << sMult << "\t";
                }
                nodeOFS << it.second.getExpMultLProb() << "\t"
                << node.getAvgCov() << "\t"
                << node.getMarginalLength() << "\t"
                << it.second.getEntropy() << "\n";
        }
        
        for (auto& it: edgeMult) {
                edgeOFS << it.first.getSrcID() << "\t"
                << it.first.getDstID() << "\t";
                for (int sMult: it.second.getExpMult())
                        edgeOFS << sMult  << "\t";
                edgeOFS << it.second.getExpMultLProb() << "\t"
                << dBG.getArc(it.first).getCov() << "\t"
                << it.second.getEntropy() << "\n";
        }
        
        cout << "Wrote file: " << settings.getEstNodeMultFilename() << endl;
        cout << "Wrote file: " << settings.getEstEdgeMultFilename() << endl;
        
        cout << "Stage 3 finished in " << Util::stopChronoStr() << endl;
}



void stageThreeCG(Settings& settings)
{
        cout << "\nEntering stage 3\n";
        cout << "================\n" << endl;

        DBGraph dBG(settings);
        if (settings.correct()){
                if (! Util::fileExists(settings.getStage2GraphFilename()))
                        throw runtime_error("File " + settings.getStage2GraphFilename() + " cannot be found, aborting." );
                dBG.loadBinary(settings.getStage2GraphFilename());
        } else {
                if (! Util::fileExists(settings.getStage1GraphFilename()))
                        throw runtime_error("File " + settings.getStage1GraphFilename() + " cannot be found, aborting." );
                dBG.loadBinary(settings.getStage1GraphFilename());
        }
        //int nStrains = settings.getInitNStrains();

        NodeID noi = settings.getVisGraphNode();
        if ((noi < -dBG.getNumNodes()) || (noi > dBG.getNumNodes()))
                throw runtime_error("Specified node id of " + to_string(noi) +
                                    " does not exist in the graph");

        if (! Util::fileExists(settings.getNodeModelFilename()) || ! Util::fileExists(settings.getEdgeModelFilename()))
                throw runtime_error("Coverage model files not found, aborting" );        
        
        CovModel nodeModel(settings.getNodeModelFilename());
        CovModel edgeModel(settings.getEdgeModelFilename());
        
        int nStrains = nodeModel.getNumStrains();

        cout << "Node model: " << nodeModel << endl;
        cout << "Edge model: " << edgeModel << endl;
        
        vector<NodeID> nodes;
        vector<pair<NodeID, NodeID> > edges;
        
	dBG.getSubgraph(noi, nodes, edges, settings.getCRFDepth());

        vector<NodeRep> nodeReps(nodes.begin(), nodes.end());
        vector<EdgeRep> edgeReps(edges.begin(), edges.end());
        
        vector<Multiplicity> nodeMult(nodes.size()), edgeMult(edges.size());
        NodeMap<Multiplicity> estNodeMult(nodes.size()); EdgeMap<Multiplicity> estEdgeMult(edges.size());
        
        if (settings.approxInf()){
                CRFSolver myCRFSolver(dBG, settings.getCRFMaxFactSize(),
				      settings.getCRFFlowStrength());
		myCRFSolver.approxSubgraphMult(noi, estNodeMult, estEdgeMult,
                                               nodeModel, edgeModel, settings.getCRFDepth(),
                                               settings.computeMAP(), settings.getNumThreads(),
                                               settings.getThreadGraphWorkSize());
	}else{
                CRFMult myCRFMult(dBG, settings.getCRFDepth(),
                          settings.getCRFMaxFactSize(),
                          settings.getCRFFlowStrength(),
                          settings.approxInf(),
                          settings.getNumThreads(),
                          settings.getThreadGraphWorkSize());

                myCRFMult.computeMult(estNodeMult, nodeModel,
                                estEdgeMult, edgeModel);
        }
        
        // If they exist, also load the true multiplicities from disk
        // Note: for efficiency, story only true multiplicities that are needed
        NodeMap<vector<int>> trueNodeMult;
        transform(nodeReps.begin(), nodeReps.end(),
                  inserter(trueNodeMult, trueNodeMult.end()),
                  [nStrains](const NodeRep& nr) { return make_pair(nr, vector<int>(nStrains,-1)); });
        
        EdgeMap<vector<int>> trueEdgeMult;
        transform(edgeReps.begin(), edgeReps.end(),
                  inserter(trueEdgeMult, trueEdgeMult.end()),
                  [nStrains](const EdgeRep& er) { return make_pair(er, vector<int>(nStrains,-1)); });
        
        readTrueStrainMultiplicities(settings,
                                     nodeReps, edgeReps,
                                     trueNodeMult, trueEdgeMult,
                                     settings.correct());
        
        dBG.writeCytoscapeGraph("cytgraph" + to_string(noi) + "nb" + to_string(settings.getCRFDepth()),
                                nodes, edges,
                                estNodeMult, estEdgeMult,
                                trueNodeMult, trueEdgeMult);
        
}

/*TODO reinstate??
 * void drawFullCytoGraph(Settings& settings, DBGraph& dBG) 
{
        vector<NodeID> nID; vector<EdgeID> eID;
        dBG.getGraph(nID, eID);
        
        vector<NodeRep> nodes;
        vector<EdgeRep> edges;
        NodeMap<Multiplicity> estNodeMult;
        EdgeMap<Multiplicity> estEdgeMult;
        
        ifstream ifs(settings.getEstNodeMultFilename());
        if (!ifs){
                cerr << "Could not find estimated nodes multiplicities file...\n"
                "Cannot write cytoscape graph\n";
                return;
        }
        
        while (ifs) {
                NodeID nodeID; int m; double lp; double d1; size_t d2;
                ifs >> nodeID >> m >> lp >> d1 >> d2;
                if (!ifs)
                        break;
                nodes.push_back(NodeRep(nodeID));
                estNodeMult[NodeRep(nodeID)] = Multiplicity(m, {lp,log(1.0-exp(lp))});
        }
        ifs.close();
        
        ifs.open(settings.getEstEdgeMultFilename());
        if (!ifs){
                cerr << "Could not find estimated edge multiplicities file...\n"
                "Cannot write cytoscape graph\n";
                return;
        }
        
        while (ifs) {
                NodeID srcID; NodeID dstID; int m; double lp; double d1;
                ifs >> srcID >> dstID >> m >> lp >> d1;
                if (!ifs)
                        break;
                edges.push_back(EdgeRep(srcID, dstID));
                estEdgeMult[EdgeRep(srcID, dstID)] = Multiplicity(m, {lp, log(1.0-exp(lp))});
        }
        ifs.close();
        
        // If they exist, also load the true multiplicities from disk
        // Note: for efficiency, story only true multiplicities that are needed
        NodeMap<vector<int>> trueNodeMult;
        transform(nodeReps.begin(), nodeReps.end(),
                inserter(trueNodeMult, trueNodeMult.end()),
                [nStrains](const NodeRep& nr) { return make_pair(nr, vector<int>(nStrains,-1)); });

        EdgeMap<vector<int>> trueEdgeMult;
        transform(edgeReps.begin(), edgeReps.end(),
                inserter(trueEdgeMult, trueEdgeMult.end()),
                [nStrains](const EdgeRep& er) { return make_pair(er, vector<int>(nStrains,-1)); });

        ifstream ifs(settings.getTrueNodeStrainFilename());

        if (!ifs)
            cerr << "Could not find true nodes multiplicities file...\n"
                    "True node strain multiplicities will be set to -1 in "
                    "resulting Cytoscape graph\n";

        while (ifs) {
                NodeID nodeID; int m; vector<int> strainM;
                ifs >> nodeID;
                for (int i=0; i < nStrains; i++){
                ifs >> m;
                strainM.push_back(m);
                }
                if (!ifs)
                        break;
                auto it = trueNodeMult.find(nodeID);
                if (it != trueNodeMult.end())
                        it->second = strainM;
        }
        ifs.close();

        ifs.open(settings.getTrueEdgeStrainFilename());
        if (!ifs)
                cerr << "Could not find true edges multiplicities file...\n"
                    "True edge multiplicities will be set to -1 in "
                    "resulting Cytoscape graph\n";

        while (ifs) {
                NodeID srcID, dstID; int m; vector<int> strainM;
                ifs >> srcID >> dstID;
                for (int i=0; i < nStrains; i++){
                    ifs >> m;
                    strainM.push_back(m);
                }
                if (!ifs)
                        break;
                auto it = trueEdgeMult.find(EdgeRep(srcID, dstID));
                if (it != trueEdgeMult.end())
                        it->second = strainM;
        }
        ifs.close();
        
        cout << "Writing Cytoscape graph... " << endl;
        dBG.writeCytoscapeGraph("Cytograph.full",
                                nID, eID,
                                estNodeMult, estEdgeMult,
                                trueNodeMult, trueEdgeMult);
}*/


void stageFourPhasedGraphs(Settings& settings)
{
        cout << "\nEntering stage 4\n";
        cout << "================\n" << endl;

        Util::startChrono();
        
        CovModel nodeModel(settings.getNodeModelFilename());
        CovModel edgeModel(settings.getEdgeModelFilename());
        
        assert(nodeModel.getNumStrains() == edgeModel.getNumStrains());
        
        int nStrains = nodeModel.getNumStrains();

        cout << "\nGetting haplotype aware graphs for " << nStrains << " strains...\n" << endl;

        ifstream nodeIFS(settings.getEstNodeMultFilename());
        ifstream edgeIFS(settings.getEstEdgeMultFilename());

        std::map<NodeRep, std::vector<int>> nodeMult;
        std::map<EdgeRep, std::vector<int>> edgeMult;

	NodeID nodeID;
        while (nodeIFS >> nodeID) {
                int m; vector<int> strainM; double temp;
                for (int i=0; i < nStrains; i++){
                        nodeIFS >> m;
                        strainM.push_back(m);
                }
                nodeIFS >> temp;
                nodeIFS >> temp;
                nodeIFS >> temp;
                nodeIFS >> temp;

                nodeMult[NodeRep(nodeID)] = strainM;
        }

	NodeID srcID;
        while (edgeIFS >> srcID) {
                NodeID dstID;
                int m;
                vector<int> strainM;
                double temp;
                edgeIFS >> dstID;
                for (int i = 0; i < nStrains; i++) {
                        edgeIFS >> m;
                        strainM.push_back(m);
                }
                edgeIFS >> temp;
                edgeIFS >> temp;
                edgeIFS >> temp;

                edgeMult[EdgeRep(srcID, dstID)] = strainM;
        }

        for(int i = 1; i <= nStrains; i++){
                DBGraph dBG(settings);
                if (settings.correct())
                        dBG.loadBinary(settings.getStage2GraphFilename());
                else
                        dBG.loadBinary(settings.getStage1GraphFilename());

                vector<NodeRep> nodeReps = dBG.getNodeReps(dBG.getNumValidNodes());
                vector<EdgeRep> edgeReps = dBG.getEdgeReps(dBG.getNumValidArcs());
                assert(nodeReps.size() == nodeMult.size());
                assert(edgeReps.size() == edgeMult.size());

                vector<NodeRep> toRemove;
                for (size_t j = 0 ; j < nodeReps.size(); j++){
                        if (nodeMult[nodeReps[j]][i-1] == 0)
                                toRemove.push_back(nodeReps[j]);
                }
                dBG.removeNodes(toRemove);

                vector<EdgeRep> edgesToRemove;
                for (size_t j = 0 ; j < edgeReps.size(); j++){
                        if (edgeMult[edgeReps[j]][i-1] == 0)
                                edgesToRemove.push_back(edgeReps[j]);
                }

                dBG.removeEdges(edgesToRemove);
                dBG.concatenateNodes();
                dBG.defrag();

                cout << "Strain " << i << "'s graph has " << dBG.getNumValidNodes() << " nodes and "
                <<  dBG.getNumValidArcs() << " arcs\nWriting graph..." << endl;
                dBG.writeContigs("strain" + to_string(i) + ".phased.fasta");

        }
        cout << "Stage 4 finished in " << Util::stopChronoStr() << endl;
}


int main(int argc, char** argv)
{
        try {
                Settings settings(argc, argv);
                LibraryContainer libraries(settings.getReadFilename());

                Multiplicity::maxM = settings.getCRFMaxMult();


                cout << "Welcome to HaploDetox v" << HAPLODETOX_MAJOR_VERSION << "."
                     << HAPLODETOX_MINOR_VERSION << "." << HAPLODETOX_PATCH_LEVEL;
		cout << " (commit " << GIT_SHA1 << ")";

#ifdef DEBUG
                cout << " (debug mode)" << endl;
#else
                cout << " (release mode)" << endl;
#endif
                
                stageOne(settings, libraries);
                stageTwo(settings);

                if (settings.getVisGraphNode() != 0) {
                        stageThreeCG(settings);
                } else {
                        stageThreeCompMult(settings);
                        stageFourPhasedGraphs(settings);
                }


        } catch (exception& e) {
                cerr << "Fatal error: " << e.what() << endl;
                return EXIT_FAILURE;
        }

        cout << "Exiting... bye!" << endl;
        return EXIT_SUCCESS;
}
