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

#include <numeric>
#include <thread>
#include <iomanip>

#include "crfmult.h"
#include "dbgraph.h"
#include "pgm/pgminference.h"
#include "util.h"
#include "libdai/dai/alldai.h"  // Include main libDAI header file
#include "libdai/dai/jtree.h"
#include "libdai/dai/bp.h"
#include "libdai/dai/parallelbp.h"
#include "libdai/dai/factor.h"
#include "libdai/dai/factorgraph.h"
#include "libdai/dai/var.h"
#include "libdai/dai/varset.h"
#include "ssnode.h"
#include "arc.h"

using namespace std;

// ============================================================================
// CONDITIONAL RANDOM FIELDS SOLVER (PER THREAD)
// ============================================================================

void CRFSolver::getSubgraph(int maxDepth)
{
        while (!pq.empty()) {
                // get and erase the current node
                NodeDFS currTop = pq.top();
                pq.pop();
                NodeID thisID = currTop.nodeID;
                int thisDepth = currTop.depth;

                // if the node was already handled, skip
                if (bitvec[abs(thisID)])
                        continue;

                nodes.push_back(NodeRep(thisID));
                SSNode n = dBG.getSSNode(thisID);

                // process the right arcs
                for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++) {
                        NodeID rightID = it->getNodeID();
                        if (bitvec[abs(rightID)])       // edge already added?
                                continue;

                        edges.push_back(EdgeRep(thisID, rightID));
                        if (thisDepth < maxDepth)
                                pq.push(NodeDFS(rightID, thisDepth+1));
                }

                // process the left arcs
                for (ArcIt it = n.leftBegin(); it != n.leftEnd(); it++) {
                        NodeID leftID = it->getNodeID();
                        if (bitvec[abs(leftID)])        // edge already added?
                                continue;
                        if (leftID == thisID)   // edge already added as right arc
                                continue;

                        edges.push_back(EdgeRep(leftID, thisID));
                        if (thisDepth < maxDepth)
                                pq.push(NodeDFS(leftID, thisDepth + 1));
                }

                // mark this node as handled
                bitvec[abs(thisID)] = true;
        }

        // reset all flags to false
        for (auto it : nodes)
                bitvec[it.getNodeID()] = false;
}

void CRFSolver::getSubgraph(NodeRep seedNode, int maxDepth)
{
        nodes.clear(); edges.clear();

        if (maxDepth == 0) {
                nodes.push_back(seedNode);
                return;
        }

        // add the seed node to the priority queue
        pq.push(NodeDFS(seedNode, 0));

        getSubgraph(maxDepth);
}

void CRFSolver::getSubgraph(EdgeRep seedEdge, int maxDepth)
{
        nodes.clear(); edges.clear();

        // special case (maxDepth == 0)
        if (maxDepth == 0) {
                edges.push_back(seedEdge);
                return;
        }

        // other cases (maxDepth > 0)
        pq.push(NodeDFS(seedEdge.getSrcID(), 1));
        pq.push(NodeDFS(seedEdge.getDstID(), 1));

        getSubgraph(maxDepth);
}

vector<NodeID> CRFSolver::getFullConnectedComponent(NodeRep seedNode)
{
        nodes.clear(); edges.clear();
        vector<NodeID> usedIDs;
        usedIDs.clear();
        
        // add the seed node to the priority queue
        pq.push(NodeDFS(seedNode, 0));
        
        // Same functionallity as getSubgraph but no maxDepth check
        while (!pq.empty()) {
                // get and erase the current node
                NodeDFS currTop = pq.top();
                pq.pop();
                NodeID thisID = currTop.nodeID;
                int thisDepth = currTop.depth;
                
                // if the node was already handled, skip
                if (bitvec[abs(thisID)])
                        continue;
                
                nodes.push_back(NodeRep(thisID));
                usedIDs.push_back(thisID);
                SSNode n = dBG.getSSNode(thisID);
                
                // process the right arcs
                for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++) {
                        NodeID rightID = it->getNodeID();
                        if (bitvec[abs(rightID)])       // edge already added?
                                continue;
                        
                        edges.push_back(EdgeRep(thisID, rightID));
                        pq.push(NodeDFS(rightID, thisDepth+1));
                }
                
                // process the left arcs
                for (ArcIt it = n.leftBegin(); it != n.leftEnd(); it++) {
                        NodeID leftID = it->getNodeID();
                        if (bitvec[abs(leftID)])        // edge already added?
                                continue;
                        if (leftID == thisID)   // edge already added as right arc
                                continue;
                        
                        edges.push_back(EdgeRep(leftID, thisID));
                        pq.push(NodeDFS(leftID, thisDepth + 1));
                }
                
                // mark this node as handled
                bitvec[abs(thisID)] = true;
        }
        
        // reset all flags to false
        for (auto it : nodes)
                bitvec[it.getNodeID()] = false;
        
        return usedIDs;
}

Factor CRFSolver::createSingletonFactor(vector<int> varIDs, vector<int> varCards,
                                        const CovModel& covModel,
                                        const double numObs, const int numIndepObs) //for 1 CRF variable per strain
{
        size_t numVal = accumulate(varCards.begin(), varCards.end(), 1ull, multiplies<size_t>());
        vector<double> val(numVal);
        int maxMult = Multiplicity::maxM;
        
        Assignment assignment(varCards);
        if (numObs > maxMult*covModel.getSosLambda()){
                for (size_t i = 0; i < numVal; i++) {
                        val[i] = (find(assignment.begin(), assignment.end(),maxMult) == assignment.end())? covModel.getLogProb(numObs, assignment): covModel.getMaxMultP();
                        assignment.increment();
                }
        } else {
                for (size_t i = 0; i < numVal; i++) {
                        val[i] = covModel.getLogProb(numObs, assignment); //*numIndepObs; // TODO is this valuable here??? 
                        assignment.increment();
                }
        }
        return Factor(move(varIDs), move(varCards), move(val));
}

dai::Factor CRFSolver::createLDUniformSingletonFactor(vector<int> varIDs, 
                                                      vector<int> varCards)
{
        dai::VarSet vars;
        for(int i = 0; i < varIDs.size(); i++)
                vars.insert(dai::Var(varIDs[i], varCards[i]));
        
        size_t numVal = accumulate(varCards.begin(), varCards.end(), 1ull, multiplies<size_t>());
        vector<double> val(numVal);
        
        Assignment assignment(varCards);
        for (size_t i = 0; i < numVal; i++) {
                val[i] = 1.0;
                assignment.increment();
        }
        
        return dai::Factor(move(vars), move(val));
}


dai::Factor CRFSolver::createLDSingletonFactor(vector<int> varIDs, vector<int> varCards,
                                               const CovModel& covModel,
                                               const double numObs, const int numIndepObs)
{
        dai::VarSet vars;
        for(int i = 0; i < varIDs.size(); i++)
                vars.insert(dai::Var(varIDs[i], varCards[i]));
        
        size_t numVal = accumulate(varCards.begin(), varCards.end(), 1ull, multiplies<size_t>());
        vector<double> val(numVal);
        int maxMult = Multiplicity::maxM;
        
        double m = numeric_limits<double>::lowest();
        Assignment assignment(varCards);
        if (numObs > maxMult*covModel.getSosLambda()){
                for (size_t i = 0; i < numVal; i++) {
                        val[i] = (find(assignment.begin(), assignment.end(),maxMult) == assignment.end())? covModel.getLogProb(numObs, assignment): covModel.getMaxMultP();
                        m = max(val[i], m);
                        assignment.increment();
                }
        } else {
                for (size_t i = 0; i < numVal; i++) {
                        val[i] = covModel.getLogProb(numObs, assignment); //*numIndepObs; // TODO is this valuable here??? 
                        m = max(val[i], m);
                        assignment.increment();
                }
        }
        
        for (int i  = 0; i< numVal; i++)
                val[i] = max(exp(val[i] - m), numeric_limits<double>::min());

        return dai::Factor(move(vars), move(val));
}


Factor CRFSolver::createFlowFactor(int sumVarID, int sumCard,
                                   const vector<int>& termVarID,
                                   const vector<int>& termCard,
                                   const vector<bool>& palindromic,
                                   double flowStrength)
{
        assert(termVarID.size() == termCard.size());
        int maxMult = Multiplicity::maxM;

        const double sumOK = flowStrength;
        const double sumNOK = 1.0;

        // create the variable vector
        vector<int> var = { sumVarID };
        var.insert(var.end(), termVarID.begin(), termVarID.end());

        // create the cardinality vector
        vector<int> card = { sumCard };
        card.insert(card.end(), termCard.begin(), termCard.end());

        // create the value vector
        size_t numVal = accumulate(card.begin(), card.end(), 1ull, multiplies<size_t>());
        vector<double> val(numVal);

        Assignment assignment(card);
        for (size_t i = 0; i < val.size(); i++) {
                int sumTerms = 0;
                for (size_t j = 0; j < termCard.size(); j++) {
                        // count palindromic arcs double
                        int c = (palindromic[j]) ? 2 : 1;
                        sumTerms += c * assignment[j+1];
                }

                bool flag = (assignment[0] == sumTerms) || ((sumTerms >= maxMult) && (assignment[0] == maxMult));
                val[i] = (flag) ? log(sumOK) : log(sumNOK);
                assignment.increment();
        }

        return Factor(move(var), move(card), move(val));
}

dai::Factor CRFSolver::createLDFlowFactor(int sumVarID, int sumCard,
                                          const vector<int>& termVarID,
                                          const vector<int>& termCard,
                                          const vector<bool>& palindromic,
                                          double flowStrength)

{
        assert(termVarID.size() == termCard.size());
        int maxMult = Multiplicity::maxM;
        
        const double sumOK = flowStrength;
        const double sumNOK = 1.0;
	
	// sort edges (and their cards, firstmult and palindromic indicator) based on varID (libdai sorts them internally, but the value vector is created based on the order we pass here). NodeID will always be the smallest so that's okay
	vector<size_t> idx(termVarID.size());
	iota(idx.begin(), idx.end(), 0);
	sort(idx.begin(), idx.end(),
	     [&termVarID](size_t i1, size_t i2) {return termVarID[i1] < termVarID[i2];});
	vector<int> termVarIDsort, termCardSort;
	vector<bool> palindromicSort;
	
	for (size_t i: idx){
		termVarIDsort.push_back(termVarID[i]);
		termCardSort.push_back(termCard[i]);
		palindromicSort.push_back(palindromic[i]);
	}

        // create the variable vector
        vector<int> var = { sumVarID };
        var.insert(var.end(), termVarIDsort.begin(), termVarIDsort.end());

        // create the cardinality vector
        vector<int> card = { sumCard };
        card.insert(card.end(), termCardSort.begin(), termCardSort.end());
        
        dai::VarSet vars(dai::Var(sumVarID,sumCard));
        for (int i=0; i< termVarID.size(); i++) {
                vars.insert(dai::Var(termVarIDsort[i], termCardSort[i]));
        }

        // create the value vector
        size_t numVal = accumulate(card.begin(), card.end(), 1ull, multiplies<size_t>());
        vector<double> val(numVal);

        Assignment assignment(card);
        for (size_t i = 0; i < val.size(); i++) {
                int sumTerms = 0;
                for (size_t j = 0; j < termCardSort.size(); j++) {
                        // count palindromic arcs double
                        int c = (palindromic[j]) ? 2 : 1;
                        sumTerms += c * assignment[j+1];
                }
                
                bool flag = (assignment[0] == sumTerms) || ((sumTerms >= maxMult) && (assignment[0] == maxMult));
                val[i] = (flag) ? sumOK : sumNOK;
                assignment.increment();
        }

        return dai::Factor(move(vars), move(val));
}

bool CRFSolver::solveSubgraph(NodeMap<vector<int>>& node2var,
                              EdgeMap<vector<int>>& edge2var,
                              const CovModel& nodeCovModel,
                              const CovModel& edgeCovModel,
                              Multiplicity& resMult, vector<int> targetCRFVar,
                              size_t& factorSize) const
{
         size_t nStrains = nodeCovModel.getNumStrains();
         int maxMult = Multiplicity::maxM;

         // create the singleton factors
         list<Factor> factorList;
         vector<int> card(nStrains * (node2var.size() + edge2var.size()), maxMult + 1);

         for (const auto it : node2var) {        // first the nodes
                 SSNode node = dBG.getSSNode(it.first.getNodeID());
                 const vector<int>& nodeVarIDs = it.second;

                 vector<int> cards(nStrains, maxMult +1);

                 // add a singleton node factor in case:
                 // - the node is bigger than 2k (less correlation with edges) OR
                 // - there are no edges (in case of isolated nodes) OR
                 // - there is a single incoming edge OR
                 // - there is a single outgoing edge
                 if ((node.getMarginalLength() < 2*Kmer::getK()) &&
                         (node.numLeftArcs() > 1) && (node.numRightArcs() > 1) && ! (edge2var.size() == 0))
                         continue;
                 
                 /*// only add a node singleton in case:
                 //      - the node is bigger than 2k (less correlation with edges)
                 //      - there are no edges (in case of isolated nodes)
                 if (node.getMarginalLength() < 2*Kmer::getK() && !edge2var.empty())
                         continue;*/
                 
                 int nIndptObs = max<int>(1, node.getMarginalLength() / (2*Kmer::getK()));
                 Factor F = createSingletonFactor(nodeVarIDs, cards,
                                                  nodeCovModel, node.getAvgCov(), nIndptObs);
                 factorList.push_back(F);

         }

         for (const auto it : edge2var) {        // then the edges
                 const EdgeRep& edge = it.first;
                 const vector<int>& edgeVarIDs = it.second;

                 SSNode src = dBG.getSSNode(edge.getSrcID());
                 Arc* arc = src.rightArc(edge.getDstID());

                 vector<int> cards(nStrains, maxMult+1);
                 
                 // add a singleton edge factor in case:
                 // - an adjacent node is not in the CRF (extremal edge) OR
                 // - both adjacent nodes have multiple edges on that side
                 if ((node2var.find(NodeRep(edge.getSrcID())) != node2var.end()) &&
                         (node2var.find(NodeRep(edge.getDstID())) != node2var.end()) && (! (node2var.size() == 0))) {
                        SSNode src = dBG.getSSNode(edge.getSrcID());
                        if (src.numRightArcs() == 1)
                                continue;
                        SSNode dst = dBG.getSSNode(edge.getDstID());
                        if (dst.numLeftArcs() == 1)
                                continue;
                 }

                 Factor F = createSingletonFactor(edgeVarIDs, cards,
                                                  edgeCovModel, arc->getCov(), 1);
                 factorList.push_back(F);
         }

         // create the flow conservation factors
         if ( !edge2var.empty() ){
         for (const auto it : node2var) {
                 NodeID currID = it.first.getNodeID();
                 const vector<int> nodeVarIDs = it.second;
                 SSNode node = dBG.getSSNode(currID);

                 // left flow conservation factor
                 vector<vector<int>> lEdgeVarID(nStrains); //one for each strain
                 vector<vector<int>> lEdgeVarCard(nStrains); //one for each strain
                 vector<bool> lEdgePalindromic;
                 for (ArcIt lIt = node.leftBegin(); lIt != node.leftEnd(); lIt++) {
                         EdgeRep edge(lIt->getNodeID(), currID);
                         vector<int> edgeVarIDs(edge2var[edge].begin(), edge2var[edge].end());

                         for (size_t strain = 0; strain<nStrains; strain++) {
                             lEdgeVarID[strain].push_back(edgeVarIDs[strain]);
                             lEdgeVarCard[strain].push_back(card[edgeVarIDs[strain]]);
                         }
                         lEdgePalindromic.push_back(lIt->getNodeID() == -currID);
                 }

                 if (!lEdgeVarID[0].empty()) {
                         for (size_t strain = 0; strain < nStrains; strain++){ //one flow factor per strain
                             Factor F = createFlowFactor(nodeVarIDs[strain], card[nodeVarIDs[strain]],
							 lEdgeVarID[strain], lEdgeVarCard[strain],
                                                         lEdgePalindromic, flowStrength);
                             factorList.push_back(F);
                         }

                 }


                 // right flow conservation factor
                 vector<vector<int>> rEdgeVarID(nStrains); //one for each strain
                 vector<vector<int>> rEdgeVarCard(nStrains); //one for each strain
                 vector<bool> rEdgePalindromic;
                 for (ArcIt rIt = node.rightBegin(); rIt != node.rightEnd(); rIt++) {
                         EdgeRep edge(currID, rIt->getNodeID());
                         vector<int> edgeVarIDs(edge2var[edge].begin(), edge2var[edge].end());

                         for (size_t strain=0; strain < nStrains; strain++){
                             rEdgeVarID[strain].push_back(edgeVarIDs[strain]);
                             rEdgeVarCard[strain].push_back(card[edgeVarIDs[strain]]);
                         }
                         rEdgePalindromic.push_back(currID == -rIt->getNodeID());
                 }

                 if (!rEdgeVarID.empty()) {
                         for (size_t strain = 0; strain < nStrains; strain++) { //one flow factor per strain
                             Factor F = createFlowFactor(nodeVarIDs[strain], card[nodeVarIDs[strain]],
							 rEdgeVarID[strain], rEdgeVarCard[strain],
							 rEdgePalindromic, flowStrength);
                             factorList.push_back(F);
                         }
                 }
         }
         }

        // Variable elimination

        //joint assignment over all strain-variables of the target dBG node/arc
        factorSize = maxFactorSize;
        bool retVal = PGMInference::solveVE(factorList, targetCRFVar, factorSize);
        if (!retVal)
                return false;

        vector<int> cardTarget(nStrains);
        for (size_t strain = 0; strain < nStrains; strain++){
            cardTarget[strain] = card[targetCRFVar[strain]];
        }
        resMult = Multiplicity(cardTarget,
                               factorList.front().getVal());
        resMult.normalize();

        //cout << "RESULT: " << resMult << endl;

        return true;
}

bool CRFSolver::checkFlow(NodeRep node, const CovModel& nodeCovModel,
                          const CovModel& edgeCovModel, int graphDepth)
{
        getSubgraph(node, graphDepth);

        for (auto e : nodes) {
                SSNode n = dBG.getSSNode(e);
                vector<int> nodeMult = nodeCovModel.getExpMult(n.getAvgCov());

                vector<int> sumRight(nodeCovModel.getNumStrains(),0);
                for (ArcIt r = n.rightBegin(); r != n.rightEnd(); r++) {
                        // count palindromic arcs double
                        int c = (r->getNodeID() == -e) ? 2 : 1;
                        vector<int> arcMult = edgeCovModel.getExpMult(r->getCov());
                        for(int i = 0; i < nodeCovModel.getNumStrains(); i++)
                                sumRight[i] += c*arcMult[i];
                }

                if (n.numRightArcs() > 0)
                        for (int i = 0; i < nodeCovModel.getNumStrains(); i++){
                                if(nodeMult[i] != sumRight[i])
                                        return false;
                        }

                vector<int> sumLeft(nodeCovModel.getNumStrains(),0);
                for (ArcIt l = n.leftBegin(); l != n.leftEnd(); l++) {
                        // count palindromic arcs double
                        int c = (l->getNodeID() == -e) ? 2 : 1;
                        vector<int> arcMult = edgeCovModel.getExpMult(l->getCov());
                        for(int i = 0; i < nodeCovModel.getNumStrains(); i++)
                                sumLeft[i] += c*arcMult[i];
                }

                if (n.numRightArcs() > 0)
                        for (int i = 0; i < nodeCovModel.getNumStrains(); i++){
                                if(nodeMult[i] != sumLeft[i])
                                        return false;
                        }
        }

        return true;
}

Multiplicity CRFSolver::getNodeMultiplicity(NodeRep node,
                                            const CovModel& nodeCovModel,
                                            const CovModel& edgeCovModel,
                                            int graphDepth)
{
        getSubgraph(node, graphDepth);

        // create a mapping between nodes/edges and CRF variables
        // 1 CRF variable per strain
        int nStrains = nodeCovModel.getNumStrains();
        NodeMap<vector<int>> node2var;
        EdgeMap<vector<int>> edge2var;
        int varCRF = 0;
        for (const auto it : nodes){
            vector<int> vars(nStrains);
            for (size_t strain=0; strain < nStrains; strain++)
                vars[strain] = varCRF++;
            node2var[it] = vars;
        }
        for (const auto it : edges) {
            vector<int> vars(nStrains);
            for (size_t strain = 0; strain < nStrains; strain++)
                vars[strain] = varCRF++;
            edge2var[it] = vars;
        }
        Multiplicity result;
        vector<int> targetCRFVar = node2var[node];

        size_t factorSize;
        if (solveSubgraph(node2var, edge2var, nodeCovModel,
                          edgeCovModel, result, targetCRFVar, factorSize)) {
                perfCounter.totNumCRF++;
                perfCounter.totNumNodes += nodes.size();
                perfCounter.totNumEdges += edges.size();
                perfCounter.totActNBSize += graphDepth;
                perfCounter.totMaxCRFFactor += factorSize;

                return result;
        }

        // fall-back to smaller subgraph if necessary
        return getNodeMultiplicity(node, nodeCovModel, edgeCovModel, graphDepth-1);
}

void CRFSolver::approxSubgraphMult(NodeRep node,
                                   NodeMap<Multiplicity>& nodeMult,
                                   EdgeMap<Multiplicity>& edgeMult,
                                   const CovModel& nodeCovModel,
                                   const CovModel& edgeCovModel,
                                   int graphDepth, bool MAP, size_t numThreads, size_t threadWork)
{
        getSubgraph(node, graphDepth);
        
        string modelIdentifier = to_string(node.getNodeID()) + ".nb" + to_string(graphDepth);
        approxMult(nodeMult, nodeCovModel, edgeMult, edgeCovModel, modelIdentifier, MAP, numThreads, threadWork);
}
        
Multiplicity CRFSolver::getEdgeMultiplicity(EdgeRep edge,
                                            const CovModel& nodeCovModel,
                                            const CovModel& edgeCovModel,
                                            int graphDepth)
{
        getSubgraph(edge, graphDepth);
        
        // create a mapping between nodes/edges and CRF variables
        // 1 CRF variable per strain
        int nStrains = edgeCovModel.getNumStrains();
        NodeMap<vector<int>> node2var;
        EdgeMap<vector<int>> edge2var;
        int varCRF = 0;

        for (const auto it : nodes) {
                vector<int> vars(nStrains);
                for (size_t strain = 0; strain < nStrains; strain++) {
                    vars[strain] = varCRF++;
                }
                node2var[it] = vars;
        }
        for (const auto it : edges){
                vector<int> vars(nStrains);
                for (size_t strain=0; strain < nStrains; strain++){
                    vars[strain] = varCRF++;
                }
                edge2var[it] = vars;
        }

        Multiplicity result;
        vector<int> targetCRFVar = edge2var[edge];
        size_t factorSize = maxFactorSize;
        if (solveSubgraph(node2var, edge2var, nodeCovModel,
                          edgeCovModel, result, targetCRFVar, factorSize)) {

                perfCounter.totNumCRF++;
                perfCounter.totNumNodes += nodes.size();
                perfCounter.totNumEdges += edges.size();
                perfCounter.totActNBSize += graphDepth;
                perfCounter.totMaxCRFFactor += factorSize;

                return result;
        }
        // fall-back to smaller subgraph if necessary
        return getEdgeMultiplicity(edge, nodeCovModel, edgeCovModel, graphDepth-1);
}

dai::FactorGraph CRFSolver::getLibdaiFG(NodeMap<std::vector<int> >& node2var, 
                                        EdgeMap<std::vector<int> >& edge2var,
                                        NodeMap<size_t>& node2fact,
                                        EdgeMap<size_t>& edge2fact,
                                        const CovModel& nodeCovModel, const CovModel& edgeCovModel)
{
        size_t nStrains = nodeCovModel.getNumStrains();
        int maxMult = Multiplicity::maxM;
        
        std::vector<dai::Factor> factorList;
        vector<int> card(nStrains * (node2var.size() + edge2var.size()), maxMult + 1);
        int nodeSingletonF = 0;
        int arcSingletonF = 0;
        int factorID = 0;
        int flowF = 0;
        
        // first the nodes
        for (const auto it : node2var) {        
                SSNode node = dBG.getSSNode(it.first.getNodeID());
                const vector<int>& nodeVarIDs = it.second;
                
                vector<int> cards(nStrains, maxMult +1);
                
                size_t numNbArcs = node.numLeftArcs() + node.numRightArcs();
                dai::Factor F;
                // add a singleton node factor in case:
                // - the node is bigger than 2k (less correlation with edges) OR
                // - there are no edges (in case of isolated nodes) OR
                // - there is a single incoming edge OR
                // - there is a single outgoing edge
                if ((node.getMarginalLength() < 2*Kmer::getK()) &&
                        (node.numLeftArcs() > 1) && (node.numRightArcs() > 1 && ! (edge2var.size() == 0))){
                        
                        F = createLDUniformSingletonFactor(nodeVarIDs, cards);
                } else {
                        int nIndptObs = max<int>(1, node.getMarginalLength() / (2*Kmer::getK()));
                        
                        F = createLDSingletonFactor(nodeVarIDs, cards,
                                                                nodeCovModel, node.getAvgCov(),
                                                                nIndptObs);
                }        
                nodeSingletonF++;
                node2fact[it.first] = factorID;
                factorID++;
                factorList.push_back(F);
        }
        
        for (const auto it : edge2var) {        // then the edges
                const EdgeRep& edge = it.first;
                const vector<int>& edgeVarIDs = it.second;
                
                SSNode src = dBG.getSSNode(edge.getSrcID());
                Arc* arc = src.rightArc(edge.getDstID());
                
                vector<int> cards(nStrains, maxMult+1);
                
                dai::Factor F;
                SSNode dst = dBG.getSSNode(edge.getDstID());
                // add a singleton edge factor in case:
                // - an adjacent node is not in the CRF (extremal edge) OR
                // - both adjacent nodes have multiple edges on that side
                if ((node2var.find(NodeRep(edge.getSrcID())) != node2var.end()) &&
                        (node2var.find(NodeRep(edge.getDstID())) != node2var.end()) && (! (node2var.size() == 0)) && (src.numRightArcs() == 1 || dst.numLeftArcs() == 1)) {
                        F = createLDUniformSingletonFactor(edgeVarIDs, cards);
                } else {
                        
                        F = createLDSingletonFactor(edgeVarIDs, cards,
                                                        edgeCovModel, arc->getCov(), 1);
                }
                factorList.push_back(F);
                edge2fact[it.first] = factorID;
                factorID++;
                arcSingletonF++;
        }
        
        // create the flow conservation factors
        if ( !edge2var.empty() ){
                for (const auto it : node2var) {
                        NodeID currID = it.first.getNodeID();
                        const vector<int> nodeVarIDs = it.second;
                        SSNode node = dBG.getSSNode(currID);
                        
                        // left flow conservation factor
                        vector<vector<int>> lEdgeVarID(nStrains); //one for each strain
                        vector<vector<int>> lEdgeVarCard(nStrains); //one for each strain
                        vector<bool> lEdgePalindromic;
                        for (ArcIt lIt = node.leftBegin(); lIt != node.leftEnd(); lIt++) {
                                EdgeRep edge(lIt->getNodeID(), currID);
                                vector<int> edgeVarIDs(edge2var[edge].begin(), edge2var[edge].end());
                                
                                for (size_t strain = 0; strain<nStrains; strain++) {
                                        lEdgeVarID[strain].push_back(edgeVarIDs[strain]);
                                        lEdgeVarCard[strain].push_back(card[edgeVarIDs[strain]]);
                                }
                                lEdgePalindromic.push_back(lIt->getNodeID() == -currID);
                        }
                        
                        if (!lEdgeVarID[0].empty()) {
                                for (size_t strain = 0; strain < nStrains; strain++){ //one flow factor per strain
                                        dai::Factor F = createLDFlowFactor(nodeVarIDs[strain], card[nodeVarIDs[strain]], 
                                                                           lEdgeVarID[strain], lEdgeVarCard[strain], lEdgePalindromic, flowStrength);
                                        factorList.push_back(F);
                                        flowF++;
                                }
                        }
                        
                        
                        // right flow conservation factor
                        vector<vector<int>> rEdgeVarID(nStrains); //one for each strain
                        vector<vector<int>> rEdgeVarCard(nStrains); //one for each strain
                        vector<bool> rEdgePalindromic;
                        for (ArcIt rIt = node.rightBegin(); rIt != node.rightEnd(); rIt++) {
                                EdgeRep edge(currID, rIt->getNodeID());
                                vector<int> edgeVarIDs(edge2var[edge].begin(), edge2var[edge].end());
                                
                                for (size_t strain=0; strain < nStrains; strain++){
                                        rEdgeVarID[strain].push_back(edgeVarIDs[strain]);
                                        rEdgeVarCard[strain].push_back(card[edgeVarIDs[strain]]);
                                }
                                rEdgePalindromic.push_back(currID == -rIt->getNodeID());
                        }
                        
                        if (!rEdgeVarID.empty()) {
                                for (size_t strain = 0; strain < nStrains; strain++) { //one flow factor per strain
                                        dai::Factor F = createLDFlowFactor(nodeVarIDs[strain], card[nodeVarIDs[strain]],
                                                                           rEdgeVarID[strain], rEdgeVarCard[strain], rEdgePalindromic, flowStrength);
                                        factorList.push_back(F);
                                        flowF++;
                                }
                        }
                }
        }
        
        cout << "libDai Factor graph contains " << nodeSingletonF << " node singleton factors, "
        << arcSingletonF << " arc singleton factors, and " 
        << flowF << " flow factors" << endl;
        
        return dai::FactorGraph(factorList);
}

void CRFSolver::approxMult(NodeMap<Multiplicity>& nodeMult,
                           const CovModel& nodeCovModel,
                           EdgeMap<Multiplicity>& edgeMult,
                           const CovModel& edgeCovModel,
                           string modelID,
                           bool MAP, size_t numThreads, size_t threadWork)
{
        //CRFSolver solver(dBG, maxFactorSize, flowStrength);
        
        // Compute a CRF of all nodes and arcs in the de Bruijn graph
        int nStrains = nodeCovModel.getNumStrains();
        NodeMap<vector<int>> node2var;
        NodeMap<size_t> node2fact;
        EdgeMap<vector<int>> edge2var;
        EdgeMap<size_t> edge2fact;
        //vector<NodeRep> allNodes = dBG.getNodeReps(dBG.getNumNodes());
        //vector<EdgeRep> allEdges = dBG.getEdgeReps(dBG.getNumArcs());
        
        int varCRF = 0;
        for (const auto it : nodes){
                vector<int> vars(nStrains);
                //cout << it << ":";
                for (size_t strain=0; strain < nStrains; strain++){
                        vars[strain] = varCRF++;
                //        cout << "\t" << varCRF-1;
                }
                //cout << "\n";
                node2var[it] = vars;
                node2fact[it] = -1;
        }
        for (const auto it : edges) {
                vector<int> vars(nStrains);
                for (size_t strain = 0; strain < nStrains; strain++)
                        vars[strain] = varCRF++;
                edge2var[it] = vars;
                edge2fact[it] = -1;
        }
        
        //cout << "Getting libDAI FG representation ...." << endl;
        //Util::startChrono();
        dai::FactorGraph fullFG = getLibdaiFG(node2var, edge2var,
                                              node2fact, edge2fact,
                                              nodeCovModel, edgeCovModel);
        //cout << "... done (" << Util::stopChronoStr() << ")" << endl;
        
        
        // Properties for belief propagation TODO: lift this to settings or main.cpp and check appropriate values
        dai::PropertySet opts;
        if (Util::fileExists("libdai.props")) {
                ifstream options;
                options.open("libdai.props");
                string line;
                while(getline(options,line)) {
                        if (line.compare(0,1,"%") == 0)
                                continue;
                        size_t tab_pos = line.find("\t");
                        opts.set(line.substr(0,tab_pos),line.substr(tab_pos + 1, line.length()));
                }
        } else {
                opts.set("maxiter", (size_t) 500);
                opts.set("tol", 1e-3);
                opts.set("logdomain", true);
                opts.set("verbose", (size_t) 1);
                opts.set("maxtime", string("1800"));
                opts.set("updates", string("SEQMAX0L"));
                //opts.set("splashsize", (size_t) 5);
                opts.set("weightdecay", true);
                opts.set("forcestop", true);
                //opts.set("nthreads", (size_t) 1);
                //opts.set("resinit", string("UNIFORM"));
                //opts.set("damping", 0.2);
        }
        
        opts.set("modelcount", modelID);
        
        //TODO fix factor graph visualisation
        if (opts.getStringAs<size_t>("verbose") == 4){
                //map<size_t, size_t> nodevars;
                //for (const auto it : nodes)
                //        nodevars[node2var[it]] = it.getNodeID();
                cout << "writing FG representation ...." << endl;
                Util::startChrono();
                string fgname = "factorgraph." + modelID + ".fg";
                fullFG.WriteToFile(fgname.c_str(), 4);
                //ofstream cytnodes, cytarcs;
                //cytnodes.open("factorgraph.cyt.nodes");
                //cytarcs.open("factorgraph.cyt.arcs");
                //fullFG.dBGfg2Cytoscape(cytarcs, cytnodes,nodevars);
                //cytnodes.close();
                //cytarcs.close();
                cout << "... done (" << Util::stopChronoStr() << ")" << endl;
        }
        map<int, dai::Var> label2var;
        double finalMaxDiff;
        vector<size_t> maxProbAssignment;
        map<int, int> label2idx;
        //cout << "Starting libDAI LBP ...." << endl;
        //Util::startChrono();
        //if( opts.getStringAs<size_t>("nthreads") > 1 ){
        if ( ! opts.getStringAs<std::string>("updates").compare("SPLASH")){
                dai::ParallelBP bp;
                bp = dai::ParallelBP(fullFG, opts); 
                bp.init();
                finalMaxDiff = bp.run();
                
                // Get all beliefs by variable ID
                for (size_t i=0; i < fullFG.nrVars(); ++i) {
                        dai::Var v = fullFG.var(i);
                        int label = v.label();
                        label2var[label] = v;
                        label2idx[label] = i;
                }
                
                if(MAP) {
                        // Extract multiplicities for the nodes
                        for (int id = 0; id < nodes.size(); id ++) {
                                NodeRep n = nodes[id];
                                vector<int> label = node2var[n];
                                vector<int> assignment;
                                assignment.reserve(label.size());
                                for (int l: label)
                                        assignment.push_back(maxProbAssignment[label2idx[l]]);
                                Multiplicity mult(assignment);
                                nodeMult[n] = mult;
                        }
                        
                        // Extract multiplicities for the edges
                        for (int id = 0; id < edges.size(); id ++) {
                                EdgeRep n = edges[id];
                                vector<int> label = edge2var[n];
                                vector<int> assignment;
                                assignment.reserve(label.size());
                                for (int l: label)
                                        assignment.push_back(maxProbAssignment[label2idx[l]]);
                                Multiplicity mult(assignment);
                                edgeMult[n] = mult;
                        }
                        
                } else {
                        for (auto &it : nodes){
                                NodeRep n = it;
                                nodeMult[n] = getMultFromFG(node2var.find(n)->second,
                                                            label2var,
                                                            node2fact.find(n)->second,
                                                            bp,
                                                            nodeCovModel);
                        }
                        
                        for (auto &it : edges){
                                EdgeRep e = it;
                                edgeMult[e] = getMultFromFG(edge2var.find(e)->second,
                                                            label2var,
                                                            edge2fact.find(e)->second,
                                                            bp,
                                                            edgeCovModel);
                        }
                }
                //cout << "... done (" << Util::stopChronoStr() << ")" << endl;
        } else {
                dai::BP bp;
                if(MAP){
                        opts.set("updates", string("SEQMAX"));
                        bp = dai::BP(fullFG, opts("inference",string("MAXPROD")));
                } else {
                        bp = dai::BP(fullFG, opts); 
                }
                bp.init();
                finalMaxDiff = bp.run();
                /*if (finalMaxDiff > opts.getAs<double>("tol"))
                 *               for (const auto it : nodes){
                 *                       cout << varCRF -1 << "\t" << it.getNodeID() << endl;
        }
        */
                if(MAP)
                        maxProbAssignment = bp.findMaximum();
                
                //cout << "... LBP finished in " << Util::stopChronoStr() << endl;
                cout << "LogZ(x): " << bp.logZ() << endl;
                logZ = bp.logZ();
                // Get all beliefs by variable ID
                //cout << "Retrieving all beliefs ...." << endl;
                //Util::startChrono();
                // Get all Var by variable ID
                for (size_t i=0; i < fullFG.nrVars(); ++i) {
                        dai::Var v = fullFG.var(i);
                        int label = v.label();
                        label2var[label] = v;
                        label2idx[label] = i;
                }
                
                
                if(MAP) {
                        // Extract multiplicities for the nodes
                        for (int id = 0; id < nodes.size(); id ++) {
                                NodeRep n = nodes[id];
                                vector<int> label = node2var[n];
                                vector<int> assignment;
                                assignment.reserve(label.size());
                                for (int l: label)
                                        assignment.push_back(maxProbAssignment[label2idx[l]]);
                                Multiplicity mult(assignment);
                                nodeMult[n] = mult;
                        }
                        
                        // Extract multiplicities for the edges
                        for (int id = 0; id < edges.size(); id ++) {
                                EdgeRep n = edges[id];
                                vector<int> label = edge2var[n];
                                vector<int> assignment;
                                assignment.reserve(label.size());
                                for (int l: label)
                                        assignment.push_back(maxProbAssignment[label2idx[l]]);
                                Multiplicity mult(assignment);
                                edgeMult[n] = mult;
                        }
                        
                } else {
                        for (auto &it : nodes){
                                NodeRep n = it;
                                if (nodeMult.find(n) != nodeMult.end())
                                        nodeMult[n] = getMultFromFG(node2var.find(n)->second,
                                                                        label2var,
                                                                        node2fact.find(n)->second,
                                                                        bp,
                                                                        nodeCovModel);
                        }
                        
                        for (auto &it : edges){
                                EdgeRep e = it;
                                if (edgeMult.find(e) != edgeMult.end())
                                        edgeMult[e] = getMultFromFG(edge2var.find(e)->second,
                                                                label2var,
                                                                edge2fact.find(e)->second,
                                                                bp,
                                                                edgeCovModel);
                        }
                }
                //cout << "... done (" << Util::stopChronoStr() << ")" << endl;
        }
                
}

// ============================================================================
// ESTIMATE NODE/EDGE MULTIPLICITIES USING CONDITIONAL RANDOM FIELDS
// ============================================================================

void CRFMult::checkNodeFlow(CRFSolver& solver, WorkLoadBalancer& wlb,
                            const vector<NodeRep>& nodes,
                            const CovModel& nodeCovModel,
                            const CovModel& edgeCovModel,
                            vector<bool>& flowOK) const
{
        size_t idxBegin, idxEnd;
        while (wlb.getChunk(idxBegin, idxEnd))
                for (size_t id = idxBegin; id < idxEnd; id++)
                        flowOK[id] = solver.checkFlow(nodes[id],
                                                      nodeCovModel,
                                                      edgeCovModel,
                                                      maxGraphDepth);
}

void CRFMult::computeNodeMult(CRFSolver& solver,
                              WorkLoadBalancer& wlb,
                              const vector<NodeRep>& nodes,
                              const CovModel& nodeCovModel,
                              const CovModel& edgeCovModel,
                              NodeMap<Multiplicity>& nodeMult) const
{
        size_t idxBegin, idxEnd;
        while (wlb.getChunk(idxBegin, idxEnd))
                for (size_t id = idxBegin; id < idxEnd; id++){
                        NodeRep n = nodes[id];
                        nodeMult[n] = solver.getNodeMultiplicity(n,
                                                                  nodeCovModel,
                                                                  edgeCovModel,
                                                                  maxGraphDepth);
                }
}

void CRFMult::computeEdgeMult(CRFSolver& solver,
                              WorkLoadBalancer& wlb,
                              const vector<EdgeRep>& edges,
                              const CovModel& nodeCovModel,
                              const CovModel& edgeCovModel,
                              EdgeMap<Multiplicity>& edgeMult) const
{
        size_t idxBegin, idxEnd;
        while (wlb.getChunk(idxBegin, idxEnd))
                for (size_t id = idxBegin; id < idxEnd; id++){
                        const EdgeRep& e = edges[id];
                        edgeMult[e] = solver.getEdgeMultiplicity(e,
                                                                  nodeCovModel,
                                                                  edgeCovModel,
                                                                  maxGraphDepth);
                }
}

/*void CRFSolver::retrieveNodeMultFromFG(WorkLoadBalancer& wlb, 
                            const std::map<int, dai::Var>& label2var,
                            const NodeMap<std::vector<int>> node2var,
                            const std::vector<NodeRep> nodes,
                            const dai::BP& bp,
                            NodeMap<Multiplicity>& nodeMult,
                            const CovModel& nodeCovModel) const
{
        size_t idxBegin, idxEnd;
        while (wlb.getChunk(idxBegin, idxEnd))
                for (size_t id = idxBegin; id < idxEnd; id++){
                        NodeRep n = nodes[id];
                        nodeMult[n] = getMultFromFG(node2var.find(n)->second, 
                                                         label2var,
                                                         bp,
                                                         nodeCovModel);
                }
}

void CRFSolver::retrieveEdgeMultFromFG(WorkLoadBalancer& wlb, 
                                     const std::map<int, dai::Var>& label2var,
                                     const EdgeMap<std::vector<int>> edge2var,
                                     const std::vector<EdgeRep> edges,
                                     const dai::BP& bp,
                                     EdgeMap<Multiplicity>& edgeMult,
                                     const CovModel& edgeCovModel) const
{
        size_t idxBegin, idxEnd;
        while (wlb.getChunk(idxBegin, idxEnd))
                for (size_t id = idxBegin; id < idxEnd; id++){
                        EdgeRep n = edges[id];
                        edgeMult[n] = getMultFromFG(edge2var.find(n)->second,
                                                     label2var,
                                                     bp,
                                                     edgeCovModel);
                }
}

Multiplicity CRFSolver::getMultFromFG(const std::vector<int>& labels,
                                        const std::map<int, dai::Var>& label2var,
                                        const dai::BP& bp,
                                        const CovModel& cm) const
{
        dai::VarSet vs;
        for (int label: labels)
                vs.insert(label2var.find(label)->second);
        dai::Factor f = bp.belief(vs);
        vector<double> logprob = f.log(true).p().p();
        for (int i = 0; i < logprob.size(); i++)
                if(!logprob[i])
                        logprob[i] = log(DOUBLE_SMALL);
        Multiplicity mult(vector<int>(cm.getNumStrains(), cm.getMaxMult()+1), logprob);
        return mult;
}*/

Multiplicity CRFSolver::getMultFromFG(const std::vector<int>& labels,
                                      const std::map<int, dai::Var>& label2var,
                                      const size_t factorID,
                                      const dai::BP& bp,
                                      const CovModel& cm) const
{
        dai::VarSet vs;
        for (int label: labels)
                vs.insert(label2var.find(label)->second);
        dai::Factor f = bp.beliefF(factorID).marginal(vs);
        vector<double> logprob = f.log(true).p().p();
        for (int i = 0; i < logprob.size(); i++)
                if(!logprob[i])
                        logprob[i] = log(DOUBLE_SMALL);
                Multiplicity mult(vector<int>(cm.getNumStrains(), cm.getMaxMult()+1), logprob);
        return mult;
}

Multiplicity CRFSolver::getMultFromFG(const std::vector<int>& labels,
                                      const std::map<int, dai::Var>& label2var,
                                      const size_t factorID,
                                      const dai::ParallelBP& bp,
                                      const CovModel& cm) const
{
        dai::VarSet vs;
        for (int label: labels)
                vs.insert(label2var.find(label)->second);
        dai::Factor f = bp.beliefF(factorID).marginal(vs);
        vector<double> logprob = f.log(true).p().p();
        for (int i = 0; i < logprob.size(); i++)
                if(!logprob[i])
                        logprob[i] = log(DOUBLE_SMALL);
                Multiplicity mult(vector<int>(cm.getNumStrains(), cm.getMaxMult()+1), logprob);
        return mult;
}

void CRFMult::checkFlow(const vector<NodeRep>& nodes,
                        vector<bool>& flowOK,
                        const CovModel& nodeCovModel,
                        const CovModel& edgeCovModel) const
{
        cout << "\tChecking flow for " << nodes.size() << " nodes (subgraph "
                "depth: " << maxGraphDepth << ")" << endl;
        vector<thread> wt(numThreads);
        vector<CRFSolver> solver(numThreads, CRFSolver(dBG,
                                                       maxFactorSize,
                                                       flowStrength));

        // assign a multiplicity to the nodes
        WorkLoadBalancer nodeWlb(0, nodes.size(), threadWork, "\tProcessing nodes");
        for (size_t i = 0; i < wt.size(); i++)
                wt[i] = thread(&CRFMult::checkNodeFlow, this, ref(solver[i]),
                               ref(nodeWlb), cref(nodes), cref(nodeCovModel),
                               cref(edgeCovModel), ref(flowOK));

        // wait for worker threads to finish
        for_each(wt.begin(), wt.end(), mem_fn(&thread::join));
}

void CRFMult::computeMult(NodeMap<Multiplicity>& nodeMult,
                          const CovModel& nodeCovModel,
                          EdgeMap<Multiplicity>& edgeMult,
                          const CovModel& edgeCovModel) const
{
        totPerfCounter.reset();
        
        cout << "Computing multiplicity for " << nodeMult.size()
        << " nodes and " << edgeMult.size() << " edges (subgraph "
        "depth: " << maxGraphDepth << ")" << endl;

        vector<thread> wt(numThreads);
        vector<CRFSolver> solver(numThreads, CRFSolver(dBG, maxFactorSize, flowStrength));

        // assign a multiplicity to the nodes
        vector<NodeRep> nodes; nodes.reserve(nodeMult.size());
        for (const auto& it : nodeMult)
                nodes.push_back(it.first);
        WorkLoadBalancer nodeWlb(0, nodes.size(), threadWork, "\tProcessing nodes");
        for (size_t i = 0; i < wt.size(); i++)
                wt[i] = thread(&CRFMult::computeNodeMult, this, ref(solver[i]),
                               ref(nodeWlb), cref(nodes), cref(nodeCovModel),
                               cref(edgeCovModel), ref(nodeMult));

        // wait for worker threads to finish
        for_each(wt.begin(), wt.end(), mem_fn(&thread::join));
        nodes.clear();

        // assign a multiplicity to the edges
        vector<EdgeRep> edges; edges.reserve(edgeMult.size());
        for (const auto& it : edgeMult)
                edges.push_back(it.first);
        WorkLoadBalancer edgeWlb(0, edges.size(), threadWork, "\tProcessing edges");
        for (size_t i = 0; i < wt.size(); i++)
                wt[i] = thread(&CRFMult::computeEdgeMult, this, ref(solver[i]),
                               ref(edgeWlb), cref(edges), cref(nodeCovModel),
                               cref(edgeCovModel), ref(edgeMult));

        // wait for worker threads to finish
        for_each(wt.begin(), wt.end(), mem_fn(&thread::join));

        // merge the performance counters
        for (const auto& it : solver)
                totPerfCounter.merge(it.perfCounter);
        
        /*cout << "Average number of nodes in neighbourhood CRF: " << totPerfCounter.getAvgNumNodes() << endl;
        cout << "Average number of edges in neighbourhood CRF: " << totPerfCounter.getAvgNumEdges() << endl;
        cout << "Average neighbourhood size actually used: " << totPerfCounter.getAvgActNBSize() << endl;
        cout << "Average maximum factor size: "  << totPerfCounter.getAvgMaxCRFFactor() << endl;*/
}

void CRFMult::approxMultAll(NodeMap<Multiplicity>& nodeMult,
                            const CovModel& nodeCovModel,
                            EdgeMap<Multiplicity>& edgeMult,
                            const CovModel& edgeCovModel,
                            double& logZ,
                            bool MAP, bool singleCRF) const
{
        size_t ctr = 0;
        logZ = -1.0;
        if (singleCRF) {
                vector<NodeRep> nodes = dBG.getNodeReps(dBG.getNumValidNodes());
                vector<EdgeRep> edges = dBG.getEdgeReps(dBG.getNumValidArcs());
                CRFSolver solver(dBG, maxFactorSize, flowStrength, nodes, edges);
                solver.approxMult(nodeMult, nodeCovModel,
                                  edgeMult, edgeCovModel,
                                  to_string(ctr), MAP, numThreads, threadWork);
                logZ = solver.getCurrentLogZ();
        } else {
                Bitvec handled(dBG.getNumNodes()+1);
                for(const auto& it : nodeMult){
                        NodeRep nr = it.first;
                        if (handled[abs(nr.getNodeID())])
                                continue;
                        CRFSolver solver(dBG, maxFactorSize, flowStrength);
                        vector<NodeID> inSubgraph = solver.getFullConnectedComponent(nr);
                        if(inSubgraph.size() == 1){
                                for (EdgeRep er: solver.getEdgeReps()){
                                        // Selfloop / palindromic
                                        edgeMult[er] = solver.getEdgeMultiplicity(er, nodeCovModel, edgeCovModel, 0);
                                }
                                nodeMult[nr] = solver.getNodeMultiplicity(nr, nodeCovModel, edgeCovModel, 0);
                                continue;
                        }
                        solver.approxMult(nodeMult, nodeCovModel,
                                        edgeMult, edgeCovModel,
                                        to_string(ctr), MAP, numThreads, threadWork);
                        double lz = solver.getCurrentLogZ();
                        //if(lz > logZ)
                                logZ += lz;
                        for(NodeID id : inSubgraph)
                                handled[abs(id)] = true;
                        ctr++;
                }
                // Check if there are edges in subgraphs that were not handles yet 
                for(const auto& it : edgeMult){
                        EdgeRep er = it.first;
                        if (handled[abs(er.getSrcID())])
                                continue;
                        NodeRep nr(er.getSrcID());
                        CRFSolver solver(dBG, maxFactorSize, flowStrength);
                        vector<NodeID> inSubgraph = solver.getFullConnectedComponent(nr);
                        if(inSubgraph.size() == 1){
                                for (EdgeRep er: solver.getEdgeReps()){
                                        // Selfloop / palindromic
                                        edgeMult[er] = solver.getEdgeMultiplicity(er, nodeCovModel, edgeCovModel, 0);
                                }
                                nodeMult[nr] = solver.getNodeMultiplicity(nr, nodeCovModel, edgeCovModel, 0);
                                continue;
                        }
                        solver.approxMult(nodeMult, nodeCovModel,
                                          edgeMult, edgeCovModel,
                                          to_string(ctr), MAP, numThreads, threadWork);
                        double lz = solver.getCurrentLogZ();
                        //if(lz > logZ)
                                logZ += lz;
                        for(NodeID id : inSubgraph)
                                handled[abs(id)] = true;
                        ctr++;
                }
        }
        
}

bool CRFMult::Estep(NodeMap<Multiplicity>& nodeMult,
                    const CovModel& nodeCovModel,
                    EdgeMap<Multiplicity>& edgeMult,
                    const CovModel& edgeCovModel,
                    double epsilon, double converged, double& prevLogZ,
                    bool approxInf, bool map) const
{

        bool changes = false;
        size_t numNodeChange = 0;
        size_t numEdgeChange = 0;
        // store the old expected multiplicities to check for convergence
        NodeMap<vector<int>> oldNodeMult(nodeMult.size());
        for (const auto& it : nodeMult)
                oldNodeMult[it.first] = it.second.getExpMult();
        EdgeMap<vector<int>> oldEdgeMult(edgeMult.size());
        for (const auto& it : edgeMult)
                oldEdgeMult[it.first] = it.second.getExpMult();
        

        double logZ = -1.0;
        // compute (normalized) multiplicities given the model
        if (approxInf)
                approxMultAll(nodeMult, nodeCovModel,
                              edgeMult, edgeCovModel,logZ, map);
        else
                computeMult(nodeMult, nodeCovModel,
                    edgeMult, edgeCovModel);

        for (auto it : nodeMult)
                if (it.second.getExpMult() != oldNodeMult[it.first])
                        numNodeChange++;

        for (auto it : edgeMult)
                if (it.second.getExpMult() != oldEdgeMult[it.first])
                        numEdgeChange++;
                

        // convergence
        cout << "number of node changes: " << numNodeChange << endl;
        cout << "number of edge changes: " << numEdgeChange << endl;
        
        //double converged = epsilon;
        if (approxInf)
                converged *= 2.5;
        
        double lr = 1.0 - prevLogZ / logZ;
        
        cout << "current 1-LR: " << std::setprecision(6) << lr << endl;
        prevLogZ = max(0.0, logZ);

        // convergence
        return ((((double)numNodeChange/nodeMult.size() + (double)numEdgeChange/edgeMult.size()) / 3.0 > converged) && (approxInf ? (lr > epsilon) : true));
}

bool CRFMult::Mstep(const NodeMap<Multiplicity>& nodeMult,
                    CovModel& nodeCovModel,
                    const EdgeMap<Multiplicity>& edgeMult,
                    CovModel& edgeCovModel, double epsilon, bool eqW, bool fixedZero) const
{
        bool retVal = true; // FIXME ?

        int PC = 1; //TODO: proportional to number of nodes/edges?
        int maxMult = Multiplicity::maxM;
        int nstrains = nodeCovModel.getNumStrains();

        // ====================================================================
        // Compute the node model parameters
        // ====================================================================
        vector<int> zeroMult = vector<int>(nodeCovModel.getNumStrains(),0);
        double nodeErrorLambda = nodeCovModel.getErrLambda();
        double nodeErrorODF = nodeCovModel.getErrorODF();
        double nodeErrorWeight = nodeCovModel.getWeight(zeroMult);
        
        if (! fixedZero){
                // We fit a negative binomial (NB) distribution to the error histogram.
                // The error histogram is always truncated: a) coverage zero is never
                // observed and b) low-coverage k-mers might have been removed by
                // a preprocessing tool like BCALM. We therefore use EM to fit the NB
                // parameters and infer the missing values from the spectrum.
                map<unsigned int, double> errorHist;

                for (const auto& it : nodeMult) {
                        SSNode node = dBG.getSSNode(it.first);
                        double p_zero = it.second[zeroMult];
                        if (p_zero > DOUBLE_SMALL) {
                                double f = node.getAvgCov() - floor(node.getAvgCov());
                                errorHist[node.getAvgCov()] += (1.0-f) * p_zero; // * node.getMarginalLength();
                                errorHist[node.getAvgCov() + 1] += f * p_zero; // * node.getMarginalLength();
                        }
                }


                // We truncate the node histogram to the provided -abundance-min value.
                // This is relevant when using q-mer counts.
                unsigned int smallSupport = (dBG.getAbundanceMin() < 0)?
                        errorHist.begin()->first : dBG.getAbundanceMin();
                for (unsigned int k = 0; k < smallSupport; k++)
                        errorHist.erase(k);
                
                /*cout << "Current NODE error hist: \n";
                 f or (unsigned int k = 0; k < 30; k++)               *
                 cout << k << "\t" << errorHist[k] << "\n";*/

                // The below EM procedure might not convergence, but the obtained NB
                // parameters should be sensible in all cases: only small support
                // values are inferred, therefore, the mean should never be larger than
                // the mean of the error histogram (and thus be close to zero).
                int maxIter = 10.0 / epsilon;
                int nIter = Util::fitTruncNegBinomEM(errorHist, nodeErrorLambda,
                                                nodeErrorODF, nodeErrorWeight,
                                                epsilon, maxIter);

                if (nIter > maxIter)
                        cout << "\tWARNING: EM algorithm to fit node error model "
                                "did not converge\n";
                else
                        cout << "\tEM algorithm to fit node error model converged "
                        << "after " << nIter << " iteration(s)" << endl;
        }
        // Compute node average and weights
        vector<double> wN(nodeCovModel.getK(), PC);
        vector<double> wN_k(nodeCovModel.getK(), PC);
        wN[0] = nodeErrorWeight;
        wN_k[0] = nodeErrorWeight;


        vector<double> nodeLambdas(nodeCovModel.getNumStrains(), 0.0);

        vector<double> mtm(nodeCovModel.getNumStrains()*nodeCovModel.getNumStrains(), 0.0);
        vector<double> mtX(nodeCovModel.getNumStrains(), 0.0);

        // for each (1,0,0), (0,1,0), (0,0,1) add PC,  make sure mtm is positive definite
        for (int strain = 0; strain < nodeCovModel.getNumStrains(); strain++){
            mtm[strain * (nodeCovModel.getNumStrains()+1)] += 1.0;
            mtX[strain] += nodeCovModel.getFraction(strain) * nodeCovModel.getSosLambda(); //prev lambda
            wN[pow(nodeCovModel.getC(), strain)] += 1.0;
            wN_k[pow(nodeCovModel.getC(), strain)] += 1.0;
        }

        for (const auto& it : nodeMult) {
                SSNode node = dBG.getSSNode(it.first);
                vector<int> expMult = it.second.getExpMult();
                double cov = node.getAvgCov(); //get average kmer coverage

                if (accumulate(expMult.begin(), expMult.end(), 0) == 0 || cov > (maxMult)*nodeCovModel.getSosLambda()) //error node
                    continue;

                //counts for weight matrix
                // soft assignments: TODO: Compute intensive procedure?
                Assignment mult(vector<int>(nodeCovModel.getNumStrains(), nodeCovModel.getC()));
                for (int idx=1; idx < nodeCovModel.getK(); idx ++){ //don't compute wNode[0], already computed
                        mult.increment();
                        wN[idx] += it.second[mult] ;//* node.getMarginalLength(); //soft assignment
                        wN_k[idx] += it.second[mult] * node.getMarginalLength(); //soft assignment
                }

                if (*max_element(expMult.begin(), expMult.end()) > 1)
                        continue;
                // fill OLS vector and matrix
                for (size_t strain = 0; strain < nodeCovModel.getNumStrains(); strain++) {
                	mtX[strain] += expMult[strain] * node.getCov();// * node.getMarginalLength(); //totalCov of all kmers
                	for (size_t j = strain; j < nodeCovModel.getNumStrains(); j++)
                        	mtm[strain * nodeCovModel.getNumStrains() + j] += expMult[strain] * expMult[j] * node.getMarginalLength();// * node.getMarginalLength();
                }

        }

        for (size_t j = 0; j < nodeCovModel.getNumStrains()-1; j++){
            for (size_t i = j+1; i < nodeCovModel.getNumStrains(); i++) {
                mtm[i*nodeCovModel.getNumStrains() + j] = mtm[j*nodeCovModel.getNumStrains() + i];
            }
        }

        Util::choleskyLinSolve(mtm,mtX);
        /*cout << "NodeLambdas based on linear system of all observations: " << "[ ";
        for (size_t strain=0; strain < nodeCovModel.getNumStrains(); strain++){
            cout << mtX[strain] << " ";
        }
        cout  <<  "]" << endl;*/

        for (int strain = 0; strain < nodeCovModel.getNumStrains(); strain++){
            if (mtX[strain] < 0){
		retVal=false;
                //assign previous smallest lambda
             	vector<double> frac = nodeCovModel.getFractions();
                mtX[strain] = *min_element(frac.begin(), frac.end()) * nodeCovModel.getSosLambda();
                cerr << "Warning: encountered negative lambda in node LambdaEstimates \n"
                        "set to smallest lambda in previous EM iteration \n"
			"will stop EM after this iteration, possibly less strains are present.\n";
            }
        }

        //Use lambda estimate based on OLS of all observations
        nodeLambdas = mtX;
	sort(nodeLambdas.begin(), nodeLambdas.end());
        double sumNodeLambdas = accumulate(nodeLambdas.begin(), nodeLambdas.end(), 0.0);
	// Update lambda_err to avoid it getting too big.
        nodeErrorLambda = min(sumNodeLambdas*0.025, nodeErrorLambda);

        // ====================================================================
        // Compute the edge model parameters
        // ====================================================================

        double edgeErrorLambda = edgeCovModel.getErrLambda();
        double edgeErrorODF = edgeCovModel.getErrorODF();
        double edgeErrorWeight = edgeCovModel.getWeight(zeroMult);
        if (! fixedZero){
                map<unsigned int, double> errorHist; 
                for (const auto& it : edgeMult) {
                        Arc& arc = dBG.getArc(dBG.getArcID(it.first)) ;
                        double p_zero = it.second[zeroMult];
                        if (p_zero > DOUBLE_SMALL) {
                                double f = arc.getCov() - floor(arc.getCov());
                                errorHist[arc.getCov()] += (1.0-f) * p_zero;
                                errorHist[arc.getCov() + 1] += f * p_zero;
                        }
                }

                // We truncate the edge histogram on the same value as the node
                // histogram. Even when removing all k-mers with coverage < T, a few
                // arcs might still have a coverage < T. We do not want keep those
                // in our histogram as this will interfere with the EM algorithm.
                unsigned int smallSupport = (dBG.getAbundanceMin() < 0)?
                errorHist.begin()->first : dBG.getAbundanceMin();
                for (unsigned int k = 0; k < smallSupport; k++)
                        errorHist.erase(k);

                /*cout << "Current EDGE error hist: \n";
                 f or (unsigned int k = 0; k < 30; k++)               *
                 cout << k << "\t" << errorHist[k] << "\n";*/
                
                // The below EM procedure might not convergence, but the obtained NB
                // parameters should be sensible in all cases: only small support
                // values are inferred, therefore, the mean should never be larger than
                // the mean of the error histogram (and thus be close to zero).
                int maxIter = 10.0 / epsilon;
                int nIter = Util::fitTruncNegBinomEM(errorHist, edgeErrorLambda,
                                                edgeErrorODF, edgeErrorWeight,
                                                epsilon, maxIter);


                if (nIter > maxIter)
                        cout << "\tWARNING: EM algorithm to fit edge error model "
                                "did not converge\n";
                else
                        cout << "\tEM algorithm to fit edge error model converged after " << nIter << " iterations" << endl;
        }

        // Compute edge average and weights
        vector<double> wE(edgeCovModel.getK(), PC);
        wE[0] = edgeErrorWeight;

        vector<double> edgeLambdas(edgeCovModel.getNumStrains(),0.0);
	vector<double> mtmE(edgeCovModel.getNumStrains()*edgeCovModel.getNumStrains(), 0.0);
	vector<double> mtXE(edgeCovModel.getNumStrains(), 0.0);

        // for each (1,0,0), (0,1,0), (0,0,1) add PC, make sure mtm is positive definite
        for (int strain = 0; strain < edgeCovModel.getNumStrains(); strain++){
            mtmE[strain * (edgeCovModel.getNumStrains()+1)] += 1.0;
            mtXE[strain] += edgeCovModel.getFraction(strain) * edgeCovModel.getSosLambda(); //prev lambda
            wE[pow(edgeCovModel.getC(), strain)] += 1.0;
        }

        //vector<double> totWEdgeCov(edgeCovModel.getK(),0.0);
        for (const auto& it : edgeMult) {
                Arc& arc = dBG.getArc(dBG.getArcID(it.first)); // TODO
                vector<int> expMult = it.second.getExpMult();
                double cov = arc.getCov();

                if (accumulate(expMult.begin(), expMult.end(), 0) == 0 || cov > (maxMult) * edgeCovModel.getSosLambda()) //error edge
                    continue;

                // soft assignments
		Assignment mult(vector<int>(edgeCovModel.getNumStrains(), edgeCovModel.getC()));
                for (int idx=1; idx < edgeCovModel.getK(); idx ++){ //don't compute wE[0], already computed
                    mult.increment();
                    wE[idx] += it.second[mult]; //soft assignment
                }

                if (*max_element(expMult.begin(), expMult.end()) > 1)
                        continue;

                // fill mtm and X for least squares estimates
                for (size_t strain = 0; strain < edgeCovModel.getNumStrains(); strain++) {
                        mtXE[strain] += expMult[strain] * cov;
                        for (size_t j = strain; j < edgeCovModel.getNumStrains(); j++)
                                mtmE[strain*edgeCovModel.getNumStrains() + j] += expMult[strain] * expMult[j];

                }
        }

        for (size_t j = 0; j < edgeCovModel.getNumStrains()-1; j++){
            for (size_t i = j+1; i < edgeCovModel.getNumStrains(); i++) {
                mtmE[i*edgeCovModel.getNumStrains() + j] = mtmE[j*edgeCovModel.getNumStrains() + i];
            }
        }

        Util::choleskyLinSolve(mtmE,mtXE);

        /*cout << "EdgeLambdas based on linear system of the observations: " << "[ ";
        for (size_t strain=0; strain < edgeCovModel.getNumStrains(); strain++){
            cout << mtXE[strain] << " ";
        }
        cout  <<  "]" << endl;*/

        for (int strain = 0; strain < edgeCovModel.getNumStrains(); strain++){
            if (mtXE[strain] < 0){
        	retVal=false;
	        //assign previous smallest lambda
                vector<double> frac = edgeCovModel.getFractions();
                mtXE[strain] = *min_element(frac.begin(), frac.end()) * edgeCovModel.getSosLambda();
                cerr << "Warning: encountered negative lambda in edge LambdaEstimates \n"
                        "set to smallest lambda in previous EM iteration \n";
        		"will stop EM after this iteration possibly less strains are present.\n";
            }
        }

        //Use lambda estimate based on OLS of all observations
        edgeLambdas = mtXE;
	sort(edgeLambdas.begin(), edgeLambdas.end());
        double sumEdgeLambdas = accumulate(edgeLambdas.begin(), edgeLambdas.end(), 0.0);
        // Update lambda_err to avoid it getting too big.
	edgeErrorLambda = min(sumEdgeLambdas*0.025, edgeErrorLambda);


        // ====================================================================
        // Compute the shared model parameters
        // ====================================================================

        for (int strain = 0; strain < nodeCovModel.getNumStrains(); strain++){
            double nodeFrac = nodeLambdas[strain]/sumNodeLambdas;
            double edgeFrac =  edgeLambdas[strain]/sumEdgeLambdas;
            double avgFrac = (nodeFrac+edgeFrac)/2;
            nodeLambdas[strain] = avgFrac * sumNodeLambdas;
            edgeLambdas[strain] = avgFrac * sumEdgeLambdas;
        }
	cout << "Final NodeLambdas: " << "[ ";
        for (size_t strain=0; strain < nodeCovModel.getNumStrains(); strain++){
            cout << nodeLambdas[strain] << " ";
        }
        cout  <<  "]" << endl;

	cout << "Final EdgeLambdas: " << "[ ";
        for (size_t strain=0; strain < edgeCovModel.getNumStrains(); strain++){
            cout << edgeLambdas[strain] << " ";
        }
        cout  <<  "]" << endl;


        // ====================================================================
        // Compute the remaining node parameters
        // ====================================================================

        // Compute node variance and ODF
        double nom = 0.0;
        for (const auto& it : nodeMult) {
                SSNode node = dBG.getSSNode(it.first);
                vector<int> expMult = it.second.getExpMult();
                if (node.getAvgCov()> maxMult* sumNodeLambdas)
                        continue;
                if (*max_element(expMult.begin(), expMult.end()) > 1)
                        continue;

                //soft assignments
                Assignment mult(vector<int>(nodeCovModel.getNumStrains(), nodeCovModel.getC()));
                for (int idx=1; idx < nodeCovModel.getK(); idx ++){ //don't compute errorODF
                        mult.increment();
                        double expCov = 0;
                        for (size_t strain = 0; strain < nodeCovModel.getNumStrains(); strain++) {
                        expCov += mult[strain] * nodeLambdas[strain];
                        }
                        double delta = node.getAvgCov() - expCov;
                        nom += it.second[mult] * delta * delta * node.getMarginalLength()/ expCov;
                }
        }

        double nodeODF = nom/accumulate(next(wN_k.begin()), wN_k.end(), 0);
        nodeODF = max(1.0, nodeODF);

        if (eqW){
                // =====================================================================
                // update WN
                // =====================================================================

                int nCond = 0;
                for (int i = 1; i < nstrains; i++){
                        int comb = Util::factorial(nstrains)/(Util::factorial(i)*Util::factorial(nstrains-i)); //choose i out of nstrains
                        nCond += Util::factorial(comb)/(2*Util::factorial(comb-2)); //choose 2 out of comb
                }
                //create Aw to model relationships between weights
                int nr = pow(2,nstrains)-2;
                int nc = pow(2,nstrains) -2 +nCond;
                double Aw_init[nr*nc] = {};
                //each component equal to itself
                for (int i = 0; i<nr; i++){
                        for (int j = 0; j<nc; j++){
                                if (i==j)
                                        Aw_init[i*nc + j] = 1.0;
                                else
                                        Aw_init[i*nc + j] = 0.0;
                        }
                }
                //extra conditions
                BinAssignment ass1(nstrains);
                int col = pow(2,nstrains)-2;
                for (int i = 0; i < pow(2,nstrains)-2; i++){
                        ass1.increment();
                        for (int sum = 1; sum < nstrains; sum++){
                                if (accumulate(ass1.begin(), ass1.end(),0) == sum){
                                        BinAssignment ass2(ass1);

                                        for (int j = i+1; j < pow(2,nstrains)-2; j++){
                                                ass2.increment();
                                                if (accumulate(ass2.begin(), ass2.end(), 0) == sum){
                                                        Aw_init[i * nc + col] = 1.0;
                                                        Aw_init[j * nc + col] = -1.0;
                                                        col ++;
                                                }
                                        }
                                }
                        }
                }

		for (int m = 1; m < nodeCovModel.getC(); m++) {
                	double Aw[nr * nc];
                	std::copy(Aw_init, Aw_init + nr * nc, Aw);
                	double bw[nc] = {}; //initialize all elements to 0.0
               		BinAssignment ass(nstrains);
                	Assignment Wass(vector<int>(nodeCovModel.getNumStrains(), nodeCovModel.getC()));
                	for (int i = 0; i< nr; i++){
                        	ass.increment();
				vector<int> temp(ass);
				transform(temp.begin(), temp.end(), temp.begin(), [m](int &entry){return m*entry;});
                        	bw[i] = wN[Wass.assignmentToIndex(temp)];
                	}
                	double wkopt;
                	double *work;
                	int info, lwork = -1, nrhs = 1;
                    char c = 'N';
                	dgels_(&c, &nc, &nr, &nrhs, Aw, &nc, bw, &nc, &wkopt, &lwork, &info);
                	lwork = (int) wkopt;
                	work = new double[lwork];
                	dgels_(&c, &nc, &nr, &nrhs, Aw, &nc, bw, &nc, work, &lwork, &info);
                	BinAssignment assOut(nstrains);
                	for (int i = 0; i < nr; i++){
                        	assOut.increment();
                        	size_t idx = 0;
                        	size_t numComp = nodeCovModel.getC();
                        	for (size_t j = 0; j < assOut.size(); j++) {
                                	idx += ((m*assOut[j] < numComp) ? m*assOut[j] : numComp-1) * pow(numComp,j);
                        	}
                        	wN[idx] = max(1.0,bw[i]);
                	}
                	delete [] work;
		}
		
		// ---------------
		// Equalise other weights 
		// ---------------
		if(nodeCovModel.getC() == 3){ // For now only reasoned for maxmult = 2 (extra conditions needed in further cases)
                        // conditions without zero mult (and 1 diff between other mults)
                        for (int m = 1; m < (nodeCovModel.getC()-1); m++) {
                                double Awn[nr * nc];
                                std::copy(Aw_init, Aw_init + nr * nc, Awn);
                                double bwn[nc] = {}; //initialize all elements to 0.0
                                BinAssignment assn(nstrains);
                                Assignment Wassn(vector<int>(nodeCovModel.getNumStrains(), nodeCovModel.getC()));
                                for (int i = 0; i< nr; i++){
                                        assn.increment();
                                        vector<int> temp(assn);
                                        transform(temp.begin(), temp.end(), temp.begin(), [m](int &entry){return m+entry;});
                                        bwn[i] = wN[Wassn.assignmentToIndex(temp)];
                                }
                                double wkopte;
                                double *worke;
                                int infoe, lworke = -1, nrhse = 1;
                                char c = 'N';
                                dgels_(&c, &nc, &nr, &nrhse, Awn, &nc, bwn, &nc, &wkopte, &lworke, &infoe);
                                lworke = (int) wkopte;
                                worke = new double[lworke];
                                dgels_(&c, &nc, &nr, &nrhse, Awn, &nc, bwn, &nc, worke, &lworke, &infoe);
                                BinAssignment assOut(nstrains);
                                for (int i = 0; i < nr; i++){
                                        assOut.increment();
                                        vector<int> temp(assOut);
                                        transform(temp.begin(), temp.end(), temp.begin(), [m](int &entry){return m+entry;});
                                        wN[Wassn.assignmentToIndex(temp)] = max(1.0,bwn[i]);
                                }
                                delete [] worke;
                        }
                        // Condtions with at least one m_i = 0
                        {
                                Assignment tempAss(vector<int>(nodeCovModel.getNumStrains(), nodeCovModel.getC()));
                                for (int n0 = 1; n0 <= nstrains - 2; n0 ++){
                                        for(int n1 = 1; n1 <= nstrains - n0 - 1; n1 ++){
                                                //                                       for(int n2 = 1; n2 < nstrains - n0 - n1; n2 ++){
                                                int nr = Util::factRatio(nstrains,n0) / Util::factorial(nstrains-n0) * Util::factRatio(nstrains-n0, n1)/Util::factorial(nstrains-n0-n1);
                                                int nc = 2*nr;
                                                double A[nr*nc] = {};
                                                double b[nc] = {};
                                                vector<int> ass(nstrains);
                                                fill(ass.begin(),ass.begin()+n0,0);
                                                fill(ass.begin()+n0, ass.begin()+n0+n1,1);
                                                fill(ass.begin()+n0+n1, ass.end(),2);
                                                int i=0;
                                                do {
                                                        A[i*nc + i] = 1.0;
                                                        A[i*nc + (nr + i)] = 1.0;
                                                        A[((i+1)%nr)*nc + (nr + i)] = -1.0;
                                                        b[i] = wN[tempAss.assignmentToIndex(ass)];
                                                        i++;
                                                } while(std::next_permutation(ass.begin(), ass.end()));
                                                //                                        }
                                                double wkopt;
                                                double *work;
                                                int info, lwork = -1, nrhs = 1;
                                                char c = 'N';
                                                dgels_(&c, &nc, &nr, &nrhs, A, &nc, b, &nc, &wkopt, &lwork, &info);
                                                lwork = (int) wkopt;
                                                work = new double[lwork];
                                                dgels_(&c, &nc, &nr, &nrhs, A, &nc, b, &nc, work, &lwork, &info);
                                                std::sort(ass.begin(), ass.end());
                                                i=0;
                                                do {
                                                        wN[tempAss.assignmentToIndex(ass)] = max(1.0,b[i]);
                                                        i++;
                                                } while(std::next_permutation(ass.begin(), ass.end()));
                                                delete [] work;
                                        }
                                }
                        }
                }
        }

	// ====================================================================
        // Compute the remaining edge parameters
        // ====================================================================

        // Compute edge variance and ODF
        nom = 0.0;
        for (const auto& it : edgeMult) {
                Arc& arc = dBG.getArc(dBG.getArcID(it.first)) ;
                vector<int> expMult = it.second.getExpMult();
		if(arc.getCov() > maxMult * sumEdgeLambdas)
	               	continue;
                if (*max_element(expMult.begin(), expMult.end()) > 1)
                        continue;

                //soft assignments
                Assignment mult(vector<int>(edgeCovModel.getNumStrains(), edgeCovModel.getC()));
                for (int idx=1; idx < edgeCovModel.getK(); idx ++){ //don't compute errorODF
                    mult.increment();
                    double expCov = 0;
                    for (size_t strain = 0; strain < edgeCovModel.getNumStrains(); strain++) {
                        expCov += mult[strain] * edgeLambdas[strain];
                    }
                    double delta = arc.getCov() - expCov;
                    nom += it.second[mult] * delta * delta / expCov;
                }
        }

        double edgeODF = nom / accumulate(next(wE.begin()), wE.end(), 0);
        edgeODF = max(1.0, edgeODF);

        if(eqW){
                // ========================================================
                // OLS to equalise weights
                // =======================================================
                int nCond = 0;
                for (int i = 1; i < nstrains; i++){
                        int comb = Util::factorial(nstrains)/(Util::factorial(i)*Util::factorial(nstrains-i)); //choose i out of nstrains
                        nCond += Util::factorial(comb)/(2*Util::factorial(comb-2)); //choose 2 out of comb
                }
                //create Aw to model relationships between weights
                int nr = pow(2,nstrains)-2;
                int nc = pow(2,nstrains) -2 +nCond;
                double Aw_init[nr*nc] = {};
                //each component equal to itself
                for (int i = 0; i<nr; i++){
                        for (int j = 0; j<nc; j++){
                                if (i==j)
                                        Aw_init[i*nc + j] = 1.0;
                                else
                                        Aw_init[i*nc + j] = 0.0;
                        }
                }
                //extra conditions
                BinAssignment ass1(nstrains);
                int col = pow(2,nstrains)-2;
                for (int i = 0; i < pow(2,nstrains)-2; i++){
                        ass1.increment();
                        for (int sum = 1; sum < nstrains; sum++){
                                if (accumulate(ass1.begin(), ass1.end(),0) == sum){
                                        BinAssignment ass2(ass1);

                                        for (int j = i+1; j < pow(2,nstrains)-2; j++){
                                                ass2.increment();
                                                if (accumulate(ass2.begin(), ass2.end(), 0) == sum){
                                                        Aw_init[i * nc + col] = 1.0;
                                                        Aw_init[j * nc + col] = -1.0;
                                                        col ++;
                                                }
                                        }
                                }
                        }
                }
                for (int m = 1; m < edgeCovModel.getC(); m++) {
                	double Awe[nr * nc];
                	std::copy(Aw_init, Aw_init + nr * nc, Awe);
                	double bwe[nc] = {}; //initialize all elements to 0.0
               		BinAssignment asse(nstrains);
                	Assignment Wasse(vector<int>(edgeCovModel.getNumStrains(), edgeCovModel.getC()));
                	for (int i = 0; i< nr; i++){
                        	asse.increment();
				vector<int> temp(asse);
				transform(temp.begin(), temp.end(), temp.begin(), [m](int &entry){return m*entry;});
                        	bwe[i] = wE[Wasse.assignmentToIndex(temp)];
                	}
                	double wkopte;
                	double *worke;
                	int infoe, lworke = -1, nrhse = 1;
                    char c = 'N';
                	dgels_(&c, &nc, &nr, &nrhse, Awe, &nc, bwe, &nc, &wkopte, &lworke, &infoe);
                	lworke = (int) wkopte;
                	worke = new double[lworke];
                	dgels_(&c, &nc, &nr, &nrhse, Awe, &nc, bwe, &nc, worke, &lworke, &infoe);
                	BinAssignment assOut(nstrains);
                	for (int i = 0; i < nr; i++){
                        	assOut.increment();
                        	size_t idx = 0;
                        	size_t numCompe = edgeCovModel.getC();
                        	for (size_t j = 0; j < assOut.size(); j++) {
                                	idx += ((m*assOut[j] < numCompe) ? m*assOut[j] : numCompe-1) * pow(numCompe,j);
                        	}
                        	wE[idx] = max(1.0,bwe[i]);
                	}
                	delete [] worke;
		}
		// ---------------
		// Equalise other weights 
		// ---------------
		if(edgeCovModel.getC() == 3){ // For now only reasoned for maxmult = 2 (extra conditions needed in further cases)
                        // conditions without zero mult (and 1 diff between other mults)
                        for (int m = 1; m < (edgeCovModel.getC()-1); m++) {
                                double Awe[nr * nc];
                                std::copy(Aw_init, Aw_init + nr * nc, Awe);
                                double bwe[nc] = {}; //initialize all elements to 0.0
                                BinAssignment asse(nstrains);
                                Assignment Wasse(vector<int>(edgeCovModel.getNumStrains(), edgeCovModel.getC()));
                                for (int i = 0; i< nr; i++){
                                        asse.increment();
                                        vector<int> temp(asse);
                                        transform(temp.begin(), temp.end(), temp.begin(), [m](int &entry){return m+entry;});
                                        bwe[i] = wE[Wasse.assignmentToIndex(temp)];
                                }
                                double wkopte;
                                double *worke;
                                int infoe, lworke = -1, nrhse = 1;
                                char c = 'N';
                                dgels_(&c, &nc, &nr, &nrhse, Awe, &nc, bwe, &nc, &wkopte, &lworke, &infoe);
                                lworke = (int) wkopte;
                                worke = new double[lworke];
                                dgels_(&c, &nc, &nr, &nrhse, Awe, &nc, bwe, &nc, worke, &lworke, &infoe);
                                BinAssignment assOut(nstrains);
                                for (int i = 0; i < nr; i++){
                                        assOut.increment();
                                        //size_t idx = 0;
                                        //size_t numCompe = edgeCovModel.getC();
                                        //for (size_t j = 0; j < assOut.size(); j++) {
                                        //        idx += ((m*assOut[j] < numCompe) ? m*assOut[j] : numCompe-1) * pow(numCompe,j);
                                        //}
                                        vector<int> temp(assOut);
                                        transform(temp.begin(), temp.end(), temp.begin(), [m](int &entry){return m+entry;});
                                        wE[Wasse.assignmentToIndex(temp)] = max(1.0,bwe[i]);
                                        //wE[idx] = max(1.0,bwe[i]);
                                }
                                delete [] worke;
                        }
                        // Condtions with at least one m_i = 0
                        {
                        Assignment tempAss(vector<int>(edgeCovModel.getNumStrains(), edgeCovModel.getC()));
                        for (int n0 = 1; n0 <= nstrains - 2; n0 ++){
                                for(int n1 = 1; n1 <= nstrains - n0 - 1; n1 ++){
//                                       for(int n2 = 1; n2 < nstrains - n0 - n1; n2 ++){
                                                int nr = Util::factRatio(nstrains,n0) / Util::factorial(nstrains-n0) * Util::factRatio(nstrains-n0, n1)/Util::factorial(nstrains-n0-n1);
                                                int nc = 2*nr;
                                                double A[nr*nc] = {};
                                                double b[nc] = {};
                                                vector<int> ass(nstrains);
                                                fill(ass.begin(),ass.begin()+n0,0);
                                                fill(ass.begin()+n0, ass.begin()+n0+n1,1);
                                                fill(ass.begin()+n0+n1, ass.end(),2);
                                                int i=0;
                                                do {
                                                        A[i*nc + i] = 1.0;
                                                        A[i*nc + (nr + i)] = 1.0;
                                                        A[((i+1)%nr)*nc + (nr + i)] = -1.0;
                                                        b[i] = wE[tempAss.assignmentToIndex(ass)];
                                                        i++;
                                                } while(std::next_permutation(ass.begin(), ass.end()));
//                                        }
                                                double wkopt;
                                                double *work;
                                                int info, lwork = -1, nrhs = 1;
                                                char c = 'N';
                                                dgels_(&c, &nc, &nr, &nrhs, A, &nc, b, &nc, &wkopt, &lwork, &info);
                                                lwork = (int) wkopt;
                                                work = new double[lwork];
                                                dgels_(&c, &nc, &nr, &nrhs, A, &nc, b, &nc, work, &lwork, &info);
                                                std::sort(ass.begin(), ass.end());
                                                i=0;
                                                do {
                                                        wE[tempAss.assignmentToIndex(ass)] = max(1.0,b[i]);
                                                        i++;
                                                } while(std::next_permutation(ass.begin(), ass.end()));
                                                delete [] work;
                                }
                        }
                        }
                }

        }

	//=======================================================================
	// Only keep weights for multiplicities with 0 and one mult
	// Derive the rest as (small) fractions
	// =====================================================================

	/*Assignment ass(vector<int>(nodeCovModel.getNumStrains(), nodeCovModel.getC()));
        int numComp = nodeCovModel.getC();
        for (int i = 0; i < pow(numComp,nodeCovModel.getNumStrains()); i++){
                auto mm = minmax_element(ass.begin(), ass.end());
                if (*mm.first == *mm.second){
			ass.increment();
			continue;
		}else {
                        vector<int> test(ass);
                        sort(test.begin(), test.end());
                        std::vector<int>::iterator it = unique_copy(test.begin(), test.end(), test.begin());
                        it = remove_if(test.begin(), it, [](int a){return a==0;});
                        test.resize(distance(test.begin(), it));
                        if(test.size() != 1){
                                vector<int> temp(ass);
				for (int j=0; j< temp.size(); j++)
					temp[j] = (temp[j]==0)? 0 : *mm.second;
				wN[i] = wN[ass.assignmentToIndex(temp)]*pow(10, -1 * log10(nodeMult.size())/2); //TODO validate and/or adapt this value!
				wE[i] = wE[ass.assignmentToIndex(temp)]*pow(10, -1 * log10(edgeMult.size())/2);
			}
                }
                ass.increment();
        }*/



        //======================================================================
        // Make new covmodel
        //======================================================================
        nodeCovModel = CovModel(nodeErrorLambda, nodeErrorODF,
                                nodeLambdas, nodeODF, wN, nodeCovModel.getMaxMult());

        Assignment mult(vector<int>(nodeCovModel.getNumStrains(), Multiplicity::maxM+1));
        double maxLProb = std::numeric_limits<double>::lowest();
        for (int i = 0; i < nodeCovModel.getK(); i++) {
            if (find(mult.begin(), mult.end(), Multiplicity::maxM) == mult.end()){
                mult.increment();
                continue;
            }
            double multP = nodeCovModel.getLogProb(((double)Multiplicity::maxM) * nodeCovModel.getSosLambda(), mult);
            if (multP > maxLProb)
                maxLProb = multP;
            mult.increment();
        }

        nodeCovModel.setMaxMultP(maxLProb);

        cout << "estimated maxMProb: " << maxLProb << endl;

        //cout << "estimated ODF: " << nodeErrorODF << "(err), " << nodeODF << "(true)" << endl;

        cout << "\tNode spectrum: " << nodeCovModel << endl;


        edgeCovModel = CovModel(edgeErrorLambda, edgeErrorODF,
				edgeLambdas, edgeODF, wE, edgeCovModel.getMaxMult());

        Assignment edgeM(vector<int>(edgeCovModel.getNumStrains(), Multiplicity::maxM+1));
        double maxLProbE = std::numeric_limits<double>::lowest();;
        for (int i = 0; i < edgeCovModel.getK(); i++) {
            if (find(edgeM.begin(), edgeM.end(), Multiplicity::maxM) == edgeM.end()){
                edgeM.increment();
                continue;
            }
            double multP = edgeCovModel.getLogProb(Multiplicity::maxM * edgeCovModel.getSosLambda(), edgeM);
            if (multP > maxLProbE)
                maxLProbE = multP;
            edgeM.increment();
        }
        edgeCovModel.setMaxMultP(maxLProbE);

	cout << "estimated maxMProb: " << maxLProbE << endl;

        //cout << "estimated ODF: " << edgeErrorODF << "(err), " << edgeODF << "(true)" << endl;

        cout << "\tEdge spectrum: " << edgeCovModel << endl;

	if (nodeLambdas[0] < 1 || edgeLambdas[0] < 1){
                cout << "smallest strain coverage < 1, stopping EM" << endl;
                retVal=false;
        }

        return retVal;
}

void CRFMult::MstepWithTrueMult(const std::vector<NodeRep>& nodes,
                               CovModel& nodeCovModel,
                               const std::vector<EdgeRep>& edges,
                               CovModel& edgeCovModel, double epsilon)
{
        NodeMap<Multiplicity> nodeMult(nodes.size()); EdgeMap<Multiplicity> edgeMult(edges.size());
        for (int i=0; i < nodes.size(); i++)
                nodeMult[nodes[i]] = Multiplicity(dBG.getTrueNodeMult(nodes[i]));

        for (int j = 0; j < edges.size(); j++)
                edgeMult[edges[j]] = Multiplicity(dBG.getTrueEdgeMult(edges[j]));

        Mstep(nodeMult, nodeCovModel,
              edgeMult, edgeCovModel, epsilon,false);
}


int CRFMult::computeMultEM(NodeMap<Multiplicity>& nodeMult,
                           CovModel& nodeCovModel,
                           EdgeMap<Multiplicity>& edgeMult,
                           CovModel& edgeCovModel,
                           double epsilon, double maxChange, int maxIter,
                           bool eqW,
                           bool approxInf, bool map, bool fixedZero)
{
        int iter;
        double logZ = -1.0;
        for (iter = 1; iter <= maxIter; iter++) {
                cout << "Iteration " << iter << endl;
                // infer multiplicities given the model
                if (!Estep(nodeMult, nodeCovModel,
                           edgeMult, edgeCovModel, epsilon, maxChange,
                           logZ, approxInf, map)) 
                        return iter;
                cout << "logZ in E-step: " << logZ << endl;
                // update the model given the multiplicities
                if (!Mstep(nodeMult, nodeCovModel,
                           edgeMult, edgeCovModel, max(epsilon, maxChange), eqW, fixedZero)) 
                        return iter;

        }

        return (iter);
}
