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

#ifndef CRFMULT_H
#define CRFMULT_H

#include <vector>
#include <cstdlib>
#include <map>
#include <queue>
#include <atomic>

#include "coverage.h"
#include "pgm/factor.h"
#include "dbgraph.h"
#include "bitvec.h"
#include "libdai/dai/bp.h"
#include "libdai/dai/parallelbp.h"
#include "libdai/dai/factorgraph.h"
#include "libdai/dai/factor.h"

// ============================================================================
// CLASS PROTOTYPES
// ============================================================================

class DBGraph;
class WorkLoadBalancer;
class NodeRep;
class EdgeRep;

// ============================================================================
// DIJKSTRA AUXILIARY CLASSES
// ============================================================================

class PathDFS {

public:
        NodeID nodeID;          // current node identifier
        NodeID prevID;          // previous node identifier
        size_t length;          // length to current node

        /**
         * Default constructor
         * @param nodeID node identifier
         * @param prevID previous node identifier
         * @param length length to the current node
         */
        PathDFS(NodeID nodeID, NodeID prevID, size_t length) :
                nodeID(nodeID), prevID(prevID), length(length) {};
};

struct PathDFSComp {

        /**
         * Compare two paths (by length) (strict weak ordering (irreflexive))
         */
        bool operator()(const PathDFS& f, const PathDFS& s) {
                return f.length > s.length;
        }
};

// ============================================================================
// CRF PERFORMANCE COUNTER
// ============================================================================

class CRFPerfCounter {

public:
        size_t totNumCRF;
        size_t totNumNodes;
        size_t totNumEdges;
        size_t totActNBSize;
        size_t totMaxCRFFactor;

        /**
         * Initialization constructor
         */
        CRFPerfCounter() {
                reset();
        }

        /**
         * Reset the counter
         */
        void reset() {
                totNumCRF = totNumNodes = totNumEdges = 0;
                totActNBSize = totMaxCRFFactor = 0;
        }

        /**
         * Merge counts
         * @param rhs Right hand size
         **/
        void merge(const CRFPerfCounter& rhs) {
                totNumCRF += rhs.totNumCRF;
                totNumNodes += rhs.totNumNodes;
                totNumEdges += rhs.totNumEdges;
                totActNBSize += rhs.totActNBSize;
                totMaxCRFFactor += rhs.totMaxCRFFactor;
        }

        /**
         * Get the number of CRF solved
         * @return The number of CRF solved
         */
        size_t getNumCRF() const {
                return totNumCRF;
        }

        /**
         * Get the average number of nodes in a CRF neighborhood
         * @return The average number of edges in a CRF neighborhood
         **/
        double getAvgNumNodes() const {
                return (double)totNumNodes / totNumCRF;
        }

        /**
         * Get the average number of edges in a CRF neighborhood
         * @return The average number of edges in a CRF neighborhood
         **/
        double getAvgNumEdges() const {
                return (double)totNumEdges / totNumCRF;
        }

        /**
         * Get the average actual CRF neighborhood size
         * @return The average actual CRF neighborhood size
         **/
        double getAvgActNBSize() const {
                return (double)totActNBSize / totNumCRF;
        }

        /**
         * Get the average maximum CRF factor size
         * @return The average maximum CRF factor size
         **/
        double getAvgMaxCRFFactor() const {
                return (double)totMaxCRFFactor / totNumCRF;
        }
};

// ============================================================================
// CONDITIONAL RANDOM FIELDS SOLVER (PER THREAD)
// ============================================================================

/**
 * This class is designed for speed. Most of the containers are allocated as
 * class members instead of local variables to avoid repeated (de)allocation
 * of these containers over many multiplicity computations.
 * Therefore, some functions may have side effects: they modify class members
 * containers. This class is therefore not thread-safe. The idea is to run
 * a single class object within each thread.
 */

class CRFSolver {

private:
        const DBGraph& dBG;             // const-ref to de Bruijn graph
        

        //int multMargin;                 // number of alt mult (one-sided)
        size_t maxFactorSize;           // maximum size of intermediate factor
        double flowStrength;            // value in flow-cpd when mults do not agree

        // containers below are class members to allow their reuse over many
        // multiplicity computations, thus avoiding repeated (de)allocation
        std::vector<NodeRep> nodes;     // nodes in subgraph
        std::vector<EdgeRep> edges;     // edges in subgraph
        std::priority_queue<NodeDFS, std::vector<NodeDFS>, NodeDFSComp> pq;
        Bitvec bitvec;
        
        double logZ;

        /**
         * Get a singleton factor with a multinomial over multiplicities (1 singleton factor over all strain-CRF-nodes
         * @param varID Variable identifiers
         * @param card Variable cardinalities
         * @param firstMults Multiplicity of the first value for each strain
         * @param covModel Coverage model
         * @param numObs Number of observations
         * @param numIndepOb FIXME !!!
         * @return Singleton multiplicity factor
         */
        static Factor createSingletonFactor(std::vector<int> varIDs, std::vector<int> cards,
                                            const CovModel& covModel,
                                            const double numObs, const int numIndepObs);
        
        static dai::Factor createLDSingletonFactor(std::vector<int> varIDs, 
                                                   std::vector<int> varCards,
                                                   const CovModel& covModel,
                                                   const double numObs, const int numIndepObs);
        
        /**
         * Create a uniform singleton factor 
         * such that each node's/arc's variables will be in a factor together
         * TODO: is this really the best way to extract?
         */
        static dai::Factor createLDUniformSingletonFactor(std::vector<int> varIDs,
                                                          std::vector<int> varCards);

        /**
         * Get a flow-conservation factor
         * @param sumVarID Variable identifier of the sum
         * @param sumCard Cardinality of the sum variable
         * @param sumFirstMult Multiplicity of the first sum value
         * @param termVarID Variable identifiers of the terms
         * @param termCard Cardinality of the term variables
         * @param termFirstMult Multiplicity of the first term values
         * @param palindromic Is the edge palindromic
         * @return Flow-conservation factor
         */
        static Factor createFlowFactor(int sumVarID, int sumCard,
                                       const std::vector<int>& termVarID,
                                       const std::vector<int>& termCard,
                                       const std::vector<bool>& palindromic,
                                       double flowStrength);
        

        static dai::Factor createLDFlowFactor(int sumVarID, int sumCard,
                                       const std::vector<int>& termVarID,
                                       const std::vector<int>& termCard,
                                       const std::vector<bool>& palindromic,
                                       double flowStrength);
        

        /**
         * Internal function to find a CRF-compliant subgraph.
         * Output is written to nodes and edges as a side-effect.
         * @param maxDepth Maximum depth (in terms of number of nodes)
         */
        void getSubgraph(int maxDepth);

        /**
         * Get a subgraph of the de Bruijn graph centered around a seed node
         * Output is written to nodes and edges as a side-effect
         * @param seedNode Seed node representative
         * @param maxDepth Maximum depth (in terms of number of nodes)
         */
        void getSubgraph(NodeRep seedNode, int maxDepth);

        /**
         * Get a subgraph of the de Bruijn graph centered around a seed edge
         * Output is written to nodes and edges as a side-effect
         * @param seedNode Seed node representative
         * @param maxDepth Maximum depth (in terms of number of nodes)
         */
        void getSubgraph(EdgeRep seedEdge, int maxDepth);

        /**
        * Compute the multiplicity using a CRF model for multiple strains
        * @param node2var Mapping between dBG nodes and CRF variables
        * @param edge2var Mapping between dBG edges and CRF variables
        * @param nodeCovModel Node coverage model
        * @param edgeCovModel Node coverage model
        * @param resMult Resulting multiplicity (output)
        * @param targetCRFVar CRF variable to keep
        */
        bool solveSubgraph(NodeMap<std::vector<int>>& node2var,
                           EdgeMap<std::vector<int>>& edge2var,
                           const CovModel& nodeCovModel,
                           const CovModel& edgeCovModel,
                           Multiplicity& resMult, std::vector<int> targetCRFVar,
                           size_t& factorSize) const;
                           
        /*void retrieveNodeMultFromFG(WorkLoadBalancer& wlb, 
                                const std::map<int, dai::Var>& label2var,
                                const NodeMap<std::vector<int>> node2var,
                                const std::vector<NodeRep> nodes,
                                const dai::BP& bp,
                                NodeMap<Multiplicity>& nodeMult,
                                const CovModel& nodeCovModel) const;
                                                       
                                                       
        void retrieveEdgeMultFromFG(WorkLoadBalancer& wlb, 
                                        const std::map<int, dai::Var>& label2var,
                                        const EdgeMap<std::vector<int>> edge2var,
                                        const std::vector<EdgeRep> edges,
                                        const dai::BP& bp,
                                        EdgeMap<Multiplicity>& edgeMult,
                                        const CovModel& edgeCovModel) const;
                                                                                   
        Multiplicity getMultFromFG(const std::vector<int>& labels,
                                const std::map<int, dai::Var>& label2var,
                                const dai::BP& bp,
                                const CovModel& cm) const;*/
                                
        Multiplicity getMultFromFG(const std::vector<int>& labels,
                                const std::map<int, dai::Var>& label2var,
                                const size_t factorID,
                                const dai::BP& bp,
                                const CovModel& cm) const;
                                
        Multiplicity getMultFromFG(const std::vector<int>& labels,
                                        const std::map<int, dai::Var>& label2var,
                                        const size_t factorID,
                                        const dai::ParallelBP& bp,
                                        const CovModel& cm) const;

public:
        CRFPerfCounter perfCounter;     // performance counter

        /**
         * Constructor
         * @param dBG Const-reference to de Bruijn graph
         * @param multMargin Number of alternative multiplicities (one-sided)
         * @param maxFactorSize Maximum size of intermediate factor
         * @param numThreads Number of threads
         * @param threadWork Amount of work (nodes/edges) per work chunk
         */
        CRFSolver(const DBGraph& dBG, size_t maxFactorSize,
                  double flowStrength) : dBG(dBG),
                  maxFactorSize(maxFactorSize), flowStrength(flowStrength),
                  bitvec(dBG.getNumNodes()+1), logZ(-1.0) {}
                  
        CRFSolver(const DBGraph& dBG, size_t maxFactorSize,
                  double flowStrength, const std::vector<NodeRep>& n, const std::vector<EdgeRep>& e) : dBG(dBG),
                            maxFactorSize(maxFactorSize), flowStrength(flowStrength), nodes(n), edges(e),
                            bitvec(dBG.getNumNodes()+1), logZ(-1.0) {}
                            
        size_t getNodesSize(){
                return nodes.size();
        }
        size_t getEdgesSize(){
                return edges.size();
        }
                                            

        /**
         * Compute the node multiplicity using a CRF model
         * Check whether flow is conserved based on the prior beliefs
         * TODO Adapt when we want to use this again +
         * TODO is this still valuable if CovModels's getExpMult does not calculate everything?
         * @param node Node
         * @param nodeCovModel Node coverage model
         * @param edgeCovModel Node coverage model
         * @param graphDepth CRF graph depth (in terms of number of nodes)
         * @return true of false (when flow is not conserved)
         */
        bool checkFlow(NodeRep node, const CovModel& nodeCovModel,
                       const CovModel& edgeCovModel, int graphDepth);

        /**
         * Compute the node multiplicity using a CRF model
         * @param node Node
         * @param nodeCovModel Node coverage model
         * @param edgeCovModel Node coverage model
         * @param graphDepth CRF graph depth (in terms of number of nodes)
         */
        Multiplicity getNodeMultiplicity(NodeRep node,
                                         const CovModel& nodeCovModel,
                                         const CovModel& edgeCovModel,
                                         int graphDepth);

        /**
         * Compute the edge multiplicity using a CRF model
         * @param edge Edge
         * @param nodeCovModel Node coverage model
         * @param edgeCovModel Node coverage model
         * @param graphDepth CRF graph depth (in terms of number of nodes)
         */
        Multiplicity getEdgeMultiplicity(EdgeRep edge,
                                         const CovModel& nodeCovModel,
                                         const CovModel& edgeCovModel,
                                         int graphDepth);
        
        /**
         * Create a factor graph for the full dBG
         * 
         */
        dai::FactorGraph getLibdaiFG(NodeMap<std::vector<int> >& node2var, 
                                     EdgeMap<std::vector<int> >& edge2var, 
                                     NodeMap<size_t>& node2fact,
                                     EdgeMap<size_t>& edge2fact,
                                     const CovModel& nodeCovModel, const CovModel& edgeCovModel);
        
	/*
	 * Approximate inference on a subgraph to be used for the cytoscape graph
	 * Pass along the nodes and edges ID vector to determine correct order for node and edge Mult
	 * because we need the outer edges for the crf, the subgraph will be extracted anew: 
	 * thus central node and required graph depth need to be passed again.
	 */
	void approxSubgraphMult(NodeRep node,
				NodeMap<Multiplicity>& nodeMult,
				EdgeMap<Multiplicity>& edgeMult,
				const CovModel& nodeCovModel,
				const CovModel& edgeCovModel,
                         int graphDepth, bool MAP, size_t numThreads, size_t threadWork);
        /**
         * Get subgraph from the de Bruijn graph 
         * that contains all nodes reachable from seedNode
         * @param seedNode Seed node representative
         */
        std::vector<NodeID> getFullConnectedComponent(NodeRep seedNode);
        
        std::vector<EdgeRep> getEdgeReps(){
                return edges;
        }
        
        /**
         * Use approximate inference to compute te multiplicities of the passed nodes
         * and arcs. Assumes that these form a connected subgraph
         */
        void approxMult(NodeMap<Multiplicity>& nodeMult,
                        const CovModel& nodeCovModel,
                        EdgeMap<Multiplicity>& edgeMult,
                        const CovModel& edgeCovModel, 
                        std::string modelCount,
                        bool MAP, size_t numThreads, size_t threadWork);
        
        double getCurrentLogZ() {
                return logZ;
        }
};

// ============================================================================
// ESTIMATE NODE/EDGE MULTIPLICITIES USING CONDITIONAL RANDOM FIELDS
// ============================================================================

class CRFMult {

private:
        const DBGraph& dBG;             // const-ref to de Bruijn graph
        
        bool approximate;               // use approximate inference on whole dBG / use exact inference on subgraphs
        int maxGraphDepth;              // maximum depth of the subgraph
        //int multMargin;                 // number of alt mult (one-sided)
        //int maxMult;			// maximum multiplicity in singleton factors
        size_t maxFactorSize;           // maximum size of intermediate factor
        double flowStrength;            // value in flow-cpd when mults do not agree

        size_t numThreads;              // number of threads
        size_t threadWork;              // amount of work per thread

        mutable CRFPerfCounter totPerfCounter;  // total performance counter

        /**
         * Compute the node conservation of flow
         * @param solver Per-thread solver
         * @param wlb Reference to the workload balancer
         * @param nodes Nodes to handle
         * @param nodeCovModel Node coverage model
         * @param edgeCovModel Node coverage model
         * @param flowOK Flag to indicate whether or not the flow is OK
         */
        void checkNodeFlow(CRFSolver& solver,
                           WorkLoadBalancer& wlb,
                           const std::vector<NodeRep>& nodes,
                           const CovModel& nodeCovModel,
                           const CovModel& edgeCovModel,
                           std::vector<bool>& flowOK) const;

        /**
         * Compute the node multiplicities (thread entry)
         * @param wlb Reference to the workload balancer
         * @param nodes Nodes to handle
         * @param nodeCovModel Node coverage model
         * @param edgeCovModel Node coverage model
         * @param nodeMult Estimated node multiplicity (output)
         */
        void computeNodeMult(CRFSolver& solver,
                             WorkLoadBalancer& wlb,
                             const std::vector<NodeRep>& nodes,
                             const CovModel& nodeCovModel,
                             const CovModel& edgeCovModel,
                             NodeMap<Multiplicity>& nodeMult) const;

        /**
         * Compute the edge multiplicities (thread entry)
         * @param wlb Reference to the workload balancer
         * @param edges Edges to handle
         * @param nodeCovModel Node coverage model
         * @param edgeCovModel Node coverage model
         * @param edgeMult Estimated edge multiplicity (output)
         */
        void computeEdgeMult(CRFSolver& solver,
                             WorkLoadBalancer& wlb,
                             const std::vector<EdgeRep>& edges,
                             const CovModel& nodeCovModel,
                             const CovModel& edgeCovModel,
                             EdgeMap<Multiplicity>& edgeMult) const;
                             

        /**
         * Compute the node/edge multiplicities and check for changes
         * @param nodeMult Estimated node multiplicity (input/output)
         * @param nodeCovModel Node coverage model
         * @param edgeMult Estimated edge multiplicity (input/output)
         * @param edgeCovModel Edge coverage model
         * @param epsilon Tolerance on fraction of altered nodes/edges
         * @return True if a fraction of > epsilon nodes or edges changed
         */
        bool Estep(NodeMap<Multiplicity>& nodeMult,
                   const CovModel& nodeCovModel,
                   EdgeMap<Multiplicity>& edgeMult,
                   const CovModel& edgeCovModel, 
                   double epsilon, double converged, double& prevLogZ,
                   bool approxInf = false,
                   bool map = false) const;


        /**
         * Update the node/edge models given the multiplicities
         * @param nodeMult Estimated node multiplicity
         * @param nodeCovModel Node coverage model (input/output)
         * @param edgeMult Estimated edge multiplicity
         * @param edgeCovModel Edge coverage model (input/output)
         * @param epsilon Tolerance on the lambda_1 estimate
         * @return True if the lambda_1 component change is > epsilon
         */
        bool Mstep(const NodeMap<Multiplicity>& nodeMult,
                   CovModel& nodeCovModel,
                   const EdgeMap<Multiplicity>& edgeMult,
                   CovModel& edgeCovModel, double epsilon, 
                   bool eqW, bool fixedZero = false) const;
public:
        /**
         * Constructor
         * @param dBG Const-reference to de Bruijn graph
         * @param maxGraphDepth Maximum depth of a subgraph
         * @param multMargin Number of alternative multiplicities (one-sided)
         * @param maxFactorSize Maximum size of intermediate factor
         * @param numThreads Number of threads
         * @param threadWork Amount of work (nodes/edges) per work chunk
         */
        CRFMult(const DBGraph& dBG, int maxGraphDepth,
                size_t maxFactorSize, double flowStrength,
                bool approximate = false,
                size_t numThreads = 1, size_t threadWork = 1000) : dBG(dBG), maxGraphDepth(maxGraphDepth), maxFactorSize(maxFactorSize), flowStrength(flowStrength),
                approximate(approximate), numThreads(numThreads), threadWork(threadWork) {}

        /**
         * Check whether the flow around a set of nodes is OK
         * @param nodes Nodes to handle
         * @param flowOK Flag to indicate whether or not the flow is OK
         * @param nodeCovModel Node coverage model
         * @param edgeCovModel Edge coverage model
         */
        void checkFlow(const std::vector<NodeRep>& nodes,
                       std::vector<bool>& flowOK,
                       const CovModel& nodeCovModel,
                       const CovModel& edgeCovModel) const;

        /**
         * Compute the normalized node/edge multiplicities using exact inference
         * @param nodeMult Estimated node multiplicity (input/output)
         * @param nodeCovModel Node coverage model
         * @param edgeMult Estimated edge multiplicity (input/output)
         * @param edgeCovModel Edge coverage model
         */
        void computeMult(NodeMap<Multiplicity>& nodeMult,
                         const CovModel& nodeCovModel,
                         EdgeMap<Multiplicity>& edgeMult,
                         const CovModel& edgeCovModel) const;

        /**
         * @param nodeMult entries to save the node multiplicities
         * @param nodeCovModel coverage model for the nodes
         * @param edgeMult entries to save the edge multiplicities
         * @param edgeCovModel coverage model for the edges
         * @param logZ logZ of the CRF model
         * Uses approximate inference to compute all multiplicities in the dBG 
         * Only saves those multiplicities for NodeReps and EdgeReps in nodeMult resp. edgeMult
         * If the dBG is disconnected, unnessecary dBG parts are skipped
         */
        void approxMultAll(NodeMap<Multiplicity>& nodeMult,
                           const CovModel& nodeCovModel,
                           EdgeMap<Multiplicity>& edgeMult,
                           const CovModel& edgeCovModel,
                           double& logZ,
                           bool MAP=false, bool singleCRF=false) const;
                           
        /**
         * TODO Remove!
         */
        void MstepWithTrueMult(const std::vector<NodeRep>& nodes,
                               CovModel& nodeCovModel,
                               const std::vector<EdgeRep>& edges,
                               CovModel& edgeCovModel, double epsilon);


        /**
         * Compute the multiplicities and models using expectation-maximization
         * @param nodeMult Estimated normalized node multiplicity (input/output)
         * @param nodeCovModel Node coverage model (input/output)
         * @param edgeMult Estimated normalized edge multiplicity (input/output)
         * @param edgeCovModel Edge coverage model (input/output)
         * @param epsilon Tolerance on the lambda_1 estimate
         * @param maxIter Maximum number of iterations
         */
        int computeMultEM(NodeMap<Multiplicity>& nodeMult,
                          CovModel& nodeCovModel,
                          EdgeMap<Multiplicity>& edgeMult,
                          CovModel& edgeCovModel,
                          double epsilon, double maxChange, int maxIter,
                          bool eqW = false,
                          bool approxInf = false, bool map = false,
                          bool fixedZero = false);

        /**
         * Get the CRF performance counter
         * @return The CRF performance counter
         */
        CRFPerfCounter getPerfCounter() const {
                return totPerfCounter;
        }
};

#endif
