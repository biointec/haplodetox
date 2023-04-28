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

#include <sstream>
#include <random>
#include <utility>

#include "coverage.h"
#include "dbgraph.h"
#include "settings.h"
#include "util.h"

DSNode* SSNode::nodes = NULL;

using namespace std;

std::ostream &operator<<(std::ostream &out, const NodeChain& nc)
{
        for (size_t i = 0; i < nc.size(); i++)
                out << nc[i] << "\t";
        out << "(" << nc.getCount() << ")";
        return out;
}

bool DBGraph::consecutiveNPP(const NodePosPair& left,
                             const NodePosPair& right,
                             size_t offset) const
{
        // return false if one of the npps is invalid
        if (!left.isValid() || !right.isValid())
                return false;

        // left and right belong to the same node?
        if (left.getNodeID() == right.getNodeID())
                if ((left.getPosition() + offset) == right.getPosition())
                        return true;

        // left and right belong to connected nodes?
        SSNode lNode = getSSNode(left.getNodeID());
        if (lNode.rightArc(right.getNodeID()) == NULL)
                return false;

        // make sure the distance is consistent
        size_t dist = lNode.getMarginalLength() - left.getPosition() + right.getPosition();
        return (offset == dist);
}

void DBGraph::writeBCalm (const string& filename)
{
        ofstream ofs(filename.c_str());
        if (!ofs.good())
                cerr << "Cannot open file " << filename
                     << " for writing" << endl;

        for (NodeID id = 1; id <= numNodes; id++)
        {
                SSNode node = getSSNode(id);
                if (! node.isValid())
                        continue;

                ofs << ">" << id << " LN:i:" <<
                node.getMarginalLength() + Kmer::getK() -1 << " KC:f:" <<
                node.getCov() << " km:f:" << node.getAvgCov();
                for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++) {
                        NodeID rightID = it->getNodeID();
                        ofs << "\tL:+:" << abs(rightID) << ":" << (rightID > 0 ? "+" : "-");
                }
                for (ArcIt it = node.leftBegin(); it != node.leftEnd(); it ++) {
                        NodeID leftID = it->getNodeID();
                        ofs << "\tL:-:" << abs(leftID) << ":" << (leftID > 0 ? "-" : "+");
                }
                ofs << "\n";
                ofs << node.getSequence() << endl;
        }
        ofs.close();
}

void DBGraph::loadBCalm (const string& filename)
{
        ifstream ifs(filename.c_str());
        if (!ifs)
                throw ios_base::failure("Cannot open file " + filename);

        // figure out file size
        ifs.seekg(0, ios_base::end);
        size_t fileSize = ifs.tellg();
        ifs.clear();
        ifs.seekg(0, ios::beg);

        // first pass through the file to find out number of nodes and arcs
        Util::startChrono();
        string progressStr = "Reading file " + filename + " (first pass)";

        numNodes = numArcs = 0;

        string line;
        while (getline(ifs, line)) {
                if (numNodes % 1024 == 0)
                        Util::progress(progressStr, ifs.tellg(), fileSize);
                if (line.front() != '>')
                        continue;

                numNodes++;

                // every occurrence of "L:" denotes a new arc
                size_t pos = 0;
                while ((pos = line.find("L:", pos + 1, 2)) != string::npos)
                        numArcs++;
        }

        double elapsed = Util::stopChrono();
        Util::progressEnd(progressStr, elapsed);

        // +1 because index 0 is not used
        nodes = new DSNode[numNodes+1];
        SSNode::setNodePointer(nodes);

        // +2 because index 0 is not used, final index denotes 'end'.
        arcs = new Arc[numArcs+2];
        DSNode::setArcsPointer(arcs);

        // second pass through the file to store nodes and arcs
        cout << "Reading " << numNodes << " nodes and "
             << numArcs << " arcs..." << endl;

        ifs.clear();
        ifs.seekg(0, ios::beg);

        Util::startChrono();
        progressStr = "Reading file " + filename + " (second pass)";

        ArcID arcOffset = 1; NodeID nodeOffset = 1;
        while (getline(ifs, line)) {
                if (nodeOffset % 1024 == 0)
                        Util::progress(progressStr, ifs.tellg(), fileSize);
                DSNode& node = getDSNode(nodeOffset);

                // the line is a sequence
                if (line.front() != '>') {
                        node.setSequence(line);
                        nodeOffset++;
                        continue;
                }

                // the line is a FASTA descriptor line
                istringstream iss(line);

                // find "KC:i:" to figure out the k-mer coverage
                size_t pos = 0; int kmerCov = 0;
                if ((pos = line.find("KC:i:")) != string::npos) {
                        iss.seekg(pos + 5);
                        iss >> kmerCov;

                }

                node.setCov(kmerCov);

                vector<NodeID> leftArcs;
                vector<NodeID> rightArcs;

                // every occurrence of "L:" denotes a new arc
                pos = 0;
                while ((pos = line.find("L:", pos + 1, 2)) != string::npos) {
                        iss.seekg(pos + 2);
                        char c, l, r;
                        int dstID;
                        iss >> l >> c >> dstID >> c >> r >> right;
                        dstID++;        // we number from 1, not 0

                        if (l == '+')
                                rightArcs.push_back(r == '+' ? dstID : -dstID);
                        else
                                leftArcs.push_back(r == '-' ? dstID : -dstID);
                }

                node.setNumLeftArcs(leftArcs.size());
                node.setFirstLeftArcID(arcOffset);
                for (auto dstID : leftArcs)
                        arcs[arcOffset++].setNodeID(dstID);

                node.setNumRightArcs(rightArcs.size());
                node.setFirstRightArcID(arcOffset);
                for (auto dstID : rightArcs)
                        arcs[arcOffset++].setNodeID(dstID);
        }

        elapsed = Util::stopChrono();
        Util::progressEnd(progressStr, elapsed);
        ifs.close();

        // figure out the value for k
        bool autoDetectK = false;

        for (NodeID id = -numNodes; id <= numNodes; id++) {
                if (id == 0)
                        continue;

                SSNode node = getSSNode(id);

                if (node.numRightArcs() < 2)
                        continue;

                ArcIt it = node.rightBegin();
                SSNode nb1 = getSSNode(it->getNodeID());
                it++;
                SSNode nb2 = getSSNode(it->getNodeID());
                string s1 = nb1.getSequence();
                string s2 = nb2.getSequence();

                int k = 0;
                while (s1[k] == s2[k])          // overlap == k-1
                        k++;
                k++;                            // add 1 to get k

                Kmer::setWordSize(k);
                autoDetectK = true;
                break;
        }

        if (!autoDetectK)
                throw runtime_error("Cannot infer kmer size from input file");
        else
                cout << "Kmer size is " << Kmer::getK() << endl;

        numValidNodes = numNodes;
        numValidArcs = numArcs;
}

void DBGraph::writeBinary(const std::string& filename) const
{
        ofstream ofs(filename.c_str(), ios::binary);
        // write k
        size_t k = Kmer::getK();
        ofs.write((char*)&k, sizeof(k));

        // A) write nodes
        ofs.write((char*)&numNodes, sizeof(numNodes));
        for (NodeID id = 1; id <= numNodes; id++)
                getDSNode(id).write(ofs);

        // B) write arcs
        ofs.write((char*)&numArcs, sizeof(numArcs));
        for (ArcID i = 0; i < numArcs+2; i++)
                arcs[i].write(ofs);

        ofs.close();

        cout << "Wrote " << numNodes << " nodes and "
             << numArcs << " arcs" << endl;
}

void DBGraph::loadBinary(const std::string& filename)
{
        // read the metadata
        ifstream ifs(filename.c_str(), ios::binary);
        if (!ifs)
                throw ios_base::failure("Cannot open " + filename);

        size_t k;
        ifs.read((char*)&k, sizeof(k));
        Kmer::setWordSize(k);

        // A) create the nodes
        ifs.read((char*)&numNodes, sizeof(numNodes));
        // +1 because index 0 isn't used
        nodes = new DSNode[numNodes+1];
        SSNode::setNodePointer(nodes);
        for (NodeID id = 1; id <= numNodes; id++)
                getDSNode(id).read(ifs);

        // B) create the arcs
        ifs.read((char*)&numArcs, sizeof(numArcs));
        // +2 because index 0 isn't used, final index denotes 'end'
        arcs = new Arc[numArcs+2];
        DSNode::setArcsPointer(arcs);
        for (ArcID i = 0; i < numArcs+2; i++)
                arcs[i].read(ifs);

        ifs.close();

        size_t countTotal = 0;
        for (NodeID seedID = 1; seedID <= numNodes; seedID++)
                if (getSSNode(seedID).isValid())
                        countTotal++;

        numValidNodes = countTotal;
        countTotal = 0;
        for (ArcID seedID = 1; seedID <= numArcs; seedID++)
                if (getArc(seedID).isValid())
                        countTotal++;
        numValidArcs = countTotal;
}

void DBGraph::writeContigs(const std::string& filename) const
{
        ofstream ofs(filename.c_str());

        size_t contigID = 0;
        for (NodeID id = 1; id <= getNumNodes(); id++) {
                SSNode n = getSSNode(id);
                if (!n.isValid())
                        continue;

                //ofs << ">contig_" << contigID++ << "\n";
                ofs << ">contig_" << id << "\n";
                Util::writeSeqWrap(ofs, n.getSequence(), 60);
        }
}

void DBGraph::getGraph(std::vector<NodeID>& nodes, std::vector<EdgeID>& edges)
{
        nodes.reserve(2 * getNumValidNodes());
        edges.reserve(getNumValidArcs());

        for (NodeID srcID = -numNodes; srcID <= numNodes; srcID++) {
                if (srcID == 0)
                        continue;
                SSNode n = getSSNode(srcID);
                if (!getSSNode(srcID).isValid())
                        continue;
                nodes.push_back(srcID);

                for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++) {
                        NodeID dstID = it->getNodeID();
                        edges.push_back(make_pair(srcID, dstID));
                }
        }
}

void DBGraph::getSubgraph(NodeID seedNode, vector<NodeID>& nodes,
                          vector<EdgeID>& edges, size_t maxDepth) const
{
        priority_queue<NodeDFS, vector<NodeDFS>, NodeDFSComp> todo;
        todo.push(NodeDFS(seedNode, 0));
        set<NodeID> nodeSet;

        while (!todo.empty()) {
                // get and erase the current node
                NodeDFS currTop = todo.top();
                todo.pop();
                NodeID thisID = currTop.nodeID;
                size_t thisDepth = currTop.depth;

                SSNode n = getSSNode(thisID);

                // if the node was already handled, skip
                if (nodeSet.find(thisID) != nodeSet.end())
                        continue;


                // don't go any deeper if we've reached the maximum depth
                if (thisDepth >= maxDepth)
                        continue;

                // process the right arcs
                for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++) {
                        NodeID rightID = it->getNodeID();
                        if (nodeSet.find(rightID) != nodeSet.end())
                                continue;       // edge already added by rightID

                        edges.push_back(make_pair(thisID, rightID));
                        todo.push(NodeDFS(rightID, thisDepth+1));
                }

                // process the left arcs
                for (ArcIt it = n.leftBegin(); it != n.leftEnd(); it++) {
                        NodeID leftID = it->getNodeID();

                        if (nodeSet.find(leftID) != nodeSet.end())
                                continue;       // edge already added by leftID
                        if (leftID == thisID)   // edge already added as right arc
                                continue;
                                

                        edges.push_back(make_pair(leftID, thisID));
                        todo.push(NodeDFS(leftID, thisDepth + 1));
                }
                
                nodeSet.insert(thisID);
                
        }
        
        nodes = vector<NodeID>(nodeSet.begin(), nodeSet.end());
}

void DBGraph::getFullDirGraph(vector<NodeID>& nodes,
                              vector<pair<NodeID, NodeID> >& edges) const
{
        for (NodeID id = 1; id <= getNumNodes(); id ++){
                DSNode& node = getDSNode(id);
                if (!node.isValid())
                        continue;

                if (node.getFlag1())
                        continue;
                priority_queue<NodeDFS, vector<NodeDFS>, NodeDFSComp> todo;
                todo.push(NodeDFS(id, 0));

                while (!todo.empty()) {
                        // get and erase the current node
                        NodeDFS currTop = todo.top();
                        todo.pop();
                        NodeID thisID = currTop.nodeID;
                        size_t thisDepth = currTop.depth;

                        SSNode n = getSSNode(thisID);

                        // if the node was already handled, skip
                        if (n.getFlag1())
                                continue;

                        nodes.push_back(thisID);

                        // process the right arcs
                        for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++) {
                                NodeID rightID = it->getNodeID();
                                SSNode r = getSSNode(rightID);
                                if (r.getFlag1())       // edge already added by r
                                        continue;

                                edges.push_back(make_pair(thisID, rightID));
                                todo.push(NodeDFS(rightID, thisDepth+1));
                        }

                        // process the left arcs
                        for (ArcIt it = n.leftBegin(); it != n.leftEnd(); it++) {
                                NodeID leftID = it->getNodeID();
                                SSNode l = getSSNode(leftID);
                                if (l.getFlag1())       // edge already added by l
                                        continue;

                                edges.push_back(make_pair(leftID, thisID));
                                todo.push(NodeDFS(leftID, thisDepth + 1));
                        }

                        // mark this node as handled
                        n.setFlag1(true);
                }
        }

        // reset all flags to false
        for (auto it : nodes)
                getSSNode(it).setFlag1(false);
}

/*void DBGraph::getSubgraph(priority_queue<NodeRepDepth, vector<NodeRepDepth>,
                          NodeRepComp>& todo, set<NodeRep>& nodes,
                          set<EdgeRep>& edges, size_t maxDepth) const
{
        while (!todo.empty()) {
                // get and erase the current node
                NodeRepDepth currTop = todo.top();
                todo.pop();
                NodeRep thisID = currTop.nodeRep;
                size_t thisDepth = currTop.depth;

                // if the node was already handled, skip
                if (nodes.find(thisID) != nodes.end())
                        continue;

                // mark this node as handled
                nodes.insert(thisID);

                // process the right arcs
                SSNode node = getSSNode(thisID);
                for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++) {
                        NodeID rightID = it->getNodeID();
                        edges.insert(EdgeRep(thisID, rightID));
                        if (thisDepth < maxDepth)
                                todo.push(NodeRepDepth(rightID, thisDepth+1));
                }

                // process the left arcs
                for (ArcIt it = node.leftBegin(); it != node.leftEnd(); it++) {
                        NodeID leftID = it->getNodeID();
                        edges.insert(EdgeRep(leftID, thisID));
                        if (thisDepth < maxDepth)
                                todo.push(NodeRepDepth(leftID, thisDepth + 1));
                }
        }
}

void DBGraph::getSubgraph(NodeRep seedNode, set<NodeRep>& nodes,
                          set<EdgeRep>& edges, size_t maxDepth) const
{
        // a queue containing nodeIDs to handle + their depth
        priority_queue<NodeRepDepth, vector<NodeRepDepth>, NodeRepComp> todo;
        todo.push(NodeRepDepth(seedNode, 0));
        
        getSubgraph(todo, nodes, edges, maxDepth);
}
                          
void DBGraph::getSubgraph(EdgeRep seedEdge, set<NodeRep>& nodes,
                        set<EdgeRep>& edges, size_t maxDepth) const
{
        edges.insert(seedEdge);
        if (maxDepth == 0)
                return;
        
        // a queue containing nodeIDs to handle + their depth
        priority_queue<NodeRepDepth, vector<NodeRepDepth>, NodeRepComp> todo;
        todo.push(NodeRepDepth(seedEdge.getSrcID(), 1));
        todo.push(NodeRepDepth(seedEdge.getDstID(), 1));
        
        getSubgraph(todo, nodes, edges, maxDepth);
}*/


vector<NodeRep> DBGraph::getNodeReps(const string& filename) const
{
        ifstream ifs(filename.c_str());

        NodeID n;
        size_t numNodes = 0;
        while (ifs >> n)
                numNodes++;

        vector<NodeRep> nodes;
        nodes.reserve(numNodes);

        ifs.clear();
        ifs.seekg(0, ios::beg);

        while (ifs >> n) {
                if (!nodeExists(n))
                        throw runtime_error("Node with ID " + to_string(n) + " does not exist");
                nodes.push_back(NodeRep(n));
        }

        return nodes;
}

vector<NodeRep> DBGraph::getNodeReps(size_t N) const
{
        vector<NodeRep> nodes;
        nodes.reserve(getNumValidNodes());
        for (size_t i = 1; i <= numNodes; i++)
                if (getSSNode(i).isValid())
                        nodes.push_back(NodeRep(i));

        N = min(N, nodes.size());
        if (N == nodes.size())  // if you need all nodes we do not shuffle
                return nodes;

        // sample N nodes using the Fisher-Yates algoritm
        random_device rd;
        mt19937 mt(rd());
        for (size_t i = 0; i < N; i++) {
                uniform_int_distribution<size_t> dis(i, nodes.size() - 1);
                swap(nodes[i], nodes[dis(mt)]);
        }

        return vector<NodeRep>(nodes.begin(), nodes.begin() + N);
}

vector<NodeRep> DBGraph::getEvalNodes(double min, double max, size_t N, const string& filename) const
{
        vector<NodeRep> nodes;
        nodes.reserve(getNumValidNodes());
        for (size_t i = 1; i <= numNodes; i++){
                double cov = getSSNode(i).getAvgCov();
                if (getSSNode(i).isValid() && cov > min && cov < max)
                        if(getSSNode(i).numLeftArcs() > 0 && getSSNode(i).numRightArcs() > 0)
                                nodes.push_back(NodeRep(i));
        }
                
        N = std::min(N, nodes.size());
        if (N == nodes.size())  // if you need all nodes we do not shuffle
                return nodes;
        
        // sample N nodes using the Fisher-Yates algoritm
        random_device rd;
        mt19937 mt(rd());
        for (size_t i = 0; i < N; i++) {
                uniform_int_distribution<size_t> dis(i, nodes.size() - 1);
                swap(nodes[i], nodes[dis(mt)]);
        }
        
        //ofstream fs(filename.c_str());
        
        //for(int i = 0; i < N; i++)
        //        fs << nodes[i] << "\n";
        
        return vector<NodeRep>(nodes.begin(), nodes.begin() + N);
}

vector<NodeRep> DBGraph::getEvalNodes(const string& filename) const
{
        ifstream ifs(filename.c_str());
        
        NodeID n;
        size_t numNodes = 0;
        while (ifs >> n)
                numNodes++;
        
        vector<NodeRep> nodes;
        nodes.reserve(numNodes);
        
        ifs.clear();
        ifs.seekg(0, ios::beg);
        
        while (ifs >> n) {
                if (!nodeExists(n))
                        throw runtime_error("Node with ID " + to_string(n) + " does not exist");
                nodes.push_back(NodeRep(n));
        }
        
        return nodes;
}


vector<EdgeRep> DBGraph::getEdgeReps(const string& filename) const
{
        ifstream ifs(filename.c_str());

        ArcID l, r;
        size_t numEdges = 0;
        while (ifs >> l >> r)
                numEdges++;

        vector<EdgeRep> edges;
        edges.reserve(numEdges);

        ifs.clear();
        ifs.seekg(0, ios::beg);

        while (ifs >> l >> r) {
                if (!edgeExists(EdgeRep(l, r)))
                        throw runtime_error("Edge from ID " + to_string(l) + " to ID " + to_string(r) + " does not exist");
                edges.push_back(EdgeRep(l, r));
        }

        return edges;
}

vector<EdgeRep> DBGraph::getEdgeReps(size_t N) const
{
        vector<EdgeRep> edges;
        edges.reserve(getNumValidArcs() / 2);   // roughly half of the nodes
        for (NodeID id = -numNodes; id <= numNodes; id++) {
                if (id == 0)
                        continue;
                SSNode n = getSSNode(id);
                if (!n.isValid())
                        continue;
                for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++)
                        if (id <= -it->getNodeID())
                                edges.push_back(EdgeRep(id, it->getNodeID()));
        }

        N = min(N, edges.size());
        if (N == edges.size())  // if you need all edges we do not shuffle
                return edges;

        // sample N edges using the Fisher-Yates algoritm
        random_device rd;
        mt19937 mt(rd());
        for (size_t i = 0; i < N; i++) {
                uniform_int_distribution<size_t> dis(i, edges.size() - 1);
                swap(edges[i], edges[dis(mt)]);
        }

        return vector<EdgeRep>(edges.begin(), edges.begin() + N);
}

vector<EdgeRep> DBGraph::getLowCovEdges(double threshold) const
{
        vector<EdgeRep> edges;
        edges.reserve(getNumValidArcs() / 2);   // roughly half of the nodes
        for (NodeID id = -numNodes; id <= numNodes; id++) {
                if (id == 0)
                        continue;
                SSNode n = getSSNode(id);
                if (!n.isValid())
                        continue;
                for (ArcIt it = n.rightBegin(); it != n.rightEnd(); it++) {
                        if (id <= -it->getNodeID() && it->getCov() <= threshold)
                                edges.push_back(EdgeRep(id, it->getNodeID()));
                }
        }

        return edges;
}

vector<NodeRep> DBGraph::getLowCovNodes(double threshold) const
{
        vector<NodeRep> nodes;
        nodes.reserve(getNumValidNodes());
        for (size_t i = 1; i <= numNodes; i++) {
                SSNode n = getSSNode(i);
                if (!n.isValid())
                        continue;

                if (n.getAvgCov() <= threshold)
                        nodes.push_back(NodeRep(i));
        }

        return nodes;
}
vector<NodeRep> DBGraph::getLowCovSmallNodes(double threshold, int maxLength) const
{
        vector<NodeRep> nodes;
        nodes.reserve(getNumValidNodes());
        for (size_t i = 1; i <= numNodes; i++) {
                SSNode n = getSSNode(i);
                if (!n.isValid())
                        continue;

                if(n.getMarginalLength() > maxLength)
                        continue;

                if (n.getAvgCov() <= threshold)
                        nodes.push_back(NodeRep(i));
        }

        return nodes;
}


vector<NodeRep> DBGraph::getLowCovTips(double threshold, size_t maxLen) const
{
        vector<NodeRep> nodes;
        for (size_t i = 1; i <= numNodes; i++) {
                SSNode n = getSSNode(i);
                if (!n.isValid())
                        continue;       // deleted node

                if ((n.numLeftArcs() > 0) && (n.numRightArcs() > 0))
                        continue;       // not a tip

                if ((maxLen != 0) && (n.getMarginalLength() > maxLen))
                        continue;       // node too long

                if (n.getAvgCov() <= threshold)
                        nodes.push_back(NodeRep(i));
        }

        return nodes;
}

vector<NodeRep> DBGraph::getLowCovSmallTips(double threshold, int maxLength) const
{
        vector<NodeRep> nodes;
        nodes.reserve(getNumValidNodes());
        for (size_t i = 1; i <= numNodes; i++) {
                SSNode n = getSSNode(i);
                if (!n.isValid())
                        continue;       // deleted node

                if ((n.numLeftArcs() > 0) && (n.numRightArcs() > 0))
                        continue;       // not a tip

                if (n.getMarginalLength() > maxLength)
                        continue;

                if (n.getAvgCov() <= threshold)
                        nodes.push_back(NodeRep(i));
        }

        return nodes;
}

vector<NodeRep> DBGraph::getLowCovBubbles(double threshold) const
{
        vector<NodeRep> nodes;
        for (size_t i = 1; i <= numNodes; i++) {
                SSNode n = getSSNode(i);
                if (!n.isValid())
                        continue;       // deleted node

                if ((n.numLeftArcs() != 1) || (n.numRightArcs() != 1))
                        continue;       // not a bubble

                /*NodeID leftID = n.leftBegin()->getNodeID();
                NodeID rightID = n.rightBegin()->getNodeID();

                bool isBubble = false;  // find parallel path
                SSNode l = getSSNode(leftID);
                for (ArcIt it = l.rightBegin(); it != l.rightEnd(); it++) {
                        if (it->getNodeID() == i)
                                continue;
                        if (edgeExists(EdgeRep(it->getNodeID(), rightID)))
                                isBubble = true;
                }

                if (!isBubble)
                        continue;*/

                if (n.getMarginalLength() != Kmer::getK())
                        continue;       // not a bubble

                if (n.getAvgCov() <= threshold)
                        nodes.push_back(NodeRep(i));
        }

        return nodes;
}

void DBGraph::getNodeCovHist(map<unsigned int, double>& hist) const
{
        for (NodeID id = 1; id <= numNodes; id++) {
                SSNode node = getSSNode(id);
                if (!node.isValid())
                        continue;

                double f = node.getAvgCov() - floor(node.getAvgCov());
                double ML = node.getMarginalLength();

                hist[node.getAvgCov()] += (1.0-f) * ML;
                hist[node.getAvgCov() + 1] += f * ML;
        }
}

void DBGraph::writeCytoscapeGraph(const std::string& filename,
                                  vector<NodeID> nodes,
                                  vector<pair<NodeID, NodeID> > edges,
                                  const NodeMap<Multiplicity>& enm,
                                  const EdgeMap<Multiplicity>& eem,
                                  const NodeMap<vector<int>>& tnm,
                                  const EdgeMap<vector<int>>& tem) const
{
        // A) write all arcs
        ofstream ofs((filename + ".arcs").c_str());
        if (!ofs.good())
                cerr << "Cannot open file " << filename + ".arcs"
                     << " for writing" << endl;
     
        int nStrains = settings.getInitNStrains();
        ofs << "Source node\tTarget node\tCoverage";
        for (int i = 1; i <= nStrains; i++)
                ofs << "\tTrue mult. strain " << i;
        ofs << "\tTrue colour";
        for (int i = 1; i <= nStrains; i++)
                ofs << "\tEst mult. strain " << i;
        ofs << "\tEst. colour\tLOR\tBelief\tcorrect\n";

        for (const auto& edge : edges) {
                EdgeRep er(edge);
                ArcID id = getArcID(er);

                auto teIt = tem.find(er);
                vector<int> trueMult = (teIt == tem.end()) ? vector<int>(nStrains,-1) : teIt->second;

                int trueColour = 0;
                for (int i = 0; i < nStrains; i++)
                    trueColour += exp2(i) * ((trueMult[i] > 0) ? 1 : 0);

                auto eeIt = eem.find(er);
                vector<int> estMult = (eeIt == eem.end()) ? vector<int>(nStrains,-1) : eeIt->second.getExpMult();

                int estColour = 0;
                for (int i = 0; i < nStrains; i++)
                    estColour += exp2(i) * ((estMult[i] > 0) ? 1 : 0);

                double multLOR = (eeIt == eem.end()) ?
                        -1.0 : eeIt->second.getExpMultLogOR();
                double multBelief = (eeIt == eem.end()) ?
                        -1.0 : exp(eeIt->second.getExpMultLProb());

                ofs << edge.first << "\t" << edge.second << "\t"
                    << getArc(id).getCov();
                for (int i = 0; i < nStrains; i++)
                    ofs << "\t" << trueMult[i];
                ofs << "\t" << trueColour;
                for (int i = 0; i < nStrains; i++)
                    ofs << "\t" << estMult[i];
                ofs << "\t" << estColour << "\t" << multLOR << "\t" << multBelief
                    << "\t" << (trueMult==estMult? 1:0)<< "\n";
        }
        ofs.close();

        // B) write all nodes
        ofs.open((filename + ".nodes").c_str());
        if (!ofs.good())
                cerr << "Cannot open file " << filename + ".nodes"
                     << " for writing" << endl;

        ofs << "Node ID\tMarginal length\tLeft arcs\tRight arcs\tCoverage";
        for (int i = 1; i <= nStrains; i++)
                ofs << "\tTrue mult. strain " << i;
        ofs << "\tTrue colour";
        for (int i = 1; i <= nStrains; i++)
                ofs << "\tEst mult. strain " << i;
        ofs << "\tEst. colour\tLOR\tBelief\tcorrect\tSequence\n";

        for (const auto& id : nodes) {
                NodeRep nr(id);
                SSNode n = getSSNode(id);

                auto tnIt = tnm.find(nr);
                vector<int> trueMult = (tnIt == tnm.end()) ? vector<int>(nStrains,-1) : tnIt->second;

                int trueColour = 0;
                for (int i = 0; i < nStrains; i++)
                    trueColour += exp2(i) * ((trueMult[i] > 0) ? 1 : 0);

                auto enIt = enm.find(nr);
                vector<int> estMult = (enIt == enm.end()) ? vector<int>(nStrains,-1) : enIt->second.getExpMult();

                int estColour = 0;
                for (int i = 0; i < nStrains; i++)
                    estColour += exp2(i) * ((estMult[i] > 0) ? 1 : 0);

                double multLOR = (enIt == enm.end()) ?
                        -1.0 : enIt->second.getExpMultLogOR();
                double multBelief = (enIt == enm.end()) ?
                        -1.0 : exp(enIt->second.getExpMultLProb());

                ofs << id << "\t" << n.getMarginalLength() << "\t"
                    << (int)n.numLeftArcs() << "\t" << (int)n.numRightArcs()
                    << "\t" << n.getAvgCov();
                for (int i = 0; i < nStrains; i++)
                        ofs << "\t" << trueMult[i];
                ofs << "\t" << trueColour;
                for (int i = 0; i < nStrains; i++)
                        ofs << "\t" << estMult[i];
                ofs << "\t" << estColour << "\t" << multLOR << "\t" << multBelief<< "\t"
                    << (trueMult==estMult? 1:0) << "\t" << n.getSequence() << "\n";
        }
        ofs.close();
}

void DBGraph::writeStage1CytoscapeGraph(const std::string& filename,
                                  vector<NodeID> nodes,
                                  vector<pair<NodeID, NodeID>> edges,
                                  const NodeMap<std::vector<int>>& tnm,
                                  const EdgeMap<std::vector<int>>& tem) const
{
        size_t numstrains = tnm.begin()->second.size(); // assuming all nodes and edges will contain the same number of strains
        // A) write all arcs
        ofstream ofs((filename + ".arcs").c_str());
        if (!ofs.good())
                cerr << "Cannot open file " << filename + ".arcs"
                     << " for writing" << endl;

        ofs << "Source node\tTarget node\tCoverage";
        for (int i = 1; i <= numstrains; i++)
            ofs << "\tTrue mult. strain " << i;
        ofs << "\tColour\n";

        for (const auto& edge : edges) {
                EdgeRep er(edge);
                ArcID id = getArcID(er);

                auto teIt = tem.find(er);
                vector<int> trueStrainMult = (teIt == tem.end()) ? vector<int>(1,-1) : teIt->second;

                int colour = 0;
                for (int i = 0; i < numstrains; i++)
                    colour += exp2(i) * ((trueStrainMult[i] > 0) ? 1 : 0);

                // We assume a directed graph as input
                // (i.e. each bi-directional node occurs only once in node list,
                // either as its positive or as its negative nodeID)
                // However, the edge list can contain both the positive and the negative nodeID,
                // (e.g. in case of palindromes)
                // so take abs value to get a good cytoscape visualisation
                ofs << abs(edge.first) << "\t" << abs(edge.second) << "\t"
                    << getArc(id).getCov();
                for (int mult: trueStrainMult)
                        ofs << "\t" << mult;
                ofs << "\t" << colour << "\n";
        }
        ofs.close();

        // B) write all nodes
        ofs.open((filename + ".nodes").c_str());
        if (!ofs.good())
                cerr << "Cannot open file " << filename + ".nodes"
                     << " for writing" << endl;

        ofs << "Node ID\tMarginal length\tLeft arcs\tRight arcs\tCoverage";
        for (int i = 1; i <= numstrains; i++)
            ofs << "\tTrue mult. strain " << i;
        ofs << "\tColour\tSequence\n";

        for (const auto& id : nodes) {
                NodeRep nr(id);
                SSNode n = getSSNode(id);

                auto tnIt = tnm.find(nr);
                vector<int> trueStrainMult = (tnIt == tnm.end()) ? vector<int>(1,-1) : tnIt->second;

                int colour = 0;
                for (int i = 0; i < numstrains; i++)
                    colour += exp2(i) * ((trueStrainMult[i] > 0) ? 1 : 0);

                ofs << abs(id) << "\t" << n.getMarginalLength() << "\t"
                    << (int)n.numLeftArcs() << "\t" << (int)n.numRightArcs()
                    << "\t" << n.getAvgCov();
                for (int mult: trueStrainMult)
                        ofs << "\t" << mult;
                    ofs << "\t" << colour << "\t" << n.getSequence() << "\n";
        }
        ofs.close();
}

void DBGraph::writeStage1GFA(const std::string& filename,
                             vector<NodeID> nodes,
                             vector<pair<NodeID, NodeID>> edges,
                             const NodeMap<vector<int>>& tnm,
                             const EdgeMap<vector<int>>& tem) const
{
        size_t numstrains = tnm.begin()->second.size(); // assuming all nodes and edges will contain the same number of strains
        vector<string> colours = {"crimson", "goldenrod", "steelblue", "mediumseagreen", "hotpink", "lightsalmon", "darkviolet", "sienna", "linen", "khaki","lightsteelblue","lightgreen","lightpink", "peachpuff","violet","peru"};
        
        ofstream ofs((filename + ".gfa").c_str());
        
        if (!ofs.good())
                cerr << "Cannot open file " << filename + ".arcs"
                << " for writing" << endl;
        
        // B) write all nodes
        
        for (NodeID id = 1; id <= numNodes; id++) {
                NodeRep nr(id);
                SSNode n = getSSNode(id);
                if (! n.isValid())
                        continue;
                
                auto tnIt = tnm.find(nr);
                vector<int> trueStrainMult = (tnIt == tnm.end()) ? vector<int>(1,-1) : tnIt->second;
                
                int colour = 0;
                for (int i = 0; i < numstrains; i++)
                        colour += exp2(i) * ((trueStrainMult[i] > 0) ? 1 : 0);
                
                ofs << "S\t" << abs(id) <<"\t" << n.getSequence() << "\tLN:i:" << n.getMarginalLength() << 
                        "\tKC:i:" << n.getCov() << "\tCL:z:" << ((colour < colours.size()) ? colours[colour] : "darkgrey") <<  "\n";
        }
        
        // A) write all arcs
        for (NodeID id: nodes)
        {
                SSNode node = getSSNode(id);
                if (! node.isValid())
                        continue;
                
                for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++) {
                        NodeID rightID = it->getNodeID();
                        ofs << "L\t" << abs(id) << "\t" << (id > 0 ? "+" : "-") << "\t" << abs(rightID) << "\t" << (rightID > 0 ? "+" : "-") << "\t" << (Kmer::getK() - 1) <<  "M" << "\n";
                }
                /*for (ArcIt it = node.leftBegin(); it != node.leftEnd(); it ++) {
                        NodeID leftID = it->getNodeID();
                        ofs << "L\t" << abs(leftID) << "\t" << (leftID > 0 ? "-" : "+") << "\t" << id << "\t+" << "\t" << (Kmer::getK() - 1) << "M" << "\n";
                }*/
        }
        
        
        ofs.close();
}

void DBGraph::readTrueStrainMultiplicities(const Settings& settings, int nStrains, bool correct)
{
        trueNodeMult.clear();
        trueEdgeMult.clear();
        // If they exist, also load the true multiplicities from disk
        // Note: for efficiency, store only true multiplicities that are needed

//         transform(nodeReps.begin(), nodeReps.end(),
//                   inserter(trueNodeMult, trueNodeMult.end()),
//                   [](const NodeRep& nr) { return make_pair(nr, Multiplicity(-1)); });
//

//         transform(edgeReps.begin(), edgeReps.end(),
//                   inserter(trueEdgeMult, trueEdgeMult.end()),
//                   [](const EdgeRep& er) { return make_pair(er, Multiplicity(-1)); });
//
        ifstream ifs(settings.getTrueNodeStrainFilename(correct));
        if (!ifs)
                cerr << "Could not find true nodes multiplicities file...\n"
                        "True node multiplicities will be set to -1 in "
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
                trueNodeMult[nodeID] = strainM;
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
                for (int i=0; i < nStrains; i++){
                    ifs >> m;
                    strainM.push_back(m);
                }
                if (!ifs)
                        break;
                trueEdgeMult[EdgeRep(srcID, dstID)] = strainM;
        }
        ifs.close();
}


void DBGraph::defragNodes()
{
        vector<NodeID> old2new(numNodes + 1, 0);

        // defragment the node array
        for (NodeID oldID = 1, newID = 1; oldID <= numNodes; oldID++) {
                SSNode n = getSSNode(oldID);
                if (!n.isValid())
                        continue;

                old2new[oldID] = newID;

                // we use the move assignment operator for efficiency
                // (to avoid a deep copy of the TString member)
                nodes[newID] = move(nodes[oldID]);
                newID++;
        }

        // update the arcs to point to the new nodeIDs
        for (ArcID id = 1; id <= numArcs; id++) {
                NodeID oldID = arcs[id].getNodeID();
                NodeID newID = (oldID < 0) ? -old2new[-oldID] : old2new[oldID];
                arcs[id].setNodeID(newID);
        }

        numNodes = numValidNodes;
        cout << "\tDefragged nodes array: " << numNodes << " nodes\n";
}

void DBGraph::defragArcs()
{
        vector<ArcID> old2new(numArcs + 2, 0);

        numValidArcs = 0;
        for (ArcID oldID = 1; oldID <= numArcs + 1; oldID++) {
                if (!arcs[oldID].isValid())
                        continue;

                numValidArcs++;
                old2new[oldID] = numValidArcs;

                // we use the regular assignment operator
                arcs[numValidArcs] = arcs[oldID];
        }

        // update the nodes to point to the new arcIDs
        for (NodeID id = 1; id <= numNodes; id++) {
                SSNode n = getSSNode(id);
                if (!n.isValid())
                        continue;

                ArcID oldID = n.getFirstLeftArcID();
                ArcID newID = (oldID < 0) ? -old2new[-oldID] : old2new[oldID];
                n.setFirstLeftArcID(newID);

                oldID = n.getFirstRightArcID();
                newID = (oldID < 0) ? -old2new[-oldID] : old2new[oldID];
                n.setFirstRightArcID(newID);
        }

        numArcs = numValidArcs;
        cout << "\tDefragged arcs array: " << numArcs << " arcs\n";
}

void DBGraph::createKmerNPPTable(KmerNPPTable& table) const
{
        Util::startChrono();
        string progressStr("Populating <k-mer, node> table");

        // count the number of k-mers in the graph
        size_t numKmers = 0;
        for (NodeID id = 1; id <= numNodes; id++) {
                const DSNode& node = nodes[id];
                if (!node.isValid())
                        continue;
                numKmers += node.getMarginalLength();
        }

        table.resize(numKmers);

        // populate the table with kmers
        for (NodeID id = 1; id <= numNodes; id++) {
                if (id % 1024 == 0)
                        Util::progress(progressStr, id, numNodes);

                const DSNode &node = nodes[id];
                if (!node.isValid())
                        continue;
                const TString& tStr = node.getTSequence();
                Kmer kmer(tStr);
                NodePosition pos = 0;
                table.insert(kmer, NodePosPair(id, pos++));

                for (size_t i = Kmer::getK(); i < tStr.getLength(); i++) {
                        kmer.pushNucleotideRight(tStr[i]);
                        table.insert(kmer, NodePosPair(id, pos++));
                }
        }

        double elapsed = Util::stopChrono();
        Util::progressEnd(progressStr, elapsed);
}

int DBGraph::getAbundanceMin() const {
        return settings.getAbundanceMin();
}
