#ifndef NETWORK_H
#define NETWORK_H

#include <vector>
#include <set>
#include <algorithm>

#include "umatrix.h"

class NetworkEdge
{
public:
    NetworkEdge();
    NetworkEdge(int source, int destination, int payload);
    bool operator==(const NetworkEdge& edge);
    bool operator!=(const NetworkEdge& edge);

    int source;
    int destination;
    /* Test payload */
    int payload;
    //UMatrix U;
};

class NetworkVertex
{
public:
    NetworkVertex();
    int index;
    std::vector<int> arrivingEdges;
    std::vector<int> leavingEdges;

};

class Network
{
public:
    Network(int vertexCount);
    void addEdge(int source, int destination, int payload);
    std::vector<NetworkEdge> getEdges();
    std::vector<NetworkEdge> getAdjacentEdges(NetworkEdge edge);
private:
    void addAdjacentEdges(int networkNode, std::vector<NetworkEdge>& adjacencyList);
    bool contains(std::vector<NetworkEdge>& adjacencyList, NetworkEdge edge);
    std::vector<NetworkVertex> vertexList;
    std::vector<NetworkEdge> edgeList;
};

#endif // NETWORK_H
