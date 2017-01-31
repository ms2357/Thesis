#include "network.h"

NetworkEdge::NetworkEdge()
{

}

NetworkEdge::NetworkEdge(int source, int destination, int payload)
        : source(source), destination(destination), payload(payload)
{

}

NetworkVertex::NetworkVertex()
{

}

bool NetworkEdge::operator ==(const NetworkEdge& edge)
{
    return (this->source == edge.source) &
           (this->destination == edge.destination);
}

bool NetworkEdge::operator !=(const NetworkEdge& edge)
{
    return (this->source != edge.source) |
           (this->destination != edge.destination);
}

Network::Network(int vertexCount)
{
    vertexList.resize(vertexCount);
    for (int i = 0; i < vertexCount; i++) {
        vertexList[i].index = i;
    }
}

void Network::addEdge(int source, int destination, int payload)
{
    NetworkEdge newEdge(source, destination, payload);
    edgeList.push_back(newEdge);
    int newEdgeIndex = edgeList.size() - 1;
    vertexList[source].leavingEdges.push_back(newEdgeIndex);
    vertexList[destination].arrivingEdges.push_back(newEdgeIndex);
}

std::vector<NetworkEdge> Network::getEdges() { return edgeList; }

std::vector<NetworkEdge> Network::getAdjacentEdges(NetworkEdge edge)
{
    std::set<int> adjacentEdgeIndices;
    NetworkVertex &vertex = vertexList[edge.source];
    adjacentEdgeIndices.insert(vertex.arrivingEdges.begin(), vertex.arrivingEdges.end());
    adjacentEdgeIndices.insert(vertex.leavingEdges.begin(), vertex.leavingEdges.end());
    vertex = vertexList[edge.destination];
    adjacentEdgeIndices.insert(vertex.arrivingEdges.begin(), vertex.arrivingEdges.end());
    adjacentEdgeIndices.insert(vertex.leavingEdges.begin(), vertex.leavingEdges.end());

    std::vector<NetworkEdge> adjacentEdges;
    for (auto &edgeIndex : adjacentEdgeIndices) {
        NetworkEdge adjacentEdge = edgeList[edgeIndex];
        if (adjacentEdge != edge) {
            adjacentEdges.push_back(adjacentEdge);
        }
    }

    return adjacentEdges;
}

void Network::addAdjacentEdges(int networkNode, std::vector<NetworkEdge>& adjacencyList)
{
}

bool Network::contains(std::vector<NetworkEdge>& adjacencyList, NetworkEdge edge)
{
   return true;
}
