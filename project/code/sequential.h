#ifndef SEQUENTIAL_H
#define SEQUENTIAL_H

#include <vector>
using namespace std;

const float INF = 1e9;
const float CONGESTION_FACTOR = 1.0;   // Increase multiplier for overloaded edges
const int MAX_ITER = 10;               // Maximum number of re-routing iterations

// Structure for a node in the graph
struct Node {
    int id;
    float x, y;
};

// Structure for an edge in the graph
struct Edge {
    int start, end;         // Node IDs for endpoints (undirected)
    float base_weight;      // Nominal travel time/distance
    float curr_weight;      // Dynamically adjusted weight (based on congestion)
    int capacity;           // Maximum vehicles allowed concurrently
    int load;               // Current number of vehicles using the edge
};

// Graph structure containing nodes, edges, and the adjacency list.
struct Graph {
    vector<Node> nodes;
    vector<Edge> edges;
    // For each node, a vector of indices into "edges" representing incident edges.
    vector<vector<int>> adj;
};

// Vehicle query: specifies source and destination node IDs.
struct Vehicle {
    int source;
    int destination;
};

// Structure for nodes used in A* search
struct AStarNode {
    int id;
    float g;       // Cost from start
    float f;       // f = g + heuristic estimate
    int parent;    // Parent node id (for path reconstruction)
    // Operator for priority queue ordering (min-heap)
    bool operator>(const AStarNode &other) const {
        return f > other.f;
    }
};

// Function prototypes
// Heuristic: Euclidean distance between two nodes.
float heuristic(const Node &a, const Node &b);

// A* search: returns true if a path is found; the path is stored in 'path'.
bool a_star(const Graph &graph, int start, int goal, vector<int> &path);

// Update edge loads based on the current set of vehicle routes.
void update_edge_loads(Graph &graph, const vector<vector<int>> &vehicle_routes);

// Update the current weight on each edge based on its current load.
void update_edge_weights(Graph &graph);

// Generate a simple test graph.
Graph create_test_graph();

#endif // SEQUENTIAL_H
