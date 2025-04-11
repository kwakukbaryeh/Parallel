#include "sequential.h"
#include <cmath>
#include <queue>
#include <limits>
#include <algorithm>
#include <iostream>
using namespace std;

// Heuristic function: Euclidean distance between nodes a and b.
float heuristic(const Node &a, const Node &b) {
    return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

// A* search
bool a_star(const Graph &graph, int start, int goal, vector<int> &path) {
    int n = graph.nodes.size();
    vector<float> gScore(n, INF);
    vector<float> fScore(n, INF);
    vector<int> cameFrom(n, -1);
    vector<bool> closed(n, false);
    
    // Priority queue for open set.
    priority_queue<AStarNode, vector<AStarNode>, greater<AStarNode>> openSet;
    
    gScore[start] = 0.0;
    fScore[start] = heuristic(graph.nodes[start], graph.nodes[goal]);
    openSet.push({start, gScore[start], fScore[start], -1});
    
    while (!openSet.empty()) {
        AStarNode current = openSet.top();
        openSet.pop();
        
        if (current.id == goal) {
            // Reconstruct the path.
            int cur = goal;
            path.clear();
            while (cur != -1) {
                path.push_back(cur);
                cur = cameFrom[cur];
            }
            reverse(path.begin(), path.end());
            return true;
        }
        
        if (closed[current.id]) continue;
        closed[current.id] = true;
        cameFrom[current.id] = current.parent;
        
        for (int edgeIdx : graph.adj[current.id]) {
            const Edge &edge = graph.edges[edgeIdx];
            int neighbor = (edge.start == current.id) ? edge.end : edge.start;
            if (closed[neighbor]) continue;
            
            float tentative_gScore = gScore[current.id] + edge.curr_weight;
            if (tentative_gScore < gScore[neighbor]) {
                gScore[neighbor] = tentative_gScore;
                fScore[neighbor] = tentative_gScore + heuristic(graph.nodes[neighbor], graph.nodes[goal]);
                cameFrom[neighbor] = current.id;
                openSet.push({neighbor, gScore[neighbor], fScore[neighbor], current.id});
            }
        }
    }
    return false; // No path found.
}

// Reset and update edge loads based on the given vehicle routes.
void update_edge_loads(Graph &graph, const vector<vector<int>> &vehicle_routes) {
    // Reset each edge's load.
    for (auto &edge : graph.edges) {
        edge.load = 0;
    }
    // Increment load for each edge used in vehicle routes.
    for (const auto &route : vehicle_routes) {
        for (size_t i = 0; i + 1 < route.size(); i++) {
            int u = route[i];
            int v = route[i+1];
            // Find the edge connecting u and v.
            for (int edgeIdx : graph.adj[u]) {
                const Edge &e = graph.edges[edgeIdx];
                if ((e.start == u && e.end == v) || (e.start == v && e.end == u)) {
                    graph.edges[edgeIdx].load++;
                    break;
                }
            }
        }
    }
}

// Update current weights for edges based on congestion (load exceeding capacity).
void update_edge_weights(Graph &graph) {
    for (auto &edge : graph.edges) {
        if (edge.load > edge.capacity) {
            // Example: increase weight proportionally to overload.
            edge.curr_weight = edge.base_weight * (1.0 + CONGESTION_FACTOR * (edge.load - edge.capacity));
        } else {
            edge.curr_weight = edge.base_weight;
        }
    }
}

// Graph with 4 nodes in a square and diagonals.
Graph create_test_graph() {
    Graph graph;
    graph.nodes = {
        {0, 0.0, 0.0},
        {1, 1.0, 0.0},
        {2, 0.0, 1.0},
        {3, 1.0, 1.0}
    };
    // Define edges: square perimeter and two diagonals.
    // Each edge has a capacity of 1.
    graph.edges = {
        {0, 1, 1.0, 1.0, 1, 0},    // edge 0: between node 0 and 1
        {1, 3, 1.0, 1.0, 1, 0},    // edge 1: between node 1 and 3
        {3, 2, 1.0, 1.0, 1, 0},    // edge 2: between node 3 and 2
        {2, 0, 1.0, 1.0, 1, 0},    // edge 3: between node 2 and 0
        {0, 3, 1.414, 1.414, 1, 0},// edge 4: diagonal between 0 and 3
        {1, 2, 1.414, 1.414, 1, 0} // edge 5: diagonal between 1 and 2
    };
    
    // Build the adjacency list.
    graph.adj.resize(graph.nodes.size());
    for (int i = 0; i < graph.edges.size(); i++) {
        int u = graph.edges[i].start;
        int v = graph.edges[i].end;
        graph.adj[u].push_back(i);
        graph.adj[v].push_back(i);
    }
    
    return graph;
}


Graph create_large_test_graph(int rows, int cols) {
    Graph graph;
    
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            Node node;
            node.id = r * cols + c;
            node.x = (float)c;
            node.y = (float)r;
            graph.nodes.push_back(node);
        }
    }
    
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            int current = r * cols + c;
            if (c < cols - 1) {
                int right = r * cols + (c + 1);
                Edge e;
                e.start = current; e.end = right;
                e.base_weight = 1.0;
                e.curr_weight = 1.0;
                e.capacity = 1;
                e.load = 0;
                graph.edges.push_back(e);
            }
            if (r < rows - 1) {
                int down = (r + 1) * cols + c;
                Edge e;
                e.start = current; e.end = down;
                e.base_weight = 1.0;
                e.curr_weight = 1.0;
                e.capacity = 1;
                e.load = 0;
                graph.edges.push_back(e);
            }
            if (r < rows - 1 && c < cols - 1) {
                int downRight = (r + 1) * cols + (c + 1);
                Edge e;
                e.start = current; e.end = downRight;
                e.base_weight = 1.414; 
                e.curr_weight = 1.414;
                e.capacity = 1;
                e.load = 0;
                graph.edges.push_back(e);
            }
            if (r < rows - 1 && c > 0) {
                int downLeft = (r + 1) * cols + (c - 1);
                Edge e;
                e.start = current; e.end = downLeft;
                e.base_weight = 1.414;
                e.curr_weight = 1.414;
                e.capacity = 1;
                e.load = 0;
                graph.edges.push_back(e);
            }
        }
    }
    
    graph.adj.resize(graph.nodes.size());
    for (int i = 0; i < graph.edges.size(); i++) {
        int u = graph.edges[i].start;
        int v = graph.edges[i].end;
        graph.adj[u].push_back(i);
        graph.adj[v].push_back(i);
    }
    
    return graph;
}