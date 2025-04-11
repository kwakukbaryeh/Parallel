#include "sequential.h"
#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
using namespace std;

int main() {
    srand(time(nullptr));
    
    int rows = 5, cols = 5;
    Graph graph = create_large_test_graph(rows, cols);
    
    int numVehicles = 10;
    vector<Vehicle> vehicles;
    
    for (int i = 0; i < numVehicles; i++) {
        int src = rand() % (rows * cols);
        int dest = rand() % (rows * cols);
        while (dest == src)
            dest = rand() % (rows * cols);
        vehicles.push_back({src, dest});
    }
    
    cout << "Vehicle queries:" << endl;
    for (int i = 0; i < vehicles.size(); i++) {
        cout << "Vehicle " << i << ": " 
             << vehicles[i].source << " -> " << vehicles[i].destination << endl;
    }
    cout << endl;
    
    vector<vector<int>> routes(numVehicles), prev_routes(numVehicles);
    
    int iteration = 0;
    bool converged = false;
    
    while (iteration < MAX_ITER && !converged) {
        cout << "Iteration " << iteration << ":\n";
        
        for (int i = 0; i < numVehicles; i++) {
            vector<int> path;
            if (a_star(graph, vehicles[i].source, vehicles[i].destination, path)) {
                routes[i] = path;
                cout << "Vehicle " << i << " route: ";
                for (int node : path)
                    cout << node << " ";
                cout << endl;
            } else {
                cout << "Vehicle " << i << ": No route found.\n";
            }
        }
        
        update_edge_loads(graph, routes);
        update_edge_weights(graph);
        
        converged = true;
        for (int i = 0; i < numVehicles; i++) {
            if (routes[i] != prev_routes[i]) {
                converged = false;
                break;
            }
        }
        
        if (converged) {
            cout << "Routes have converged.\n";
            break;
        }
        
        prev_routes = routes;
        iteration++;
        
        cout << "\nUpdated edge weights due to congestion (showing first few edges):\n";
        int displayEdges = min((int)graph.edges.size(), 10);
        for (int i = 0; i < displayEdges; i++) {
            cout << "Edge " << i 
                 << " (" << graph.edges[i].start << "-" << graph.edges[i].end << "): "
                 << "Load = " << graph.edges[i].load 
                 << ", New Weight = " << graph.edges[i].curr_weight << "\n";
        }
        cout << "\n";
    }
    
    cout << "Final Routes after " << iteration << " iterations:\n";
    for (int i = 0; i < numVehicles; i++) {
        cout << "Vehicle " << i << " final route: ";
        for (int node : routes[i])
            cout << node << " ";
        cout << endl;
    }
    
    return 0;
}
