#include "sequential.h"
#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

int main() {
    // Create a test graph.
    Graph graph = create_test_graph();
    
    // Define test vehicle queries. For example, two vehicles with distinct routes.
    vector<Vehicle> vehicles = { {0, 3}, {1, 2} };
    int numVehicles = vehicles.size();
    
    // Allocate a vector to store routes for each vehicle.
    vector<vector<int>> routes(numVehicles);
    vector<vector<int>> prev_routes(numVehicles);  // For checking convergence.
    
    int iteration = 0;
    bool converged = false;
    
    while (iteration < MAX_ITER && !converged) {
        cout << "Iteration " << iteration << ":\n";
        // Compute routes for each vehicle using A*.
        for (int i = 0; i < numVehicles; i++) {
            vector<int> path;
            if (a_star(graph, vehicles[i].source, vehicles[i].destination, path)) {
                routes[i] = path;
                cout << "Vehicle " << i << " route: ";
                for (int node : path) {
                    cout << node << " ";
                }
                cout << "\n";
            } else {
                cout << "Vehicle " << i << ": No route found.\n";
            }
        }
        
        // Update edge loads and then update edge weights based on congestion.
        update_edge_loads(graph, routes);
        update_edge_weights(graph);
        
        // Check for convergence: if routes haven't changed.
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
        
        // Save current routes for the next iteration comparison.
        prev_routes = routes;
        iteration++;
        
        // Display updated edge information.
        cout << "Updated edge weights due to congestion:\n";
        for (int i = 0; i < graph.edges.size(); i++) {
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
        for (int node : routes[i]) {
            cout << node << " ";
        }
        cout << "\n";
    }
    
    return 0;
}
