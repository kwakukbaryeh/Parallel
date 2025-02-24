/**
 * Parallel VLSI Wire Routing via OpenMP
 * Kwaku Baryeh(kbaryeh), Parth Iyer(pniyer)
 */

#include "wireroute.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <string>
#include <vector>

#include <unistd.h>
#include <omp.h>

void print_stats(const std::vector<std::vector<int>>& occupancy) {
  int max_occupancy = 0;
  long long total_cost = 0;

  for (const auto& row : occupancy) {
    for (const int count : row) {
      max_occupancy = std::max(max_occupancy, count);
      total_cost += count * count;
    }
  }

  std::cout << "Max occupancy: " << max_occupancy << '\n';
  std::cout << "Total cost: " << total_cost << '\n';
}

void write_output(const std::vector<Wire>& wires, const int num_wires, const std::vector<std::vector<int>>& occupancy, const int dim_x, const int dim_y, const int num_threads, std::string input_filename) {
  if (std::size(input_filename) >= 4 && input_filename.substr(std::size(input_filename) - 4) == ".txt") {
    input_filename.resize(std::size(input_filename) - 4);
  }

  const std::string occupancy_filename = input_filename + "_occupancy_" + std::to_string(num_threads) + ".txt";
  const std::string wires_filename = input_filename + "_wires_" + std::to_string(num_threads) + ".txt";

  std::ofstream out_occupancy(occupancy_filename, std::fstream::out);
  if (!out_occupancy) {
    std::cerr << "Unable to open file: " << occupancy_filename << '\n';
    exit(EXIT_FAILURE);
  }

  out_occupancy << dim_x << ' ' << dim_y << '\n';
  for (const auto& row : occupancy) {
    for (const int count : row) {
      out_occupancy << count << ' ';
    }
    out_occupancy << '\n';
  }

  out_occupancy.close();

  std::ofstream out_wires(wires_filename, std::fstream:: out);
  if (!out_wires) {
    std::cerr << "Unable to open file: " << wires_filename << '\n';
    exit(EXIT_FAILURE);
  }

  out_wires << dim_x << ' ' << dim_y << '\n' << num_wires << '\n';

  for (const auto& [start_x, start_y, end_x, end_y, bend1_x, bend1_y] : wires) {
    out_wires << start_x << ' ' << start_y << ' ' << bend1_x << ' ' << bend1_y << ' ';

    if (start_y == bend1_y) {
    // first bend was horizontal

      if (end_x != bend1_x) {
        // two bends

        out_wires << bend1_x << ' ' << end_y << ' ';
      }
    } else if (start_x == bend1_x) {
      // first bend was vertical

      if(end_y != bend1_y) {
        // two bends

        out_wires << end_x << ' ' << bend1_y << ' ';
      }
    }
    out_wires << end_x << ' ' << end_y << '\n';
  }

  out_wires.close();
}

void route_wires_within(std::vector<Wire>& wires, std::vector<std::vector<int>>& occupancy, 
  int dim_x, int dim_y, int num_threads, double SA_prob, int SA_iters) {
    std::cout << "route_wires_within() is not implemented yet.\n";
}

void update_occupancy(std::vector<std::vector<int>>& occupancy, const std::vector<std::pair<int, int>>& route) {
  for (const auto& point : route) {
      int x = point.first, y = point.second;  // Correct way to access x and y
      occupancy[y][x]++;  // Increment occupancy at this point
  }
}

std::vector<std::pair<int, int>> find_best_route(const Wire& wire, 
                                                 const std::vector<std::vector<int>>& occupancy, 
                                                 int dim_x, int dim_y, double SA_prob) {
    int x = wire.start_x, y = wire.start_y;
    int end_x = wire.end_x, end_y = wire.end_y;
    
    std::vector<std::pair<int, int>> route;
    route.push_back({x, y});  // Start point
    
    // Simple greedy approach with simulated annealing
    while (x != end_x || y != end_y) {
        std::vector<std::pair<int, int>> possible_moves;

        // Generate possible moves (right, left, up, down)
        if (x < end_x) possible_moves.push_back({x + 1, y});  // Move right
        if (x > end_x) possible_moves.push_back({x - 1, y});  // Move left
        if (y < end_y) possible_moves.push_back({x, y + 1});  // Move down
        if (y > end_y) possible_moves.push_back({x, y - 1});  // Move up

        // Sort moves by occupancy to prefer less crowded paths
        std::sort(possible_moves.begin(), possible_moves.end(), [&](const auto& a, const auto& b) {
            return occupancy[a.second][a.first] < occupancy[b.second][b.first];
        });

        // Apply simulated annealing: occasionally pick a non-optimal move
        if (!possible_moves.empty()) {
            if ((rand() / (double)RAND_MAX) < SA_prob) {
                int random_index = rand() % possible_moves.size();
                x = possible_moves[random_index].first;
                y = possible_moves[random_index].second;
            } else {
                x = possible_moves[0].first;
                y = possible_moves[0].second;
            }
            route.push_back({x, y});
        } else {
            // If no valid moves, break (this should not happen in a valid grid)
            break;
        }
    }

    return route;
}

void route_wires_across(std::vector<Wire>& wires, std::vector<std::vector<int>>& occupancy, 
                        int dim_x, int dim_y, int num_threads, int batch_size, 
                        double SA_prob, int SA_iters) {
    #pragma omp parallel num_threads(num_threads)
    {
        for (int timestep = 0; timestep < SA_iters; timestep++) {
            #pragma omp for schedule(dynamic)
            for (int batch_start = 0; batch_start < wires.size(); batch_start += batch_size) {
                std::vector<std::vector<std::pair<int, int>>> batch_routes(batch_size);  // Corrected type

                // Step 1: Compute new routes for batch wires
                for (int i = 0; i < batch_size && (batch_start + i) < wires.size(); i++) {
                    int wire_idx = batch_start + i;
                    batch_routes[i] = find_best_route(wires[wire_idx], occupancy, dim_x, dim_y, SA_prob);
                }

                // Step 2: Synchronize update of occupancy matrix
                #pragma omp critical
                {
                    for (int i = 0; i < batch_size && (batch_start + i) < wires.size(); i++) {
                        update_occupancy(occupancy, batch_routes[i]); // Now correctly passing std::vector<std::pair<int, int>>
                    }
                }
            }
        }
    }
}

int main(int argc, char *argv[]) {
  const auto init_start = std::chrono::steady_clock::now();

  std::string input_filename;
  int num_threads = 0;
  double SA_prob = 0.1;
  int SA_iters = 5;
  char parallel_mode = '\0';
  int batch_size = 1;

  int opt;
  while ((opt = getopt(argc, argv, "f:n:p:i:m:b:")) != -1) {
    switch (opt) {
      case 'f':
        input_filename = optarg;
        break;
      case 'n':
        num_threads = atoi(optarg);
        break;
      case 'p':
        SA_prob = atof(optarg);
        break;
      case 'i':
        SA_iters = atoi(optarg);
        break;
      case 'm':
        parallel_mode = *optarg;
        break;
      case 'b':
        batch_size = atoi(optarg);
        break;
      default:
        std::cerr << "Usage: " << argv[0] << " -f input_filename -n num_threads [-p SA_prob] [-i SA_iters] -m parallel_mode -b batch_size\n";
        exit(EXIT_FAILURE);
    }
  }

  // Check if required options are provided
  if (empty(input_filename) || num_threads <= 0 || SA_iters <= 0 || (parallel_mode != 'A' && parallel_mode != 'W') || batch_size <= 0) {
    std::cerr << "Usage: " << argv[0] << " -f input_filename -n num_threads [-p SA_prob] [-i SA_iters] -m parallel_mode -b batch_size\n";
    exit(EXIT_FAILURE);
  }

  std::cout << "Number of threads: " << num_threads << '\n';
  std::cout << "Simulated annealing probability parameter: " << SA_prob << '\n';
  std::cout << "Simulated annealing iterations: " << SA_iters << '\n';
  std::cout << "Input file: " << input_filename << '\n';
  std::cout << "Parallel mode: " << parallel_mode << '\n';
  std::cout << "Batch size: " << batch_size << '\n';

  std::ifstream fin(input_filename);

  if (!fin) {
    std::cerr << "Unable to open file: " << input_filename << ".\n";
    exit(EXIT_FAILURE);
  }

  int dim_x, dim_y;
  int num_wires;

  /* Read the grid dimension and wire information from file */
  fin >> dim_x >> dim_y >> num_wires;

  std::vector<Wire> wires(num_wires);
  std::vector occupancy(dim_y, std::vector<int>(dim_x));

  for (auto& wire : wires) {
    fin >> wire.start_x >> wire.start_y >> wire.end_x >> wire.end_y;
    wire.bend1_x = wire.start_x;
    wire.bend1_y = wire.start_y;
  }

  /* Initialize any additional data structures needed in the algorithm */

  const double init_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - init_start).count();
  std::cout << "Initialization time (sec): " << std::fixed << std::setprecision(10) << init_time << '\n';

  const auto compute_start = std::chrono::steady_clock::now();

  /** 
   * Implement the wire routing algorithm here
   * Feel free to structure the algorithm into different functions
   * Don't use global variables.
   * Use OpenMP to parallelize the algorithm. 
   */
   if (parallel_mode == 'A') {
    route_wires_across(wires, occupancy, dim_x, dim_y, num_threads, batch_size, SA_prob, SA_iters);
  } else {
    route_wires_within(wires, occupancy, dim_x, dim_y, num_threads, SA_prob, SA_iters);
  }

  const double compute_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - compute_start).count();
  std::cout << "Computation time (sec): " << compute_time << '\n';

  /* Write wires and occupancy matrix to files */

  print_stats(occupancy);
  write_output(wires, num_wires, occupancy, dim_x, dim_y, num_threads, input_filename);
}

validate_wire_t Wire::to_validate_format(void) const {
  /* TODO(student): Implement this if you want to use the wr_checker. */
  /* See wireroute.h for details on validate_wire_t. */
  throw std::logic_error("to_validate_format not implemented.");
}
