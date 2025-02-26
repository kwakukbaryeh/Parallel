/**
 * Parallel VLSI Wire Routing via OpenMP
 * Name 1(andrew_id 1), Name 2(andrew_id 2)
 */

 #include "wireroute.h"

 #include <algorithm>
 #include <iostream>
 #include <fstream>
 #include <iomanip>
 #include <chrono>
 #include <string>
 #include <vector>
 #include <climits>
 
 #include <unistd.h>
 #include <omp.h>
 #include <stdlib.h>
 
 // Point Data Structure
 struct Point {
     int x,y;
 };
 
 bool operator==(const Point &a, const Point &b) {
    return a.x == b.x && a.y == b.y;
}

 void print_point(Point p) {
     printf("(%d,%d)", p.x, p.y);
 }
 
 void print_wire(Wire wire) {
     printf("((%d,%d),(%d,%d),(%d,%d))", wire.start_x, wire.start_y, wire.bend1_x, wire.bend1_y, wire.end_x, wire.end_y);
 }
 
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
 
 //Get Path for a wire
 std::vector<Point> get_path(Wire w) {
     std::vector<Point> path;
     Point start = {w.start_x, w.start_y};
     Point bend1 = {w.bend1_x, w.bend1_y};
     Point bend2;
 
     if (w.start_x == w.bend1_x) {
         bend2.x = w.end_x;
         bend2.y = w.bend1_y;
     } else {
         bend2.x = w.bend1_x;
         bend2.y = w.end_y;
     }
 
     Point end = {w.end_x, w.end_y};
     int x_dir = (start.x - end.x > 0) ? -1 : 1;
     int y_dir = (start.y - end.y > 0) ? -1 : 1;
 
     if (start.x == bend1.x) {
        // Vertical movement
        for (int y = start.y; (y_dir < 0 ? y >= bend1.y : y <= bend1.y); y += y_dir)
            path.push_back({start.x, y});
    } else {
        // Horizontal movement
        for (int x = start.x; (x_dir < 0 ? x >= bend1.x : x <= bend1.x); x += x_dir)
            path.push_back({x, start.y});
    }
    
    // Movement from bend1 to bend2
    if (bend1.x == bend2.x) {
        // Vertical movement
        for (int y = bend1.y; (y_dir < 0 ? y >= bend2.y : y <= bend2.y); y += y_dir)
            path.push_back({bend1.x, y});
    } else {
        // Horizontal movement
        for (int x = bend1.x; (x_dir < 0 ? x >= bend2.x : x <= bend2.x); x += x_dir)
            path.push_back({x, bend1.y});
    }
    
    // Movement from bend2 to end
    if (bend2.x == end.x) {
        // Vertical movement
        for (int y = bend2.y; (y_dir < 0 ? y >= end.y : y <= end.y); y += y_dir)
            path.push_back({bend2.x, y});
    } else {
        // Horizontal movement
        for (int x = bend2.x; (x_dir < 0 ? x >= end.x : x <= end.x); x += x_dir)
            path.push_back({x, bend2.y});
    }
  
 
     path.push_back(end);
     return path;
 }
 
 //Write the wire to the occupancy matrix
 static inline void write_wire(Wire w, int v, std::vector<std::vector<int>>& grid) {
     Point start = {w.start_x, w.start_y};
     Point bend1 = {w.bend1_x, w.bend1_y};
     Point bend2;
     if (w.start_x == w.bend1_x) {
         bend2.x = w.end_x;
         bend2.y = w.bend1_y;
     } else {
         bend2.x = w.bend1_x;
         bend2.y = w.end_y;
     }
     Point end = {w.end_x, w.end_y};
     int x_dir = (start.x - end.x > 0) ? -1 : 1;
     int y_dir = (start.y - end.y > 0) ? -1 : 1;
 
     if (start.x == bend1.x) {
         for (int y = start.y; (y_dir < 0 ? y > bend1.y : y < bend1.y); y += y_dir)
             grid[y][start.x] += v;
     } else {
         for (int x = start.x; (x_dir < 0 ? x > bend1.x : x < bend1.x); x += x_dir)
             grid[start.y][x] += v;
     }
     if (bend1.x == bend2.x) {
         for (int y = bend1.y; (y_dir < 0 ? y > bend2.y : y < bend2.y); y += y_dir)
             grid[y][bend1.x] += v;
     } else {
         for (int x = bend1.x; (x_dir < 0 ? x > bend2.x : x < bend2.x); x += x_dir)
             grid[bend1.y][x] += v;
     }
     if (bend2.x == end.x) {
         for (int y = bend2.y; (y_dir < 0 ? y > end.y : y < end.y); y += y_dir)
             grid[y][bend2.x] += v;
     } else {
         for (int x = bend2.x; (x_dir < 0 ? x > end.x : x < end.x); x += x_dir)
             grid[bend2.y][x] += v;
     }
     grid[end.y][end.x] += v;
 }
 
 //Get the cost of a wire route
 static inline int get_cost(Wire w, std::vector<std::vector<int>>& grid) {
     int cost = 0;
     int sq;
     Point start = {w.start_x, w.start_y};
     Point bend1 = {w.bend1_x, w.bend1_y};
     Point bend2;
     if (w.start_x == w.bend1_x) {
         bend2.x = w.end_x;
         bend2.y = w.bend1_y;
     } else {
         bend2.x = w.bend1_x;
         bend2.y = w.end_y;
     }
     Point end = {w.end_x, w.end_y};
     int x_dir = (start.x - end.x > 0) ? -1 : 1;
     int y_dir = (start.y - end.y > 0) ? -1 : 1;
 
     if (start.x == bend1.x) {
         for (int y = start.y; (y_dir < 0 ? y > bend1.y : y < bend1.y); y += y_dir) {
             sq = grid[y][start.x];
             cost += sq * sq;
         }
     } else {
         for (int x = start.x; (x_dir < 0 ? x > bend1.x : x < bend1.x); x += x_dir) {
             sq = grid[start.y][x];
             cost += sq * sq;
         }
     }
     if (bend1.x == bend2.x) {
         for (int y = bend1.y; (y_dir < 0 ? y > bend2.y : y < bend2.y); y += y_dir) {
             sq = grid[y][bend1.x];
             cost += sq * sq;
         }
     } else {
         for (int x = bend1.x; (x_dir < 0 ? x > bend2.x : x < bend2.x); x += x_dir) {
             sq = grid[bend1.y][x];
             cost += sq * sq;
         }
     }
     if (bend2.x == end.x) {
         for (int y = bend2.y; (y_dir < 0 ? y > end.y : y < end.y); y += y_dir) {
             sq = grid[y][bend2.x];
             cost += sq * sq;
         }
     } else {
         for (int x = bend2.x; (x_dir < 0 ? x > end.x : x < end.x); x += x_dir) {
             sq = grid[bend2.y][x];
             cost += sq * sq;
         }
     }
     sq = grid[end.y][end.x];
     cost += sq * sq;
     return cost;
 }
 
 // A simple linear search to check if a point appears in a vector
 bool point_in_vector(Point p, const std::vector<Point>& vec) {
   for (const auto& q : vec) {
       if (p.x == q.x && p.y == q.y) return true;
   }
   return false;
 }
 
 // Compute the cost of candidate wire route assuming the old route is "removed".
 // For each point in the candidate route, if it was occupied by the old route, subtract one.
 int get_cost_adjusted(Wire w, const std::vector<Point>& old_path, const std::vector<std::vector<int>>& grid) {
     int cost = 0;
     std::vector<Point> candidate_path = get_path(w);
     for (const auto& p : candidate_path) {
         int occ = grid[p.y][p.x];
         if (point_in_vector(p, old_path)) occ--;
         cost += occ * occ;
     }
     return cost;
 }
 
 // The across-wires routing function. For each SA iteration, process wires in batches.
 // Within a batch, first determine a new route for each wire (without updating the occupancy matrix),
 // then update the occupancy matrix atomically for each wire in the batch.
 void a_route(std::vector<Wire>& w, std::vector<std::vector<int>>& grid, double SA_prob, Point dim, int SA_iters, int batch_size) {
     int num_wires = w.size();
     for (int t = 0; t < SA_iters; t++) {
         // Process wires in batches in parallel.
         #pragma omp parallel for schedule(dynamic)
         for (int batch_start = 0; batch_start < num_wires; batch_start += batch_size) {
             int batch_end = std::min(batch_start + batch_size, num_wires);
             std::vector<Wire> new_candidates(batch_end - batch_start);
             // Candidate evaluation phase:
             for (int j = batch_start; j < batch_end; j++) {
                 Wire old_wire = w[j];
                 std::vector<Point> old_path = get_path(old_wire);
                 Wire best_wire = old_wire;
                 int minCost = get_cost_adjusted(old_wire, old_path, grid);
                 double r = ((double)rand()) / RAND_MAX;
                 if (r < SA_prob) {
                     // Choose a random candidate.
                     best_wire = old_wire;
                     double r2 = ((double)rand()) / RAND_MAX;
                     if (r2 > 0.5) {
                         int range = std::abs(old_wire.end_y - old_wire.start_y);
                         if (range > 0)
                             best_wire.bend1_y = old_wire.start_y + (int)(range * (((double)rand()) / RAND_MAX));
                     } else {
                         int range = std::abs(old_wire.end_x - old_wire.start_x);
                         if (range > 0)
                             best_wire.bend1_x = old_wire.start_x + (int)(range * (((double)rand()) / RAND_MAX));
                     }
                 } else {
                     // Deterministically search for a better route.
                     if (old_wire.start_x != old_wire.end_x) {
                         int sign = (old_wire.end_x - old_wire.start_x > 0) ? 1 : -1;
                         for (int x = old_wire.start_x + sign; (sign > 0 ? x <= old_wire.end_x : x >= old_wire.end_x); x += sign) {
                             Wire temp = old_wire;
                             temp.bend1_x = x;
                             temp.bend1_y = old_wire.start_y;
                             int cost = get_cost_adjusted(temp, old_path, grid);
                             if (cost < minCost) {
                                 minCost = cost;
                                 best_wire = temp;
                             }
                         }
                     }
                     if (old_wire.start_y != old_wire.end_y) {
                         int sign = (old_wire.end_y - old_wire.start_y > 0) ? 1 : -1;
                         for (int y = old_wire.start_y + sign; (sign > 0 ? y <= old_wire.end_y : y >= old_wire.end_y); y += sign) {
                             Wire temp = old_wire;
                             temp.bend1_x = old_wire.start_x;
                             temp.bend1_y = y;
                             int cost = get_cost_adjusted(temp, old_path, grid);
                             if (cost < minCost) {
                                 minCost = cost;
                                 best_wire = temp;
                             }
                         }
                     }
                 }
                 new_candidates[j - batch_start] = best_wire;
             }
             // Update phase: For each wire in the batch, remove the old route then add the new route.
             for (int j = batch_start; j < batch_end; j++) {
                 Wire old_wire = w[j];
                 Wire candidate = new_candidates[j - batch_start];
                 std::vector<Point> old_path = get_path(old_wire);
                 std::vector<Point> new_path = get_path(candidate);
                 if (old_path != new_path) {
                  // Atomically remove old route.
                  for (const auto& p : old_path) {
                      #pragma omp atomic
                      grid[p.y][p.x] -= 1;
                  }
                  // Atomically add new route.
                  for (const auto& p : new_path) {
                      #pragma omp atomic
                      grid[p.y][p.x] += 1;
                  }
                }
              w[j] = candidate;
            }
         } // end of batch loop
     } // end SA iteration loop
 }
 
 // The sequential “within wires” routing (provided as baseline)
 void s_route(std::vector<Wire> &w, std::vector<std::vector<int>> &grid, double SA_prob, Point dim, int t) {
     for (size_t i = 0; i < w.size(); i++) {
         // Remove the old route if it exists.
         if (w[i].start_x != w[i].bend1_x || w[i].start_y != w[i].bend1_y)
             write_wire(w[i], -1, grid);
 
         int minCost = INT_MAX;
         Wire minWire = w[i];
 
         if (((double)rand())/RAND_MAX < SA_prob) {
             Point dirs = {w[i].end_x - w[i].start_x, w[i].end_y - w[i].start_y};
             if (((double)rand())/RAND_MAX > 0.5) {
                 w[i].bend1_y = (int)(dirs.y * (((double)rand())/RAND_MAX)) + 1;
             } else {
                 w[i].bend1_x = (int)(dirs.x * (((double)rand())/RAND_MAX)) + 1;
             }
         } else {
             Point dir = {(w[i].end_x - w[i].start_x > 0 ? 1 : -1), (w[i].end_y - w[i].start_y > 0 ? 1 : -1)};
             int cost;
             if (dir.x > 0) {
                 for (int x = w[i].start_x + dir.x; x <= w[i].end_x; x += dir.x) {
                     Wire wire = w[i];
                     wire.bend1_x = x;
                     wire.bend1_y = w[i].start_y;
                     cost = get_cost(wire, grid);
                     if (cost < minCost) {
                         minWire = wire;
                         minCost = cost;
                     }
                 }
             } else {
                 for (int x = w[i].start_x + dir.x; w[i].end_x <= x; x += dir.x) {
                     Wire wire = w[i];
                     wire.bend1_x = x;
                     wire.bend1_y = w[i].start_y;
                     cost = get_cost(wire, grid);
                     if (cost < minCost) {
                         minWire = wire;
                         minCost = cost;
                     }
                 }
             }
             if (dir.y > 0) {
                 for (int y = w[i].start_y + dir.y; y <= w[i].end_y; y += dir.y) {
                     Wire wire = w[i];
                     wire.bend1_x = w[i].start_x;
                     wire.bend1_y = y;
                     cost = get_cost(wire, grid);
                     if (cost < minCost) {
                         minWire = wire;
                         minCost = cost;
                     }
                 }
             } else {
                 for (int y = w[i].start_y + dir.y; w[i].end_y <= y; y += dir.y) {
                     Wire wire = w[i];
                     wire.bend1_x = w[i].start_x;
                     wire.bend1_y = y;
                     cost = get_cost(wire, grid);
                     if (cost < minCost) {
                         minWire = wire;
                         minCost = cost;
                     }
                 }
             }
         }
         w[i] = minWire;
         write_wire(w[i], 1, grid);
     }
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
 
 int main(int argc, char *argv[]) {
     srandom(0);
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
             case 'f': input_filename = optarg; break;
             case 'n': num_threads = atoi(optarg); break;
             case 'p': SA_prob = atof(optarg); break;
             case 'i': SA_iters = atoi(optarg); break;
             case 'm': parallel_mode = *optarg; break;
             case 'b': batch_size = atoi(optarg); break;
             default:
                 std::cerr << "Usage: " << argv[0]
                           << " -f input_filename -n num_threads [-p SA_prob] [-i SA_iters] -m parallel_mode -b batch_size\n";
                 exit(EXIT_FAILURE);
         }
     }
 
     if (input_filename.empty() || num_threads <= 0 || SA_iters <= 0 ||
         (parallel_mode != 'A' && parallel_mode != 'W') || batch_size <= 0) {
         std::cerr << "Usage: " << argv[0]
                   << " -f input_filename -n num_threads [-p SA_prob] [-i SA_iters] -m parallel_mode -b batch_size\n";
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
 
     int dim_x, dim_y, num_wires;
     fin >> dim_x >> dim_y >> num_wires;
     std::vector<Wire> wires(num_wires);
     std::vector<std::vector<int>> occupancy(dim_y, std::vector<int>(dim_x, 0));
 
     for (auto& wire : wires) {
         fin >> wire.start_x >> wire.start_y >> wire.end_x >> wire.end_y;
         wire.bend1_x = wire.start_x;
         wire.bend1_y = wire.start_y;
     }
 
     const double init_time = std::chrono::duration_cast<std::chrono::duration<double>>(
                                  std::chrono::steady_clock::now() - init_start).count();
     std::cout << "Initialization time (sec): " << std::fixed << std::setprecision(10) << init_time << '\n';
 
     const auto compute_start = std::chrono::steady_clock::now();
 
     if (num_threads == 1) {
         for (int i = 0; i < SA_iters; i++) {
             a_route(wires, occupancy, SA_prob, {dim_x, dim_y}, SA_iters, num_threads);
         }
     } else if (parallel_mode == 'W') {
         // Within-wires approach (not shown here)
     } else { // parallel_mode == 'A'
         omp_set_num_threads(num_threads);
         a_route(wires, occupancy, SA_prob, {dim_x, dim_y}, SA_iters, batch_size);
     }
 
     const double compute_time = std::chrono::duration_cast<std::chrono::duration<double>>(
                                     std::chrono::steady_clock::now() - compute_start).count();
     std::cout << "Computation time (sec): " << compute_time << '\n';
 
     print_stats(occupancy);
     write_output(wires, num_wires, occupancy, dim_x, dim_y, num_threads, input_filename);
     return 0;
 }
 
 validate_wire_t Wire::to_validate_format(void) const {
     throw std::logic_error("to_validate_format not implemented.");
 }
 
 