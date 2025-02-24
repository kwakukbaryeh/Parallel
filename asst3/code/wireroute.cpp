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

    int x_dir = start.x - end.x > 0 ? -1 : 1;
    int y_dir = start.y - end.y > 0 ? -1 : 1;

    if (start.x == bend1.x) { 
        for (int y = start.y; y_dir < 0 ? bend1.y < y : y < bend1.y; y += y_dir)
            path.push_back({start.x, y});
    } else { 
        for (int x = start.x; x_dir < 0 ? bend1.x < x : x < bend1.x; x += x_dir)
            path.push_back({x, start.y});
    }

    if (bend1.x == bend2.x) { 
        for (int y = bend1.y; y_dir < 0 ? bend2.y < y : y < bend2.y; y += y_dir)
            path.push_back({bend1.x, y});
    } else { 
        for (int x = bend1.x; x_dir < 0 ? bend2.x < x : x < bend2.x; x += x_dir)
            path.push_back({x, bend1.y});
    }

    if (bend2.x == end.x) { 
        for (int y = bend2.y; y_dir < 0 ? end.y < y : y < end.y; y += y_dir)
            path.push_back({bend2.x, y});
    } else { 
        for (int x = bend2.x; x_dir < 0 ? end.x < x : x < end.x; x += x_dir)
            path.push_back({x, bend2.y});
    }

    path.push_back(end);
    return path;
}

//Route Wires per wire parallelism
void w_route(std::vector<Wire> &w, std::vector<std::vector<int>> &grid, double SA_prob, Point dim, int t) {
    for(size_t i = 0; i < w.size(); i++) {
        //Remove wire from grid if present
        if (w[i].start_x != w[i].bend1_x || w[i].start_y != w[i].bend1_y) {
            std::vector<Point> path = get_path(w[i]);
            for (Point p : path)
                grid[p.y][p.x]--;
        }

        int minCost = INT_MAX;
        Wire minWire = w[i];

        if (((double) random())/((double) RAND_MAX) < SA_prob) {
            Point dirs = {w[i].end_x - w[i].start_x, w[i].end_y - w[i].start_y};
            if (((double) random())/((double) RAND_MAX) > 0.5) {
                //Vertical
                w[i].bend1_y = (int) (dirs.y * ((double) random())/((double) RAND_MAX)) + 1;
            } else {
                //Horizontal
                w[i].bend1_x = (int) (dirs.x * ((double) random())/((double) RAND_MAX)) + 1;
            }
        } else {
            #pragma omp parallel num_threads(t)
            {
                //Check Horizontal
                #pragma omp for
                for(int x = std::min(w[i].start_x, w[i].end_x); x <= std::max(w[i].start_x, w[i].end_x); ++x) {
                    Wire wire = w[i];
                    wire.bend1_x = x;
                    wire.bend1_y = w[i].start_y;
                    std::vector<Point> path = get_path(wire);
                    int cost = 0;
                    for(Point p : path)
                        cost += grid[p.y][p.x] * grid[p.y][p.x];
                    if (cost < minCost) {
                        minWire = wire;
                        minCost = cost;
                    }
                }

                //Check Vertical
                #pragma omp for
                for(int y = std::min(w[i].start_y, w[i].end_y); y <= std::max(w[i].start_y, w[i].end_y); ++y) {
                    Wire wire = w[i];
                    wire.bend1_x = w[i].start_x;
                    wire.bend1_y = y;
                    std::vector<Point> path = get_path(wire);
                    int cost = 0;
                    for(Point p : path)
                        cost += grid[p.y][p.x] * grid[p.y][p.x];
                    if (cost < minCost) {
                        minWire = wire;
                        minCost = cost;
                    }
                }
            }
        }

        //Add new wire to grid and wires
        w[i] = minWire;
        std::vector<Point> path = get_path(minWire);
        for(Point p : path)
            grid[p.y][p.x]++;
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
    if (parallel_mode == 'W') {
        for (int i = 0; i < SA_iters; i++)
            w_route(wires, occupancy, SA_prob, {dim_x, dim_y}, num_threads);
        printf("DEBUG GRID\n");
        for(int i = 0; i < dim_y; i++) {
            for (int j = 0; j < dim_x; j++) {
                printf("%d", occupancy[i][j]);
            }
            printf("\n");
        }
    } else /*parallel_mode == 'A'*/ {
      
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
