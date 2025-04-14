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

#define GRID_SIZE 128
 
// Point Data Structure
struct Point {
     int x,y;
};

// Wire cost Structure
struct WireCost {
    int c;
    Wire w;
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

    int x_dir = start.x - end.x > 0 ? -1 : 1;
    int y_dir = start.y - end.y > 0 ? -1 : 1;

    if (start.x == bend1.x) { 
        for (int y = start.y; y_dir < 0 ? bend1.y < y : y < bend1.y; y += y_dir)
            grid[y][start.x] += v;

        for (int x = bend1.x; x_dir < 0 ? bend2.x < x : x < bend2.x; x += x_dir)
            grid[bend1.y][x] += v;

        for (int y = bend2.y; y_dir < 0 ? end.y < y : y < end.y; y += y_dir)
            grid[y][bend2.x] += v;
    } else { 
        for (int x = start.x; x_dir < 0 ? bend1.x < x : x < bend1.x; x += x_dir)
            grid[start.y][x] += v;

        for (int y = bend1.y; y_dir < 0 ? bend2.y < y : y < bend2.y; y += y_dir)
            grid[y][bend1.x] += v;

        for (int x = bend2.x; x_dir < 0 ? end.x < x : x < end.x; x += x_dir)
            grid[bend2.y][x] += v;
    }

    grid[end.y][end.x] += v;
}

//Write the wire to the occupancy matrix
static void lock_write_wire(Wire w, int v, std::vector<std::vector<int>>& grid, std::vector<std::vector<omp_lock_t>>& locks) { //Fix by passing locks by reference 
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

    Point lock_coords = {start.x/GRID_SIZE, start.y/GRID_SIZE};
    omp_set_lock(&locks[lock_coords.y][lock_coords.x]);
    if (start.x == bend1.x) { 
        for (int y = start.y; y_dir < 0 ? bend1.y < y : y < bend1.y; y += y_dir) {
            if (y/GRID_SIZE != lock_coords.y) {
                omp_unset_lock(&locks[lock_coords.y][lock_coords.x]);
                lock_coords.y = y/GRID_SIZE;
                omp_set_lock(&locks[lock_coords.y][lock_coords.x]);
            }
            grid[y][start.x] += v;
        }

        for (int x = bend1.x; x_dir < 0 ? bend2.x < x : x < bend2.x; x += x_dir) {
            if (x/GRID_SIZE != lock_coords.x) {
                omp_unset_lock(&locks[lock_coords.y][lock_coords.x]);
                lock_coords.x = x/GRID_SIZE;
                omp_set_lock(&locks[lock_coords.y][lock_coords.x]);
            }
            grid[bend1.y][x] += v;
        }

        for (int y = bend2.y; y_dir < 0 ? end.y < y : y < end.y; y += y_dir) {
            if (y/GRID_SIZE != lock_coords.y) {
                omp_unset_lock(&locks[lock_coords.y][lock_coords.x]);
                lock_coords.y = y/GRID_SIZE;
                omp_set_lock(&locks[lock_coords.y][lock_coords.x]);
            }
            grid[y][bend2.x] += v;
        }
    } else { 
        for (int x = start.x; x_dir < 0 ? bend1.x < x : x < bend1.x; x += x_dir) {
            if (x/GRID_SIZE != lock_coords.x) {
                omp_unset_lock(&locks[lock_coords.y][lock_coords.x]);
                lock_coords.x = x/GRID_SIZE;
                omp_set_lock(&locks[lock_coords.y][lock_coords.x]);
            }
            grid[start.y][x] += v;
        }

        for (int y = bend1.y; y_dir < 0 ? bend2.y < y : y < bend2.y; y += y_dir) {
            if (y/GRID_SIZE != lock_coords.y) {
                omp_unset_lock(&locks[lock_coords.y][lock_coords.x]);
                lock_coords.y = y/GRID_SIZE;
                omp_set_lock(&locks[lock_coords.y][lock_coords.x]);
            }
            grid[y][bend1.x] += v;
        }

        for (int x = bend2.x; x_dir < 0 ? end.x < x : x < end.x; x += x_dir) {
            if (x/GRID_SIZE != lock_coords.x) {
                omp_unset_lock(&locks[lock_coords.y][lock_coords.x]);
                lock_coords.x = x/GRID_SIZE;
                omp_set_lock(&locks[lock_coords.y][lock_coords.x]);
            }
            grid[bend2.y][x] += v;
        }
    }

    grid[end.y][end.x] += v;
    omp_unset_lock(&locks[lock_coords.y][lock_coords.x]);
}
 
//Get the cost of a horizontal wire route
static int get_cost_horizontal(Wire w, std::vector<std::vector<int>>& grid, Point dir) {
    Point bend2;

    if (w.start_x == w.bend1_x) {
        bend2.x = w.end_x;
        bend2.y = w.bend1_y;
    } else {
        bend2.x = w.bend1_x;
        bend2.y = w.end_y;
    }

    int cost = 0;
    int sq;

    for (int y = w.start_y; dir.y < 0 ? w.bend1_y < y : y < w.bend1_y; y += dir.y){
        sq = grid[y][w.start_x] + 1;
        cost += sq * sq;
    }

    for (int x = w.bend1_x; dir.x < 0 ? bend2.x < x : x < bend2.x; x += dir.x) {
        sq = grid[w.bend1_y][x] + 1;
        cost += sq * sq;
    }

    for (int y = bend2.y; dir.y < 0 ? w.end_y < y : y < w.end_y; y += dir.y) {
        sq = grid[y][bend2.x] + 1;
        cost += sq * sq;
    }

    sq = grid[w.end_y][w.end_x] + 1;
    cost += sq * sq;
    return cost;
}

//Get the cost of a vertical wire route
static int get_cost_vertical(Wire w, std::vector<std::vector<int>>& grid, Point dir) {
    Point bend2;

    if (w.start_x == w.bend1_x) {
        bend2.x = w.end_x;
        bend2.y = w.bend1_y;
    } else {
        bend2.x = w.bend1_x;
        bend2.y = w.end_y;
    }

    int cost = 0;
    int sq;

    for (int x = w.start_x; dir.x < 0 ? w.bend1_x < x : x < w.bend1_x; x += dir.x) {
        sq = grid[w.start_y][x] + 1;
        cost += sq * sq;
    }

    for (int y = w.bend1_y; dir.y < 0 ? bend2.y < y : y < bend2.y; y += dir.y) {
        sq = grid[y][w.bend1_x] + 1;
        cost += sq * sq;
    }

    for (int x = bend2.x; dir.x < 0 ? w.end_x < x : x < w.end_x; x += dir.x) {
        sq = grid[bend2.y][x] + 1;
        cost += sq * sq;
    }

    sq = grid[w.end_y][w.end_x] + 1;
    cost += sq * sq;
    return cost;
}

int randInt(int min, int max) {
    if (min > max) {
        int t = min;
        min = max;
        max = t;
    }

    int range = max - min + 1;

    return rand() % range + min;;
}
 
//Route Wires sequentially
void s_route(std::vector<Wire> &w, std::vector<std::vector<int>> &grid, double SA_prob, Point dim, int t, int it) {
    for(size_t i = 0; i < w.size(); i++) {
        //Remove wire from grid if present
        if (w[i].start_x == w[i].end_x || w[i].start_y == w[i].end_y) {
            if(it == 0)
                write_wire(w[i], 1, grid);
            continue;
        } else if (it > 0) {
            write_wire(w[i], -1, grid);
        }

        int minCost = INT_MAX;
        Wire minWire = w[i];

        if (it > 0 && (double) random() / RAND_MAX < SA_prob) {
            if ((double) random() / RAND_MAX < SA_prob) {
                //Vertical
                w[i].bend1_x = w[i].start_x;
                w[i].bend1_y = randInt(w[i].start_y, w[i].end_y);
            } else {
                //Horizontal
                w[i].bend1_x = randInt(w[i].start_x, w[i].end_x);
                w[i].bend1_y = w[i].start_y;
            }
            write_wire(w[i], 1, grid);
            continue;
        } else {
            Point dir = {(w[i].end_x - w[i].start_x > 0 ? 1 : -1), (w[i].end_y - w[i].start_y > 0 ? 1 : -1)};
            int cost;
            //Check Vertical
            if (dir.x > 0) {
                for(int x = w[i].start_x + dir.x; x <= w[i].end_x; x += dir.x) {
                    Wire wire = w[i];
                    wire.bend1_x = x;
                    wire.bend1_y = w[i].start_y;
                    cost = get_cost_vertical(wire, grid, dir);
                    if (cost < minCost) {
                        minWire = wire;
                        minCost = cost;
                    }
                }
            } else {
                for(int x = w[i].start_x + dir.x; w[i].end_x <= x; x += dir.x) {
                    Wire wire = w[i];
                    wire.bend1_x = x;
                    wire.bend1_y = w[i].start_y;
                    cost = get_cost_vertical(wire, grid, dir);
                    if (cost < minCost) {
                        minWire = wire;
                        minCost = cost;
                    }
                }
            }

            //Check Horizontal
            if (dir.y > 0) {
                for(int y = w[i].start_y + dir.y; y <= w[i].end_y; y += dir.y) {
                    Wire wire = w[i];
                    wire.bend1_x = w[i].start_x;
                    wire.bend1_y = y;
                    cost = get_cost_horizontal(wire, grid, dir);
                    if (cost < minCost) {
                        minWire = wire;
                        minCost = cost;
                    }
                }
            } else {
                for(int y = w[i].start_y + dir.y; w[i].end_y <= y; y += dir.y) {
                    Wire wire = w[i];
                    wire.bend1_x = w[i].start_x;
                    wire.bend1_y = y;
                    cost = get_cost_horizontal(wire, grid, dir);
                    if (cost < minCost) {
                        minWire = wire;
                        minCost = cost;
                    }
                }
            }
        }

        //Add new wire to grid and wires
        w[i] = minWire;
        write_wire(w[i], 1, grid);
    }
}

//Route Wires per wire parallelism
void w_route(std::vector<Wire> &w, std::vector<std::vector<int>> &grid, std::vector<std::vector<omp_lock_t>> locks, double SA_prob, Point dim, int t, int it) {
    #pragma omp parallel for
    for(size_t i = 0; i < w.size(); i++) {
        Wire old_wire = w[i]; // Save previous iteration's wire
        if (it > 0) {
            // Remove the previous path
            lock_write_wire(old_wire, -1, grid, locks);
        }

        // Check if current wire is straight
        if (w[i].start_x == w[i].end_x || w[i].start_y == w[i].end_y) {
            // Add straight wire to grid
            lock_write_wire(w[i], 1, grid, locks);
            continue;
        }

        WireCost minCost = {INT_MAX, w[i]};

        if (it > 0 && (double) random() / RAND_MAX < SA_prob) {
            if ((double) random() / RAND_MAX < SA_prob) {
                // Vertical
                w[i].bend1_x = w[i].start_x;
                w[i].bend1_y = randInt(w[i].start_y, w[i].end_y);
            } else {
                // Horizontal
                w[i].bend1_x = randInt(w[i].start_x, w[i].end_x);
                w[i].bend1_y = w[i].start_y;
            }
            lock_write_wire(w[i], 1, grid, locks);
            continue;
        } else {
            #pragma omp declare reduction(min_cost : WireCost : omp_out = (omp_out.c < omp_in.c) ? omp_out : omp_in) initializer(omp_priv={INT_MAX, {0,0,0,0,0,0}})
            Point dir = {(w[i].end_x - w[i].start_x > 0 ? 1 : -1), (w[i].end_y - w[i].start_y > 0 ? 1 : -1)};
            //Check Vertical
            if (dir.x > 0) {
                #pragma omp parallel for reduction(min_cost:minCost)
                for(int x = w[i].start_x + dir.x; x <= w[i].end_x; x += dir.x) {
                    Wire wire = w[i];
                    wire.bend1_x = x;
                    wire.bend1_y = w[i].start_y;
                    WireCost cost = {get_cost_vertical(wire, grid, dir), wire};
                    minCost = (cost.c < minCost.c) ? cost : minCost;
                }
            } else {
                #pragma omp parallel for reduction(min_cost:minCost)
                for(int x = w[i].start_x + dir.x; w[i].end_x <= x; x += dir.x) {
                    Wire wire = w[i];
                    wire.bend1_x = x;
                    wire.bend1_y = w[i].start_y;
                    WireCost cost = {get_cost_vertical(wire, grid, dir), wire};
                    minCost = (cost.c < minCost.c) ? cost : minCost;
                }
            }

            //Check Horizontal
            if (dir.y > 0) {
                #pragma omp parallel for reduction(min_cost:minCost)
                for(int y = w[i].start_y + dir.y; y <= w[i].end_y; y += dir.y) {
                    Wire wire = w[i];
                    wire.bend1_x = w[i].start_x;
                    wire.bend1_y = y;
                    WireCost cost = {get_cost_horizontal(wire, grid, dir), wire};
                    minCost = (cost.c < minCost.c) ? cost : minCost;
                }
            } else {
                #pragma omp parallel for reduction(min_cost:minCost)
                for(int y = w[i].start_y + dir.y; w[i].end_y <= y; y += dir.y) {
                    Wire wire = w[i];
                    wire.bend1_x = w[i].start_x;
                    wire.bend1_y = y;
                    WireCost cost = {get_cost_horizontal(wire, grid, dir), wire};
                    minCost = (cost.c < minCost.c) ? cost : minCost;
                }
            }
        }

        //Add new wire to grid and wires
        w[i] = minCost.w;
        lock_write_wire(w[i], 1, grid, locks);
    }
}

//Route Wires parallel across wires
void a_route(std::vector<Wire> &w, std::vector<std::vector<int>> &grid, std::vector<std::vector<omp_lock_t>> locks, double SA_prob, Point dim, int t, int it) {
    #pragma omp parallel for schedule(dynamic, 2)
    for(size_t i = 0; i < w.size(); i++) {
        Wire old_wire = w[i]; // Save previous iteration's wire
        if (it > 0) {
            // Remove the previous path
            lock_write_wire(old_wire, -1, grid, locks);
        }

        // Check if current wire is straight
        if (w[i].start_x == w[i].end_x || w[i].start_y == w[i].end_y) {
            // Add straight wire to grid
            lock_write_wire(w[i], 1, grid, locks);
            continue;
        }

        int minCost = INT_MAX;
        Wire minWire = w[i];

        if (it > 0 && (double) random() / RAND_MAX < SA_prob) {
            if ((double) random() / RAND_MAX < SA_prob) {
                // Vertical
                w[i].bend1_x = w[i].start_x;
                w[i].bend1_y = randInt(w[i].start_y, w[i].end_y);
            } else {
                // Horizontal
                w[i].bend1_x = randInt(w[i].start_x, w[i].end_x);
                w[i].bend1_y = w[i].start_y;
            }
            lock_write_wire(w[i], 1, grid, locks);
            continue;
        } else {
            Point dir = {(w[i].end_x - w[i].start_x > 0 ? 1 : -1), (w[i].end_y - w[i].start_y > 0 ? 1 : -1)};
            int cost;
            //Check Vertical
            if (dir.x > 0) {
                for(int x = w[i].start_x + dir.x; x <= w[i].end_x; x += dir.x) {
                    Wire wire = w[i];
                    wire.bend1_x = x;
                    wire.bend1_y = w[i].start_y;
                    cost = get_cost_vertical(wire, grid, dir);
                    if (cost < minCost) {
                        minWire = wire;
                        minCost = cost;
                    }
                }
            } else {
                for(int x = w[i].start_x + dir.x; w[i].end_x <= x; x += dir.x) {
                    Wire wire = w[i];
                    wire.bend1_x = x;
                    wire.bend1_y = w[i].start_y;
                    cost = get_cost_vertical(wire, grid, dir);
                    if (cost < minCost) {
                        minWire = wire;
                        minCost = cost;
                    }
                }
            }

            //Check Horizontal
            if (dir.y > 0) {
                for(int y = w[i].start_y + dir.y; y <= w[i].end_y; y += dir.y) {
                    Wire wire = w[i];
                    wire.bend1_x = w[i].start_x;
                    wire.bend1_y = y;
                    cost = get_cost_horizontal(wire, grid, dir);
                    if (cost < minCost) {
                        minWire = wire;
                        minCost = cost;
                    }
                }
            } else {
                for(int y = w[i].start_y + dir.y; w[i].end_y <= y; y += dir.y) {
                    Wire wire = w[i];
                    wire.bend1_x = w[i].start_x;
                    wire.bend1_y = y;
                    cost = get_cost_horizontal(wire, grid, dir);
                    if (cost < minCost) {
                        minWire = wire;
                        minCost = cost;
                    }
                }
            }
        }

        //Add new wire to grid and wires
        w[i] = minWire;
        lock_write_wire(w[i], 1, grid, locks);
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

    std::vector locks((dim_y + GRID_SIZE - 1) / GRID_SIZE,
                      std::vector<omp_lock_t> ((dim_x + GRID_SIZE - 1) / GRID_SIZE));

    /* Initialize any additional data structures needed in the algorithm */
    for (auto &row : locks) {      // Use reference to avoid copying
        for (auto &lock : row) {   // Iterate by reference
            omp_init_lock(&lock);
        }
    }

    omp_set_num_threads(num_threads);

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
            w_route(wires, occupancy, locks, SA_prob, {dim_x, dim_y}, num_threads, i);
    } else /*parallel_mode == 'A'*/ {
        for (int i = 0; i < SA_iters; i++)
            a_route(wires, occupancy, locks, SA_prob, {dim_x, dim_y}, num_threads, i);
    }

    for (auto &row : locks) {
        for (auto &lock : row) {
            omp_destroy_lock(&lock);
        }
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