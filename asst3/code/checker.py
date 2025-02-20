import subprocess
import os
import re
import argparse

def run_command(command):
    """Run a shell command and return the output"""
    try:
        result = subprocess.run(command, shell=True, check=True, 
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              text=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {command}")
        print(f"Error message: {e.stderr}")
        return None

def parse_output(output):
    """Parse the wireroute output for relevant metrics"""
    if not output:
        return None
    
    metrics = {}
    
    # Extract initialization time
    init_match = re.search(r'Initialization time \(sec\): ([0-9.]+)', output)
    metrics['init_time'] = float(init_match.group(1)) if init_match else 0
    
    # Extract computation time
    comp_match = re.search(r'Computation time \(sec\): ([0-9.]+)', output)
    metrics['comp_time'] = float(comp_match.group(1)) if comp_match else 0
    
    # Calculate total time
    metrics['total_time'] = metrics['init_time'] + metrics['comp_time']
    
    # Extract total cost
    cost_match = re.search(r'Total cost: ([0-9]+)', output)
    metrics['total_cost'] = int(cost_match.group(1)) if cost_match else 0
    
    # Extract max occupancy
    occupancy_match = re.search(r'Max occupancy: ([0-9]+)', output)
    metrics['max_occupancy'] = int(occupancy_match.group(1)) if occupancy_match else 0
    
    return metrics

def main():
    # Set up command line argument parser
    parser = argparse.ArgumentParser(description='Run wire routing evaluation')
    parser.add_argument('-f', '--input_file', required=True, help='Input file path')
    parser.add_argument('-n', '--num_threads', type=int, required=True, help='Number of threads')
    parser.add_argument('-m', '--mode', choices=['W', 'A'], required=True, 
                        help='Mode: W for within, A for across')
    parser.add_argument('-b', '--batch_size', type=int, 
                        help='Batch size (default: 4 * num_threads)')
    parser.add_argument('-p', '--probability', type=float, default=0.1,
                        help='Simulated annealing probability (default: 0.1)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Verbose (default: False)')
    
    args = parser.parse_args()
    
    # Set default batch size if not provided
    if args.batch_size is None:
        args.batch_size = 4 * args.num_threads
    
    # Verify input file exists
    if not os.path.exists(args.input_file):
        print(f"Error: Input file {args.input_file} does not exist!")
        return
    
    # Clean and make
    print("Cleaning and building...")
    run_command("make clean")
    make_output = run_command("make")
    if make_output is None:
        print("Compilation failed!")
        return
    
    if args.verbose:
        print(make_output)

    # Run wireroute
    base_name = os.path.basename(args.input_file)
    print(f"\nProcessing {base_name}")
    print(f"Mode: {'Within' if args.mode == 'W' else 'Across'}")
    print(f"Threads: {args.num_threads}")
    print(f"Batch size: {args.batch_size}")
    print(f"Probability: {args.probability}")
    
    # Construct and run wireroute command
    cmd = (f"./wireroute -f {args.input_file} -n {args.num_threads} "
           f"-m {args.mode} -b {args.batch_size} -p {args.probability}")
    output = run_command(cmd)
    if args.verbose:
        print(output)
    
    if output:
        # Run validation
        input_base = os.path.splitext(args.input_file)[0]
        validate_cmd = (f"python3 validate.py "
                       f"-r {input_base}_wires_{args.num_threads}.txt "
                       f"-c {input_base}_occupancy_{args.num_threads}.txt "
                       f"{'-v' if args.verbose else ''}")
        validate_output = run_command(validate_cmd)
        
        if validate_output and "Validate succeeded" in validate_output:
            print("\nValidation successful!")
            metrics = parse_output(output)
            if metrics:
                print(f"\nResults:")
                print(f"  Total time: {metrics['total_time']:.3f} sec")
                print(f"  Initialization time: {metrics['init_time']:.3f} sec")
                print(f"  Computation time: {metrics['comp_time']:.3f} sec")
                print(f"  Total cost: {metrics['total_cost']}")
                print(f"  Max occupancy: {metrics['max_occupancy']}")
        else:
            print("\nValidation failed!")
            print(validate_output)

if __name__ == "__main__":
    main()