import os
import shutil
import subprocess
import concurrent.futures

def move_csv_to_outputs():
    # Ensure the 'outputs' directory exists
    if not os.path.exists('outputs'):
        os.makedirs('outputs')

    # Iterate over all files in the current directory
    for file in os.listdir('.'):
        if file.endswith('.csv'):
            shutil.move(file, os.path.join('outputs', file))

def run_solver(npx, npy, nx, ny):
    num_processors = npx * npy

    args = [str(arg) for arg in (npx, npy, nx, ny)]

    binary = 'build/gfortran_6FB96CAE7088C0B9/app/solver_parallel'
    command = ['mpirun', '-np', str(num_processors), binary] + args

    process = subprocess.run(command, check=True)

def main():
    params_list = [
        # tests for grid size
        (3, 4, 25, 25),
        (3, 4, 50, 50),
        (3, 4, 75, 75),
        (3, 4, 100, 100),
        # tests for domain size
        (4, 3, 50, 50),
        (4, 4, 50, 50),
        (4, 5, 50, 50),
        (4, 6, 50, 50),
        (4, 7, 50, 50),
        (4, 8, 50, 50),
        # tests for domain topology
        (1, 24, 50, 50),
        (2, 12, 50, 50),
        (3, 8, 50, 50),
        (4, 6, 50, 50),
        (6, 4, 50, 50),
        (8, 3, 50, 50),
        (12, 2, 50, 50),
        (24, 1, 50, 50),
    ]

    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = [executor.submit(run_solver, *params) for params in params_list]
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as e:
                print(f"An exception occurred: {e}")

    move_csv_to_outputs()

if __name__ == '__main__':
    main()