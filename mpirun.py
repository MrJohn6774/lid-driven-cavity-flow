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

    binary = 'build/solver_parallel'
    command = ['mpirun', '-np', str(num_processors), binary] + args

    process = subprocess.run(command, check=True)

def main():
    params_list = [
        (3, 2, 100, 100),
        (3, 2, 200, 200),       
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
