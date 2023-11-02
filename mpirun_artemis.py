import os
import shutil
import subprocess

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
        # tests for domain size
        (4, 7, 50, 50),
        (4, 8, 50, 50),
        (4, 9, 50, 50),
        (4, 10, 50, 50),
    ]

    for params in params_list:
        print(f'Running with params: {params}')
        run_solver(*params)
        move_csv_to_outputs()


if __name__ == '__main__':
    main()
