import subprocess
import concurrent.futures

def run_solver(npx, npy, nx, ny):
    num_processors = npx * npy

    args = [str(arg) for arg in (npx, npy, nx, ny)]

    binary = 'build/solver_parallel'
    command = ['mpirun', '-np', str(num_processors), binary] + args

    process = subprocess.run(command, check=True)

def main():
    params_list = [
        (2, 2, 21, 21),
        (3, 2, 21, 21),
    ]

    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = [executor.submit(run_solver, *params) for params in params_list]
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as e:
                print(f"An exception occurred: {e}")

if __name__ == '__main__':
    main()
