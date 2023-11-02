import subprocess

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

    for params in params_list:
        print(f'Running with params: {params}')
        run_solver(*params)

if __name__ == '__main__':
    main()
