program main
    implicit none

    integer, parameter :: nx = 21, ny = 21, iter_max = 30
    real(8), parameter :: lx = 1.0d0, ly = 1.0d0, Ut = 5.0d0, Ub = 0.0d0, Vl = 0.0d0, Vr = 0.0d0, beta = 1.5, nu = 0.05
    real(8) :: dx, dy, dt, rhs, tol, t, t_end, err
    real(8), allocatable :: psi(:,:), omega(:,:), psi_old(:,:), C_x(:,:), C_y(:,:), D_xy(:,:)
    integer :: i, j, iter_num

    ! Initialize parameters and arrays
    dx = lx / (nx - 1)
    dy = ly / (ny - 1)

    t = 0.0d0
    dt = 0.008
    t_end = 4 * min(0.25*dx*dx/nu, 4.0*nu/Ut/Ut)

    tol = 1e-3

    ! Print parameters
    print *, "dt = ", dt
    print *, "t_end = ", t_end
    print *, "Re = ", Ut * lx / nu

    allocate(psi(nx, ny), omega(nx, ny), psi_old(nx, ny))
    psi = 0.0d0
    omega = 0.0d0


    ! Apply initial boundary conditions
    ! Vorticity on bottom wall (stationary)
    omega(2:nx-1, 1) = -2.0d0 * psi(2:nx-1, 2) / (dy * dy) + 2.0d0 * Ub / dy

    ! Vorticity on top wall (moving at Uwall)
    omega(2:nx-1, ny) = -2.0d0 * psi(2:nx-1, ny-1) / (dy * dy) - 2.0d0 * Ut / dy

    ! Vorticity on right wall
    omega(nx, 2:ny-1) = -2.0d0 * psi(nx-1, 2:ny-1) / (dx * dx)

    ! Vorticity on left wall
    omega(1, 2:ny-1) = -2.0d0 * psi(2, 2:ny-1) / (dx * dx)

    do while (t < t_end)
        ! first find psi by solving the Poisson equation
        iter_num = 0
        err = 1e3
        do while (err > tol)
            if (iter_num > iter_max) then
                print *, "Exceed max iteration: ", iter_max, " at t =", t, " err = ", err
                exit
            end if

            psi_old = psi
            do i = 2, nx-1
                do j = 2, ny-1
                    rhs = (dx*dy)**2 * omega(i, j) + dx**2 * (psi(i, j+1) + psi(i, j-1)) + dy**2 * (psi(i+1, j) + psi(i-1, j))
                    rhs = beta * rhs / (2.0 * (dx**2 + dy**2))
                    psi(i, j) = rhs + (1 - beta) * psi(i, j)
                end do
            end do
            err = norm_frobenius(psi - psi_old, nx, ny)
            iter_num = iter_num + 1
        end do

        ! update omega
        C_x = -(psi(2:nx-1, 3:ny) - psi(2:nx-1,1:ny-2)) / 2.0 / dy * (omega(3:nx,2:ny-1) - omega(1:nx-2,2:ny-1)) / 2.0 / dx

        C_y  =  (omega(2:nx-1, 3:ny) - omega(2:nx-1,1:ny-2)) / 2.0 / dy * (psi(3:nx,2:ny-1) - psi(1:nx-2,2:ny-1)) / 2.0 / dx

        D_xy =  (omega(3:nx, 2:ny-1) - 2.0*omega(2:nx-1,2:ny-1) + omega(1:nx-2,2:ny-1)) / dx / dx + (omega(2:nx-1,3:ny) &
        -2.0 * omega(2:nx-1,2:ny-1) + omega(2:nx-1,1:ny-2)) / dy / dy

        omega(2:nx-1, 2:ny-1) = omega(2:nx-1 ,2:ny-1) + dt * (C_x + C_y + nu * D_xy)

        ! apply boundary conditions
        omega(2:nx-1, 1) = -2.0 * psi(2:nx-1, 2) / dy / dy + 2.0*Ub/dy; ! vorticity on bottom wall (stationary)
        omega(2:nx-1, ny) = -2.0 * psi(2:nx-1, ny-1) / dy / dy - 2.0*Ut/dy; ! vorticity on top wall (moving at Uwall)
        omega(nx, 2:ny-1) = -2.0 * psi(nx-1, 2:ny-1) / dx / dx; ! right wall
        omega(1, 2:ny-1) = -2.0 * psi(2, 2:ny-1) / dx / dx; ! left wall

        ! update time
        t = t + dt
    end do

    ! Save results to files
    open(unit=10, file='psi.dat', status='replace')
    open(unit=11, file='omega.dat', status='replace')

    do j = 1, ny
        do i = 1, nx
            write(10, '(E20.10)') psi(i, j)
            write(11, '(E20.10)') omega(i, j)
        end do
        write(10, *)
        write(11, *)
    end do

    close(10)
    close(11)

    deallocate(psi, omega, psi_old)

    contains
    function norm_frobenius(A, n, m) result(norm)
        implicit none
        real(8), intent(in) :: A(n, m)
        integer, intent(in) :: m, n
        real(8) :: norm
        integer :: x, y

        norm = 0.0d0
        do x = 1, m
            do y = 1, n
                norm = norm + A(x, y)**2
            end do
        end do
        norm = sqrt(norm)
    end function norm_frobenius    
end program main