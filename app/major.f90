  module mpi_mod
  use mpi
  implicit none
  integer:: rank=0, nproc=8
  integer previous, next
  integer ierr
  integer, allocatable:: work(:) !allocate taskes

  contains

  !initialize MPI
  subroutine initMPI()
  implicit none
  call MPI_INIT( ierr )
  call MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, nproc, ierr )
  allocate(work(nproc+1))
  previous = merge(MPI_PROC_NULL, rank - 1, rank==0)
  next     = merge(MPI_PROC_NULL, rank + 1, rank==nproc-1)
  end subroutine

  !allocate to processor
  pure subroutine distributeWork(work, nproc, nwork)
  ! [work(i-1), work(i)-1]
  implicit none
  integer, intent(in):: nproc, nwork
  integer, intent(out):: work(:)
  integer i, m, n
  work(1) = 2 !from row 2
  work(2:) = nwork / nproc
  i = mod(nwork,nproc) + 1
  work(2:i) = work(2:i) + 1
  do i = 2, nproc+1
    work(i) = work(i) + work(i-1)
  end do
  end subroutine

  subroutine exchangeBoundaryData(nx, ny, d)
  implicit none
  integer, intent(in):: nx, ny
  real(8), intent(inout):: d(ny,nx)
  integer status(MPI_STATUS_SIZE)
  ! send and receive boundary data
  call MPI_Sendrecv(d(:,work(rank+1)), ny, MPI_REAL8, previous, 0, d(:,work(rank+2)), ny, MPI_REAL8, next, 0, MPI_COMM_WORLD, status, ierr)
  call MPI_Sendrecv(d(:,work(rank+2)-1), ny, MPI_REAL8, next, 1, d(:,work(rank+1)-1), ny, MPI_REAL8, previous, 1, MPI_COMM_WORLD, status, ierr)
  end subroutine

  end module

  program Test
  use mpi_mod
  implicit none
  integer:: nx, ny
  real(8):: lx, ly, dx, dy
  real(8):: Ut, Ub, Vl, Vr
  real(8):: t, dt, tol, beta, nu, tend, err, rhs
  real(8), allocatable:: psi0(:,:), omega0(:,:), psi(:,:), omega(:,:)
  real(8), allocatable:: u(:,:), v(:,:)
  integer maxIt, it, i, j
  ! init MPI
  call initMPI()
  nx = 21
  ny = 21
  lx = 1.0d0
  ly = 1.0d0
  dx = lx / (nx-1)
  dy = ly / (ny-1)
  Ut = 5.0d0 ! u top wall
  Ub = 0     ! u bottom wall
  Vl = 0     ! V left wall
  Vr = 0     ! V right wall
  ! memery allocate
  allocate(psi0(ny,nx), omega0(ny,nx), psi(ny,nx), omega(ny,nx))
  psi0 = 0
  omega0 = 0
  ! Apply boundary conditions
  omega0(1,2:nx-1)  = -2.0d0 * psi0(2,2:nx-1)/dy/dy + 2.0d0*Ub/dy; ! vorticity on bottom wall (stationary)
  omega0(ny,2:nx-1) = -2.0d0 * psi0(ny-1,2:nx-1)/dy/dy - 2.0d0*Ut/dy; ! vorticity on top wall (moving at Uwall)
  omega0(2:ny-1,nx) = -2.0d0 * psi0(2:ny-1,nx-1)/dx/dx; ! right wall
  omega0(2:ny-1,1)  = -2.0d0 * psi0(2:ny-1,2)/dx/dx; ! left wall
  psi = psi0
  omega = omega0

  t = 0
  maxIt = 300
  beta = 1.5d0
  nu = 0.05d0
  tol = 1.0d-3
  dt = min(0.25d0*dx*dx/nu, 4.0d0*nu/Ut/Ut)
  tend = 4 * dt
  if(rank==0) then
    write(*,"('dt = ', f7.4)") dt
    write(*,"('Reynods = ', f9.4)") Ut*lx/nu
  end if
  dt = 0.008d0

  call distributeWork(work, nproc, nx-2)

  do while(t<tend)
    ! first find psi by solving the Poisson equation
    do it = 1, maxIt
      err = 0
      do i = work(rank+1), work(rank+2)-1
        do j = 2, ny-1
          psi0(j,i) = psi(j,i)
          rhs = (dx*dy)**2 * omega(j,i) + dx*dx * (psi(j+1,i)+psi(j-1,i)) + dy*dy * (psi(j,i+1)+psi(j,i-1))
          rhs = beta * rhs / 2.0d0 / (dx*dx+dy*dy)
          psi(j,i) = rhs + (1-beta)*psi(j,i)
          err = err + (psi(j,i)-psi0(j,i))**2
        end do
      end do
      ! send message
      call exchangeBoundaryData(nx, ny, psi)
      call MPI_Allreduce(err, rhs, 1, MPI_REAL8, MPI_SUM, MPI_Comm_World, ierr)
      err = sqrt( rhs )
      if(err<tol) exit
    end do

    ! update omega
    omega0(:,:) = omega(:,:)
    do i = work(rank+1), work(rank+2)-1
      do j = 2, ny-1
        rhs = (omega0(j,i+1)-2.0d0*omega0(j,i)+omega0(j,i-1)) / (dx*dx) + &
          (omega0(j+1,i)-2.0d0*omega0(j,i)+omega0(j-1,i)) / (dy*dy)  ! Dxy
        rhs = rhs * nu
        rhs  = rhs - (psi(j+1,i)-psi(j-1,i))/2.0d0/dy * (omega0(j,i+1)-omega0(j,i-1))/2.0d0/dx
        rhs  = rhs + (omega0(j+1,i)-omega0(j-1,i))/2.0d0/dy * (psi(j,i+1)-psi(j,i-1))/2.0d0/dx
        omega(j,i) = omega0(j,i) + dt*rhs
      end do
    end do
    ! send message
    call exchangeBoundaryData(nx, ny, omega)

    ! apply boundary conditions
    omega(1,2:nx-1)  = -2.0d0*psi(2,2:nx-1)/dy/dy + 2.0d0*Ub/dy ! vorticity on bottom wall (stationary)
    omega(ny,2:nx-1) = -2.0d0*psi(ny-1,2:nx-1)/dy/dy - 2.0d0*Ut/dy ! vorticity on top wall (moving at Uwall)
    omega(2:ny-1,nx) = -2.0d0*psi(2:ny-1,nx-1)/dx/dx ! right wall
    omega(2:ny-1,1)  = -2.0d0*psi(2:ny-1,2)/dx/dx ! left wall

    t = t + dt
  end do

  psi(:,:) = -psi(:,:)

  !allocate(u(ny,nx), v(ny,nx))
  !u = 0
  !v = 0
  !u(2:ny-1, 2:nx-1) = (psi(3:ny, 2:nx-1)-psi(1:ny-2, 2:nx-1)) / 2.0d0 / dy
  !v(2:ny-1, 2:nx-1) = -(psi(2:ny-1, 3:nx) - psi(2:ny-1, 1:nx-2)) / 2.0d0 / dx

  work(1) = 1; work(nproc+1) = nx+1
  do i = 0, nproc-1
    if(i==rank) then
      if(i==0) then
        open(10, file='res.csv')
      else
        open(10, file='res.csv',position='append')
      end if
      do j = work(rank+1), work(rank+2)-1
        write(10,"(*(g0,','))") psi(:,j)
      end do
      close(10)
    end if
    call MPI_Barrier(MPI_Comm_World, ierr)
  end do

  call MPI_FINALIZE(ierr)
  end program
