 !> Extending the work on the serial version, parallelism is implemented with openMPI,
 !> then the cartesian topology is implemented.
 !> GPT-4 was employed as a debugging technique. 
 !----------------------------------------------------------------------------
  module mpi_mod
  implicit none
  include 'mpif.h'
  integer:: rank=0, nproc, npx, npy
  integer:: COMM_CART, xyrank(2) = 0
  integer:: nx, ny
  integer:: UR = -1, DR = -1, LR = -1, RR = -1
  integer ierr
  integer ix1, ix2, iy1, iy2
  character(len=20) :: arg
  integer, allocatable:: xwork(:) !allocate tasks
  integer, allocatable:: ywork(:) 
  real(8), allocatable:: tempDat(:)
  contains
!--------------------------------------------------------------------
! subroutine init() included other two subroutine including allocating jobs to different procs
! and boundary dara exchange.
  !initialize mpi
  subroutine init()
  implicit none
  call MPI_INIT( ierr )

!------------------------------------------------------
!> This part of the code is using a get_command_argument command
!> to read the input parameters (npx, npy, nx and ny) and cast them to their respective types
!------------------------------------------------------
  ! read input parameters
  if (command_argument_count() /= 4) then
    if (rank .eq. 0) then
      print *, 'Usage: mpirun -np <num_procs> ./solver_parallel npx npy nx ny'
    end if
    call MPI_FINALIZE(ierr)
    stop 1
  end if
  call get_command_argument(1, arg)
  read(arg, *) npx
  call get_command_argument(2, arg)
  read(arg, *) npy
  call get_command_argument(3, arg)
  read(arg, *) nx
  call get_command_argument(4, arg)
  read(arg, *) ny 
  call MPI_COMM_SIZE( MPI_COMM_WORLD, nproc, ierr )
  call MPI_Comm_RANK( COMM_CART, rank, ierr ) 
  ! Panic if number of processors does not match npx and npy
  if (nproc /= npx * npy) then
    if (rank .eq. 0) then
      print *, "ERROR! npx*npy does not equal to proc"
    end if
    call MPI_FINALIZE(ierr)
    stop 1
  end if
!-------------------------------------------------------
  ! allocate work domain
! allocate xwork as (npx+1) and ywork as (npy+1) as these are the boundary number of each domain 
  allocate(xwork(npx+1), ywork(npy+1))
! use subroutine to allocate task--[xwork (nx-2, boundary data not considered)] to different processor
  call distributeWork(xwork, npx, nx-2)
  call distributeWork(ywork, npy, ny-2)
  
  
! To improve efficiency, we use temporary data storage to save the data in 
! rows (discontinuous data) and then transform it into columns (continuous data)."
  allocate(tempDat(nx)) !tempro data
  
! use cartesian topology to decomposet the domain 
  ! new doamin 
  call MPI_CART_CREATE(MPI_COMM_WORLD, 2, [npy,npx], [.false., .false.], .true., COMM_CART, ierr)
  call MPI_CART_COORDS( COMM_CART, rank, 2, xyrank, ierr )
! change to a 1 start rank
  xyrank = xyrank + 1
  ! doamin for different processor
  ix1 = xwork(xyrank(2))
  ix2 = xwork(xyrank(2)+1) - 1
  iy1 = ywork(xyrank(1))
  iy2 = ywork(xyrank(1)+1) - 1
  ! get rank
  call MPI_CART_SHIFT(COMM_CART, 0, 1, DR, UR, ierr) ! down and up
  call MPI_CART_SHIFT(COMM_CART, 1, 1, LR, RR, ierr) ! left and right
  end subroutine
!-------------------------------------------------------------------------------


! distributeWork subdomain----divid the nwork (nodes-2) by nproc and showing the
! boundary of each sub domain by work array.
  !allocate task to different processor
  pure subroutine distributeWork(work, nproc, nwork)
  ! [work(i-1), work(i)-1]
  implicit none
  integer, intent(in):: nproc, nwork
  integer, intent(out):: work(:)
  integer i, m, n
  
! work start at second node, and use nwork/ nproc determin the nodes in each proc
! if mod().NE. zero, the rest of the code will put the remainding node evenly in to the 
! first several number in work array.
! if mod() .EQ. zero, work(2:1) is invalid and will have no effect on work array.
! finally, the upper limit and lower limit in each processor is determained by add 
! work array together
  work(1) = 2 !row two start
  work(2:) = nwork / nproc
  i = mod(nwork,nproc) + 1
  work(2:i) = work(2:i) + 1
  do i = 2, nproc+1
    work(i) = work(i) + work(i-1)
  end do
  end subroutine
!---------------------------------------------------------------------------


! exchanging the data on each subdomain boundary to compute the data on boundary 
  subroutine exchangeBoundaryData(nx, ny, d)
  
  implicit none
  integer, intent(in):: nx, ny
  real(8), intent(inout):: d(ny,nx)
  integer status(MPI_STATUS_SIZE)
  ! send and receive boundary data
! the tempdat is used here to save the data in row direction 
 ! U direction
  if(UR>=0) then
    tempDat(ix1:ix2) = d(iy2,ix1:ix2)
    call MPI_Send(tempDat(ix1:ix2), ix2-ix1+1, MPI_REAL8, UR, 1, COMM_CART, ierr)
  end if
  if(DR>=0) then
    call MPI_Recv(tempDat(ix1:ix2), ix2-ix1+1, MPI_REAL8, DR, 1, COMM_CART, status, ierr)
    d(iy1-1,ix1:ix2) = tempDat(ix1:ix2)
  end if

  ! D direction
  if(DR>=0) then
    tempDat(ix1:ix2) = d(iy1,ix1:ix2)
    call MPI_Send(tempDat(ix1:ix2), ix2-ix1+1, MPI_REAL8, DR, 2, COMM_CART, ierr)
  end if
  if(UR>=0) then
    call MPI_Recv(tempDat(ix1:ix2), ix2-ix1+1, MPI_REAL8, UR, 2, COMM_CART, status, ierr)
    d(iy2+1,ix1:ix2) = tempDat(ix1:ix2)
  end if
! as data in left and right dirction is in column (continuous data) so don't need tempdat
  ! L direction
  if(LR>=0) call MPI_Send(d(iy1:iy2,ix1), iy2-iy1+1, MPI_REAL8, LR, 3, COMM_CART, ierr)
  if(RR>=0) call MPI_Recv(d(iy1:iy2,ix2+1), iy2-iy1+1, MPI_REAL8, RR, 3, COMM_CART, status, ierr)

  ! R direction
  if(RR>=0) call MPI_Send(d(iy1:iy2,ix2), iy2-iy1+1, MPI_REAL8, RR, 4, COMM_CART, ierr)
  if(LR>=0) call MPI_Recv(d(iy1:iy2,ix1-1), iy2-iy1+1, MPI_REAL8, LR, 4, COMM_CART, status, ierr)
  end subroutine

  end module
!--------------------------------------------------------------------------------



! the main coding part, the discretised equation and main loop structure is using the 
! Vorticity-Streamfunction FD Python written by Prof. Tony Saad
! (https://nbviewer.org/github/saadtony/uCFD/blob/master/Vorticity-Streamfunction%20FDM%20-%20lid-driven%20cavity.ipynb)
  program solver_parallel
  use mpi_mod
  implicit none
  real(8):: lx, ly, dx, dy
  real(8):: Ut, Ub, Vl, Vr
  real(8):: t, dt, tol, beta, nu, tend, err, rhs
  character(len=256) :: filename
  real(8), allocatable:: psi0(:,:), omega0(:,:), psi(:,:), omega(:,:)
  real(8), allocatable:: u(:,:), v(:,:)
  real(8) time1, time2
  integer maxIt, it, i, j, k
  ! init MPI
  call init()

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
  psi = 0
  omega = 0
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
! dt has been set to a upper limit value of stable result
  dt = 0.25d0*dx*dx/nu
  !dt = min(0.25d0*dx*dx/nu, 4.0d0*nu/Ut/Ut)
  tend = 400 * dt
  if(rank==0) then
    write(*,"('dt = ', f7.4)") dt
    write(*,"('Reynods = ', f9.4)") Ut*lx/nu
  end if


  time1 = mpi_wtime()
  do while(t<tend)
    ! first find psi by solving the Poisson equation
    do it = 1, maxIt
      err = 0
      do i = ix1, ix2
        do j = iy1, iy2
          psi0(j,i) = psi(j,i)
          rhs = (dx*dy)**2 * omega(j,i) + dx*dx * (psi(j+1,i)+psi(j-1,i)) + dy*dy * (psi(j,i+1)+psi(j,i-1))
          rhs = beta * rhs / 2.0d0 / (dx*dx+dy*dy)
          psi(j,i) = rhs + (1-beta)*psi(j,i)
          err = err + (psi(j,i)-psi0(j,i))**2
        end do
      end do
      ! send message
      rhs = err
! computing global error by Allreduce 
      call MPI_Allreduce(rhs, err, 1, MPI_REAL8, MPI_SUM, MPI_Comm_World, ierr)
      err = sqrt( err )
      if(err<tol) exit
      call exchangeBoundaryData(nx, ny, psi)
    end do

    ! update omega
    omega0(:,:) = omega(:,:)
    do i = ix1, ix2
      do j = iy1, iy2
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
  time2 = mpi_wtime()

  allocate(u(ny,nx), v(ny,nx))
  u = 0
  v = 0
  u(ny, 2:nx-1) = Ut
  u(2:ny-1, 2:nx-1) = (psi(3:ny, 2:nx-1)-psi(1:ny-2, 2:nx-1)) / 2.0d0 / dy
  v(2:ny-1, 2:nx-1) = -(psi(2:ny-1, 3:nx) - psi(2:ny-1, 1:nx-2)) / 2.0d0 / dx
!-------------------------------------------------------------------------




! based on the http://fcode.cn/guide-86-1.html (Fortran流文件-读写二进制任意位置), the data 
! is initially put into a stream file with binary data, and then conver to a CSV file.
  ! output data
  write(filename, '(A,I0,A,I0,A,I0,A,I0,A)') 'psi_', npx, '-', npy, '_', nx, '-', ny, '.csv'
! add boundary data to the main array
  if(xyrank(1)==1) iy1 = 1
  if(xyrank(1)==npy) iy2 = ny
  if(xyrank(2)==1) ix1 = 1
  if(xyrank(2)==npx) ix2 = nx
  do i = 0, nproc-1
    if(i==rank) then
      open(10, file=filename, access='stream')
      do j = ix1, ix2
        do k = iy1, iy2
          write(10,pos=((j-1)*ny+k-1)*8+1) psi(k,j)
        end do
      end do
      close(10)
    end if
    call MPI_Barrier(COMM_CART, ierr)
  end do
  if(rank==0) then
    open(10, file=filename, access='stream')
    read(10) psi(:,:)
    close(10,status='delete')
    open(10, file=filename)
    do j = 1, ny
      write(10,"(*(g0,','))") psi(j,:)
    end do
    close(10)
  end if

  if(rank==0) write(*,"(a,f10.6)") 'wtime = ', time2-time1
  call MPI_FINALIZE(ierr)
  end program solver_parallel
