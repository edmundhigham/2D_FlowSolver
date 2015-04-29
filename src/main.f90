! Main Program
! E. Higham Feb 2014
! Final Year Project - 2D Laminar Flow Solver
program main

#include <src/PreProcVars.f90>

use SetPrecision
use MeshGen
use Derivatives
use SimulationVars
use Poissonmod
use RKMod

implicit none

!=============================================================================!
!                            Flow Variable Declaration                        !
!=============================================================================!
real(wp), allocatable ::   psi(:,:)     ! 2D Stream Funtion Array           
real(wp), allocatable :: omega(:,:)     ! 2D Vorticity Array
real(wp), allocatable ::     u(:,:)     ! 2D Horizontal Velocity Array
real(wp), allocatable ::     v(:,:)     ! 2D Vertical Velocity Array
real(wp)              ::    Re          ! Reynolds Number

! Other Declarations ==========================================================
integer               :: i,j,k
character(len=100)    :: filename
integer               :: nrec

! Timing Variables
real(wp) :: start, finish

! Resume Last Simulation?
integer :: resume

! If resume = 1, then we need some buffer variables
real(wp), allocatable :: buf1(:), buf2(:)

! Grid Pointer to for heap deallocation at end of simulation
type(GRID), pointer :: ptr => null()

!=============================================================================!
!                             Computational Domain                            !
!=============================================================================!

! Limits of Computational Domain
xLim = (/-1.0_wp, 1.0_wp/)
yLim = (/-1.0_wp, 1.0_wp/)

! Build Computational Domain
call Mesh_MultiGrid(xLim,yLim,1,1)

! Allocate Flow Variables
allocate(psi(1:nx,1:ny),omega(1:nx,1:ny),u(1:nx,1:ny),v(1:nx,1:ny))

! PoissonSolver Variables
tol = 0.00000001_wp

!============================================================================!
!                             Set up Time Domain                             !
!============================================================================!

tLim = (/32.0_wp,34.0_wp/)
nt   = 2**11+1

call time(tLim,nt)          ! Compute Time Step

!============================================================================!
!                            Initialise Flow Domain                          !
!============================================================================!

!Re = 625                    ! Clercx et al.
Re = 1000

nrec = wp*(1+1+1+nx+ny+4*nx*ny)             ! Record Length to Write

resume = 1

SELECT CASE (resume)

    CASE (0)
    
        call Initialise(omega)      ! Initialise Vorticity Field
        call CPU_TIME(start)
        call CavityInterp(psi,omega)
        !call MG(psi,omega)          ! Compute Stream Function
        !!psi(:,:) = 0.0_wp
        call psiBC(psi,1)
        call CPU_TIME(finish)
        
        print*, '========================================================='
        print '("Time Elaspsed = ",f6.3," seconds.")',finish-start
        print*, '========================================================='
        
        call psi2u(psi,u,v)         ! Compute Velocity Field (/u,v/)
        
        omega = u2omega(u,v)        ! Recompute Vorticity Distribution from
                                    ! Velocity Field Which Satisfies Continuity
        
        call RK4_CFL(u,v)           ! Check dt Satisfies CFL analysis
        
        ! Save Initial Conditions
        write(filename,"(A16,I5,A4)") 'dat/Iteration',10000,'.bin'
        
        OPEN(UNIT=12, FILE=filename, ACTION="write", STATUS="replace", &
                      FORM='unformatted', ACCESS='direct', RECL= nrec)
         
        WRITE(12,rec=1) t, REAL(nx,wp), REAL(ny,wp), x(:), y(:), &
                           psi(:,:), omega(:,:), u(:,:), v(:,:)
        
        CLOSE(12)
    
        k = 0

    CASE (1)

        allocate(buf1(1), buf2(1))

        ! Load Last File

        k = 38912

        write(filename,"(A16,I5,A4)") 'dat/Iteration',10000+k,'.bin'
        
        OPEN(UNIT=12, FILE=filename, ACTION="read", &
                      FORM='unformatted', ACCESS='direct', RECL= nrec)
         
        READ(12,rec=1) t, buf1, buf2, x(1:nx), y(1:ny), &
                           psi(1:nx,1:ny), omega(1:nx,1:ny), &
                             u(1:nx,1:ny),     v(1:nx,1:ny)
        
        CLOSE(12)

        deallocate(buf1,buf2)

END SELECT



!==========================================================================!
!                          Time Marching                                   !
!==========================================================================!

do while(t<tLim(2))
  t = t + dt
  k = k + 1
  
  print*,'_______________________________________'
  print*,'Iteration',k
  
  call RK4_CFL(u,v)
  
  call RKsolver(omega,u,v,psi,Re)
  
  
  ! Save 100 iterations ===================================================
  
  if (MOD(k,32)==0) then
  
    write(filename,"(A16,I5,A4)") 'dat/Iteration',10000+k,'.bin'
    
    OPEN(UNIT=12, FILE=filename, ACTION="write", STATUS="replace", &
                  FORM='unformatted', ACCESS='direct', RECL= nrec)
     
    WRITE(12,rec=1) t, REAL(nx,wp), REAL(ny,wp), x(:), y(:), &
                       psi(:,:), omega(:,:), u(:,:), v(:,:)
    
    CLOSE(12)
  
  end if

end do


! Tidy Heap
deallocate(psi,omega,u,v)

# if .NOT.defined(MultiGrid)
    do i = 1,1

# elif defined(MultiGrid) && .NOT.defined(nGrids3) && .NOT.defined(nGrids4) && .NOT.defined(nGrids5) && .NOT.defined(nGrids6)
    do i = 1,2
    
# elif defined(MultiGrid) && defined(nGrids3) && .NOT.defined(nGrids4) && .NOT.defined(nGrids5) && .NOT.defined(nGrids6)
    do i = 1,3

# elif defined(MultiGrid) && defined(nGrids3) && defined(nGrids4) && .NOT.defined(nGrids5) && .NOT.defined(nGrids6)
    do i = 1,4

# elif defined(MultiGrid) && defined(nGrids3) && defined(nGrids4) && defined(nGrids5) && .NOT.defined(nGrids6)
    do i = 1,5

# elif defined(MultiGrid) && defined(nGrids3) && defined(nGrids4) && defined(nGrids5) && defined(nGrids6)
    do i = 1,6

# endif

    call GetMesh(i,ptr)

    deallocate(ptr%xi )
    deallocate(ptr%eta)

#ifdef xGridRefinement     
    ! Metrics of Transformation
    deallocate(ptr%DxiDx  )
    deallocate(ptr%DxiDxsq)
    deallocate(ptr%D2xiDx2)
#endif
#ifdef yGridRefinement
    ! Metrics of Transformation
    deallocate(ptr%DetaDy  )
    deallocate(ptr%DetaDysq)
    deallocate(ptr%D2etaDy2)
#endif

    ptr => null()
    
  end do

contains

subroutine CavityInterp(psi,omega)
implicit none

real(wp) ::   psi(:,:)
real(wp) :: omega(:,:)

real(wp), allocatable ::   psiTemp1(:,:)
real(wp), allocatable :: omegaTemp1(:,:)

real(wp), allocatable ::   psiTemp2(:,:)
real(wp), allocatable :: omegaTemp2(:,:)

! Reading buffers
real(wp), allocatable ::  bufx(:)
real(wp), allocatable ::  bufy(:)
real(wp), allocatable ::  bufu(:,:)
real(wp), allocatable ::  bufv(:,:)
real(wp), allocatable ::  buft(:)
real(wp), allocatable :: bufnx(:)
real(wp), allocatable :: bufny(:)

type(GRID), pointer :: grid_ptr =>null()

integer :: count

call GetMesh(4,grid_ptr)

allocate(bufx(1:grid_ptr%nxi),bufy(1:grid_ptr%neta),buft(1),bufnx(1),bufny(1))

allocate(    bufu(1:grid_ptr%nxi,1:grid_ptr%neta))
allocate(    bufv(1:grid_ptr%nxi,1:grid_ptr%neta))
allocate(omegaTemp1(1:grid_ptr%nxi,1:grid_ptr%neta))
allocate(psiTemp1(1:grid_ptr%nxi,1:grid_ptr%neta))

! Read in Old Simulation

count = wp*(1+1+1+grid_ptr%nxi+grid_ptr%neta+4*grid_ptr%nxi*grid_ptr%neta)

write(filename,"(A15)") 'dat/Initial.bin'


OPEN(UNIT=12, FILE=filename, ACTION="read", &
              FORM='unformatted', ACCESS='direct', RECL= count)
 
READ(12,rec=1) buft, bufnx, bufny, bufx(1:grid_ptr%nxi), bufy(1:grid_ptr%neta), &
                   psiTemp1(1:grid_ptr%nxi,1:grid_ptr%neta), omegaTemp1(1:grid_ptr%nxi,1:grid_ptr%neta), &
                       bufu(1:grid_ptr%nxi,1:grid_ptr%neta),       bufv(1:grid_ptr%nxi,1:grid_ptr%neta)


CLOSE(12)

print*, 'nx = ', bufnx
print*, 'ny = ', bufny



deallocate (bufx,bufy,bufu,bufv,buft,bufnx,bufny)

do count = 3,2,-1

call GetMesh(count,grid_ptr)

allocate(  psiTemp2(1:grid_ptr%nxi,1:grid_ptr%neta))
allocate(omegaTemp2(1:grid_ptr%nxi,1:grid_ptr%neta))

  psiTemp2 = Interp2(count+1,  psiTemp1)
omegaTemp2 = Interp2(count+1,omegaTemp1)

deallocate(psiTemp1,omegaTemp1)

allocate(  psiTemp1(1:grid_ptr%nxi,1:grid_ptr%neta))
allocate(omegaTemp1(1:grid_ptr%nxi,1:grid_ptr%neta))


  psiTemp1 = psiTemp2
omegaTemp1 = omegaTemp2

deallocate(omegaTemp2,psiTemp2)

end do


  psi = Interp2(2,psiTemp1  )
omega = Interp2(2,omegaTemp1)

grid_ptr => null()

deallocate(psiTemp1,omegaTemp1)

RETURN
end subroutine CavityInterp


end program main