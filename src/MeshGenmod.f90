! Module for Defining 2D Computational Mesh
! E. Higham Feb 2014
! Final Year Project - 2D Laminar Flow Solver
module MeshGen
! Declare the Precison to Be Used
#include <src/PreProcVars.f90>
use SetPrecision

implicit none
!===========================================================================================!
!                               Mesh Variable Declaration                                   !   
!===========================================================================================!

! Co-ordinate System and Assignment Convention

!       y
!       ^               For some arbitrary array u:
!       |                       u(x,y)
!       |
!       |
!       *------> x

!******************************************************************************************!
! Physical Domain

! Limits of Physical Domain
real(wp) :: xLim(1:2)                                   ! Min and Max Limits in x-Direction
real(wp) :: yLim(1:2)                                   ! Min and Max Limits in y-Direction

! The number of points used in any axis will be 2**n+1, where n is a natural number
! and Must Be Greater Than or Equal to 6 (Limit for 4th Order Derivatives)
integer, parameter :: m = 9
integer, parameter :: n = 9

! Phsyical Domain Properties
integer, parameter :: nx = 2**m+1                       ! Number of Points in x-Direction 
integer, parameter :: ny = 2**n+1                       ! Number of Points in y-Direction

! Axis Co-ordinates in Physical Domain
real(wp) :: x(1:nx)                                     ! Vector of x ordinates
real(wp) :: y(1:ny)                                     ! Vector of y-ordinates

!*******************************************************************************************!
! Computational Domain

! Mapping Co-ordinate System From Physical to Computational Space

!       Physical                          Computational
!
!       y                                 eta
!       ^             ==> eta(y) ==>      ^
!       |                                 |
!       |                                 |
!       |             ==>  xi(x) ==>      |
!       *------> x                        *------> xi


! Decompose Physical Mesh into 4 Uniform Computational Meshes of Differing Resolution:

! grid1 - 2**m    +1   x 2**n    +1    Points in xi- and eta- Respectively  |
!                                                                           |
! grid2 - 2**(m-1)+1   x 2**(n-1)+1    Points in xi- and eta- Respectively  | Coarsening
!                                                                           | 
! grid3 - 2**(m-2)+1   x 2**(n-2)+1    Points in xi- and eta- Respectively  |
!   .         .              .                                              |
!   .         .              .                                             \ /
!   .         .              .                                              V
!                                                                           
! gridl - 2**(m-l-1)+1 x 2**(n-l-1)+1  Points in xi- and eta- Respectively


! This is Done to Increase the Rate of Convergence of the Poisson Solver

! All 3 Meshes Have the Same Limits:
real(wp) ::  xiLim(1:2)                                    ! Min and Max Limits in  xi-Direction
real(wp) :: etaLim(1:2)                                    ! Min and Max Limits in eta-Direction

type GRID
    SEQUENCE

    ! Type Structure for Each Mesh
    integer :: nxi                                         ! Number of Points in  xi-Direction 
    integer :: neta                                        ! Number of Points in eta-Direction

    ! Spacial Step Size
    real(wp) ::  dxi                                       ! Step Between Adjacent Elements in  xi
    real(wp) :: deta                                       ! Step Between Adjacent Elements in eta

    ! Axis Co-ordinates in Computational Domain
    real(wp), allocatable ::  xi(:)                        ! Vector of  xi-ordinates
    real(wp), allocatable :: eta(:)                        ! Vector of eta-ordinates

#ifdef xGridRefinement
    !  xi(x) Transformation Metrics
    real(wp), allocatable ::  DxiDx(:)                     ! Metric of Transformation - 1st Derivative
    real(wp), allocatable ::  DxiDxsq(:)                   ! We Often Require the Square of the 1st Derivative
    real(wp), allocatable ::  D2xiDx2(:)                   ! Metrics of Transformation - 2nd Derivative
#endif
#ifdef yGridRefinement
    ! eta(y) Transformation Metric
    real(wp), allocatable :: DetaDy(:)                     ! Metric of Transformation - 1st Derivative
    real(wp), allocatable :: DetaDysq(:)                   ! We Often Require the Square of the 1st Derivative
    real(wp), allocatable :: D2etaDy2(:)                   ! Metrics of Transformation - 2nd Derivative
#endif
end type GRID

    ! Grid 1
    type(GRID), target :: grid1

#ifdef MultiGrid
    ! Grid 2
    type(GRID), target :: grid2
    
# ifdef nGrids3
    ! Grid 3
    type(GRID), target :: grid3

#  ifdef nGrids4
    ! Grid 4
    type(GRID), target :: grid4

#   ifdef nGrids5
    ! Grid 5
    type(GRID), target :: grid5

#    ifdef nGrids6
    ! Grid 6
    type(GRID), target :: grid6

#    endif
#   endif
#  endif
# endif
#endif

! Everyone Loves pi
real(wp) :: pi 

contains

subroutine Mesh_MultiGrid(xLim, yLim, xRefinement,yRefinement)
! Mesh_MultiGrid(xLim, yLim, xRefinement,yRefinement) computes the physical
! and 4-Tier computational mesh variables between the limits xLim and yLim.

! Inputs: xLim(/xMin,xMax/), where     xMin is the Lower Bound in the x-Direction
!                                      xMax is the Upper Bound in the x-Direction

!         yLim(/yMin,yMax/), where     yMin is the Lower Bound in the y-Direction
!                                      yMax is the Upper Bound in the y-Direction

!                NOTE: LIMITS MUST BE INPUT TO WORKING PRECISION

!         xRefinement is An Integer:   0 = No Refinement
!                                      1 = Near-Wall Symmetric Hyperbolic Tangent Refinement about y-Axis
!                                      2 = Tangent Refinement about y-Axis

!         yRefinement is An Integer:   0 = No Refinement
!                                      1 = Near-Wall Symmetric Hyperbolic Tangent Refinement about x-Axis
!                                      2 = Tangent Refinement about x-Axis

! Limits of Physical Domain
real(wp) :: xLim(1:2)                                   ! Min and Max Limits in x-Direction
real(wp) :: yLim(1:2)                                   ! Min and Max Limits in y-Direction

! Define x- and yRefinement Inputs
integer, optional :: xRefinement
integer, optional :: yRefinement

! And their counterparts used by the code
integer :: xR
integer :: yR

! Assignment Variables
integer :: i,j, thisGRID               ! Loop Variables
integer :: p, q                        ! Populating Variables for Intermediate and Grid 4

! Grid Refinement Variables
real(wp) :: a, b

! Pointer to point to grid levels
type(GRID), pointer :: grid_ptr => null()

pi = 4.0_wp*ATAN(1.0_wp) 

!===========================================================================================!
!                                Spot of Book Keeping                                       !
!===========================================================================================!

! Check That x- and y-Limits Have Been Specified Correctly
if     (xLim(1) > xLim(2)) then
    write(*,*) 'ERROR: LOWER BOUND OF x-DIRECTION EXCEEDS UPPER BOUND'
    stop
elseif (xLim(1) == xLim(2)) then
    write(*,*) 'ERROR: x-DIRECTION LIMITS ARE EQUAL'
    stop
end if
if (yLim(1) > yLim(2)) then
    write(*,*) 'ERROR: LOWER BOUND OF y-DIRECTION EXCEEDS UPPER BOUND'
    stop
elseif (yLim(1) == yLim(2)) then
    write(*,*) 'ERROR: y-DIRECTION LIMITS ARE EQUAL'
    stop
end if

! Assess the Requested Refinement Type
if     (present(xRefinement)) then
    xR = xRefinement
elseif (.NOT.present(xRefinement)) then
    ! Default - No Refinement
    xR = 0
end if
if (present(yRefinement)) then
    yR = yRefinement
elseif (.NOT.present(yRefinement)) then
    ! Default - No Refinement
    yR = 0
end if

! If x and or y Grid Refinement has not been enabled in PreProcessing Flags, xRefinement and
! yRefinement will be set to zero

#ifndef xGridRefinement
    xR = 0
#endif 
#ifndef yGridRefinement
    yR = 0
#endif 

! Check That x- and yRefinement Inputs Have Been Specified Correctly
if    ((xR < 0).OR.(xR > 2)) then
    write(*,*) 'ERROR: UNKNOWN REFINEMENT IN x-DIRECTION'
    stop
end if
if ((yR < 0).OR.(yR > 2)) then
    write(*,*) 'ERROR: UNKNOWN REFINEMENT IN y-DIRECTION'
    stop
end if

!==========================================================================================!
!                              Fine Computational Mesh                                     !
!==========================================================================================!

! From the definition of the number of elements each tier, the values for the intermediete
! and Grid 4s can be populated from those of the Grid 1.

! Initialise Fine Computational Grid

! Grid 1
    ! Grid 1 Has the Same Number of Grid Points As the Physical Mesh
    grid1%nxi  = nx
    grid1%neta = ny
    ! Allocate Memory for Grid 1
    ! Grid 1 ordinate Vectors
    allocate(grid1%xi (1:grid1%nxi ))
    allocate(grid1%eta(1:grid1%neta))

#ifdef xGridRefinement     
    ! Grid 1 Metrics of Transformation
    allocate(grid1%DxiDx  (1:grid1%nxi))
    allocate(grid1%DxiDxsq(1:grid1%nxi))
    allocate(grid1%D2xiDx2(1:grid1%nxi))
#endif
#ifdef yGridRefinement
    allocate(grid1%DetaDy  (1:grid1%neta))
    allocate(grid1%DetaDysq(1:grid1%neta))
    allocate(grid1%D2etaDy2(1:grid1%neta))
#endif

! Limits of Physical and Fine Computational Meshes are Equal
 xiLim = xLim
etaLim = yLim

! Establish Step Size in Each Dimension
 grid1%dxi =  (xiLim(2)- xiLim(1))/(REAL(grid1%nxi -1,wp))
grid1%deta = (etaLim(2)-etaLim(1))/(REAL(grid1%neta-1,wp))

! Populate Ordinate Vectors Between Their Respective Limits
! xi Ordinate Vector
do i = 1,grid1%nxi
     grid1%xi(i) = xiLim(1) + REAL(i-1,wp)*grid1%dxi 
end do 

! eta Ordinate Vector
do j = 1,grid1%neta
    grid1%eta(j) = etaLim(1) + REAL(j-1,wp)*grid1%deta
end do

! Compute Mesh Refinement according to User Specification

!*********************************************************************************************!
! xRefinement

SELECT CASE (xR)

    CASE (0)
    ! No Refinement
        x = grid1%xi

#ifdef xGridRefinement     
    ! Grid 1 Metrics of Transformation
        grid1%DxiDx  (:) = 1.0_wp
        grid1%DxiDxsq(:) = 1.0_wp
        grid1%D2xiDx2(:) = 0.0_wp

    CASE (1)
    ! Near-Wall Symmetric Hyperbolic Tangent Refinement about y-Axis
        b = 1.3_wp              ! Refinement Factor, b, where 0 < b < 2:
                                !     b --> 0 : Unrefined Case
                                !     b --> 2 : Points Concentrated at Boundaries 
        a = TANH(b*xLim(2))     ! Ensures xi(end) = x(end)

        do i = 1,grid1%nxi

            ! Map Physical Domain x-Ordinate
                        x(i) = TANH(b*grid1%xi(i))/a

            ! Compute Grid Metrics of Transformation, 1st derivative
              grid1%DxiDx(i) = a/(b-b*(x(i)**2.0_wp)*(a**2.0_wp))

            grid1%DxiDxsq(i) = grid1%DxiDx(i)**2.0_wp             

            ! Compute Grid Metrics of Transformation, 2nd derivative
            grid1%D2xiDx2(i) = (2.0_wp*b*x(i)*(a**3.0_wp))/( &
                               (b-b*(x(i)**2.0_wp)*(a**2.0_wp))**2.0_wp)

        end do

    CASE (2)
    ! Tangent Refinement about y-Axis
        write(*,*) 'ERROR: YOU HAVE BEEN LAZY AND NOT CODED THIS YET'
        stop

#endif

END SELECT

!********************************************************************************************!
! yRefinement

SELECT CASE (yR)
    
    CASE (0)
    ! No Refinement
        y =  grid1%eta
    
#ifdef yGridRefinement
        grid1%DetaDy  (:) = 1.0_wp
        grid1%DetaDysq(:) = 1.0_wp
        grid1%D2etaDy2(:) = 0.0_wp

    CASE (1)
    ! Near-Wall Symmetric Hyperbolic Tangent Refinement about x-Axis
        b = 1.3_wp              ! Refinement Factor, b, where 0 < b < 2:
                                !     b --> 0 : Unrefined Case
                                !     b --> 2 : Points Concentrated at Boundaries 
        a = TANH(b*yLim(2))     ! Ensures eta(end) = y(end)

        do j = 1,grid1%neta

            ! Map Physical Domain x-Ordinate
                         y(j) = TANH(b*grid1%eta(j))/a

            ! Compute Grid Metrics of Transformation, 1st derivative
              grid1%DetaDy(j) = a/(b-b*(y(j)**2.0_wp)*(a**2.0_wp))

            grid1%DetaDysq(j) = grid1%DetaDy(j)**2.0_wp
            ! Compute Grid Metrics of Transformation, 2nd derivative
            grid1%D2etaDy2(j) = (2.0_wp*b*y(j)*(a**3.0_wp))/( &
                                (b-b*(y(j)**2.0_wp)*(a**2.0_wp))**2.0_wp)

        end do

    CASE (2)
    ! Tangent Refinement about x-Axis
        write(*,*) 'ERROR: YOU HAVE BEEN LAZY AND NOT CODED THIS YET'
        stop

#endif

END SELECT

!============================================================================================!
!                                 Populate Lower Grid Parameters                             !
!============================================================================================!

#ifdef MultiGrid
# if defined(MultiGrid) && .NOT.defined(nGrids3) && .NOT.defined(nGrids4) && .NOT.defined(nGrids5) && .NOT.defined(nGrids6)
    do thisGRID = 2,2
    
# elif defined(MultiGrid) && defined(nGrids3) && .NOT.defined(nGrids4) && .NOT.defined(nGrids5) && .NOT.defined(nGrids6)
    do thisGRID = 2,3

# elif defined(MultiGrid) && defined(nGrids3) && defined(nGrids4) && .NOT.defined(nGrids5) && .NOT.defined(nGrids6)
    do thisGRID = 2,4

# elif defined(MultiGrid) && defined(nGrids3) && defined(nGrids4) && defined(nGrids5) && .NOT.defined(nGrids6)
    do thisGRID = 2,5

# elif defined(MultiGrid) && defined(nGrids3) && defined(nGrids4) && defined(nGrids5) && defined(nGrids6)
    do thisGRID = 2,6

# endif

    call GetMesh(thisGRID,grid_ptr)

! Grid 2 Allocations
    grid_ptr%nxi  = 2**(m-(thisGRID-1))+1
    grid_ptr%neta = 2**(n-(thisGRID-1))+1
    ! Allocate Memory for Grid
    ! Grid ordinate Vectors
    allocate(grid_ptr%xi (1:grid_ptr%nxi ))
    allocate(grid_ptr%eta(1:grid_ptr%neta))

#ifdef xGridRefinement     
    ! Grid x(xi) Metrics of Transformation
    allocate(grid_ptr%DxiDx  (1:grid_ptr%nxi))
    allocate(grid_ptr%DxiDxsq(1:grid_ptr%nxi))
    allocate(grid_ptr%D2xiDx2(1:grid_ptr%nxi))
#endif
#ifdef yGridRefinement
    ! Grid y(eta) Metrics of Transformation
    allocate(grid_ptr%DetaDy  (1:grid_ptr%neta))
    allocate(grid_ptr%DetaDysq(1:grid_ptr%neta))
    allocate(grid_ptr%D2etaDy2(1:grid_ptr%neta))
#endif

    ! Step Size in Grid
    grid_ptr%dxi  =  (xiLim(2)- xiLim(1))/(REAL(grid_ptr%nxi -1,wp))
    grid_ptr%deta = (etaLim(2)-etaLim(1))/(REAL(grid_ptr%neta-1,wp))
    
    ! Populate Grid Parameters
    p = 0
    q = 0
    do j = 1,grid1%neta,2**(thisGRID-1)
        ! Advance the y index of the Grid
        q = q + 1

        do i = 1,grid1%nxi,2**(thisGRID-1)
            ! Advance the x index of the Grid
            p = p + 1
        
            ! Populate Ordinate Vectors
            grid_ptr%xi (p) = grid1%xi (i)
            grid_ptr%eta(q) = grid1%eta(j)

       ! Populate Grid Metrics of Transformation
#ifdef xGridRefinement       
            grid_ptr%DxiDx  (p) =  grid1%DxiDx  (i)    ! 1st Derivative Metrics
            grid_ptr%DxiDxsq(p) =  grid1%DxiDxsq(i)    ! 1st Derivative Metrics Squared
            grid_ptr%D2xiDx2(p) =  grid1%D2xiDx2(i)    ! 2nd Derivative Metrics
#endif
#ifdef yGridRefinement
            grid_ptr%DetaDy  (q) = grid1%DetaDy  (j)    ! 1st Derivative Metrics
            grid_ptr%DetaDysq(q) = grid1%DetaDysq(j)    ! 1st Derivative Metrics Squared
            grid_ptr%D2etaDy2(q) = grid1%D2etaDy2(j)    ! 2nd Derivative Metrics
#endif
    
        end do
        ! Reset x-index counter
        p = 0
    end do 

    ! nullify pointer
    grid_ptr => null()

end do
#endif


RETURN
end subroutine Mesh_MultiGrid

subroutine GetMesh(TargetGRID, ptr)
! GetMeshValues Takes the Input TargetGRID and Returns the Address of the Requested Grid Type.

! The TargetGRID input is an integer and defines the grid on which PoissonSolver will solve:
!
!   1 --> Grid 1
!   2 --> Grid 2
!   3 --> Grid 3
!   4 --> Grid 4
!   5 --> Grid 5
!   6 --> Grid 6
!   7 --> Grid 7

implicit none

integer :: TargetGRID
type(GRID), pointer :: ptr

! Nullify Pointer before Re-Pointing
ptr => null()

SELECT CASE (TargetGRID)

    CASE (1)
    
        ptr => grid1                       ! Point to Grid 1

#ifdef MultiGrid
    CASE (2)
    
        ptr => grid2                       ! Point to Grid 2 

# ifdef nGrids3    
    CASE (3)
    
        ptr => grid3                       ! Point to Grid 3

#  ifdef nGrids4    
    CASE (4)
    
        ptr => grid4                       ! Point to Grid 4

#   ifdef nGrids5    
    CASE (5)
    
        ptr => grid5                       ! Point to Grid 5        

#    ifdef nGrids6    
    CASE (6)
    
        ptr => grid6                       ! Point to Grid 6        

#    endif
#   endif
#  endif
# endif
#endif

END SELECT

RETURN
end subroutine GetMesh

end module MeshGen