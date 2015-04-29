! Module For Poisson Solver for Stream Function
! Final Year Project
! E. Higham Feb 2014
module Poissonmod
! Declare The Precision to be Used
use SetPrecision
#include <src/PreProcVars.f90>
implicit none

! Declare Poisson Solver Parameters
real(wp) :: tol                                              ! Residual to Converge on
real(wp) :: alpha                                            ! SOR Relaxation Parameter


contains

subroutine Relax(psi,omega,TargetGRID,nIterations,SorR)
#include <src/PreProcVars.f90>
use MeshGen
use SimulationVars
implicit none

! Relax Approximates the Stream Function Numerically using 4th (DEFAULT) or 2nd Order Finite Difference
! Approximations to Derivatives and Successive Over Relaxation (SOR) with relaxation parameter
! alpha to within a residual of tol.

! SOR Method:

!                       psi_{n+1} = {1-alpha)*psi_{n} + alpha*psi_{GS}

! Where n Denotes the Iteration
!       GS Denotes The Gauss-Seidel Method

! The TargetGRID input is an integer and defines the grid on which Relax will solve:
!
!   1 --> Grid 1
!   2 --> Grid 2
!   3 --> Grid 3
!   4 --> Grid 4
!   5 --> Grid 5
!   6 --> Grid 6

!==============================================================================================!
!                                 Declaration of Common Variables                              !
!==============================================================================================!

! Inputs
integer  :: TargetGRID
integer  :: nIterations
real(wp) :: alpha
real(wp) ::   psi(:,:)
real(wp) :: omega(:,:)

! Pointer To Point to the Grid to Solve on
type(GRID), pointer :: grid_ptr => null()

! Temporary Variables To Solve Poisson Equation
! Store as Temporary Variables to Avoid Calling an Entire Data Structure
! Number of Points in the Array
integer :: nXi
integer :: nEta
! Temporary Grid Step Size
! Poisson Solver uses Reciprocal of Step Size
real(wp) :: r1_dXi                                           
real(wp) :: r1_dEta

! And Reciporcal of Step Size Squared
real(wp) :: r1_dXisq
real(wp) :: r1_dEtasq

! Temporary Stream Function Arrays
real(wp), allocatable  ::  psiTemp(:,:)

#ifdef xGridRefinement
real(wp), allocatable ::  DxiDxsq(:)        ! Temporary Transformation Metrics, 1st Derivative
                                            ! --> Note, we require the squares of the 1st derivative
real(wp), allocatable ::  D2xiDx2(:)        ! Temporary Transformation Metrics, 2nd Derivative
#endif
#ifdef yGridRefinement
real(wp), allocatable :: DetaDysq(:)        ! Temporary Transformation Metrics, 1st Derivative
                                            ! --> Note, we require the squares of the 1st derivative
real(wp), allocatable :: D2etaDy2(:)        ! Temporary Transformation Metrics, 2nd Derivative                           
#endif

! Loop Indicies
integer :: i, j, it

! Streamfunction or Residual?
integer :: SorR

real(wp) :: dif1alpha

# ifndef Order2
! 4th Order Poisson's Equation Solver - DEFAULT ==================================================!

! 4th Order Finite Difference Stensil Coefficients
! a - Central Difference Coefficients
! b - Off-Centre Coefficients

#if defined(xGridRefinement).OR.defined(yGridRefinement)
! 1st Derivative
real(wp) :: a1(1:4)
real(wp) :: b1(1:4)
#endif

!2nd Derivative
real(wp) :: a2(1:4)
real(wp) :: b2(1:5)

! Denomenator Array
real(wp), allocatable :: den(:,:)

! Coefficient Arrays
real(wp), allocatable ::   a(:)             ! Coefficients of psi(i-2,j  )
real(wp), allocatable ::   b(:)             ! Coefficients of psi(i-1,j  )
real(wp), allocatable ::   c(:)             ! Coefficients of psi(i+1,j  )
real(wp), allocatable ::   d(:)             ! Coefficients of psi(i+2,j  )
real(wp), allocatable ::   e(:)             ! Coefficients of psi(i  ,j-2)
real(wp), allocatable ::   f(:)             ! Coefficients of psi(i  ,j-1)
real(wp), allocatable ::   g(:)             ! Coefficients of psi(i  ,j+1)
real(wp), allocatable ::   h(:)             ! Coefficients of psi(i  ,j+2)
real(wp)              ::   p(1:5)           ! Forward Off-Center Coefficients in xi - 
real(wp)              ::   q(1:5)           ! Forward Off-Center Coefficients in eta - 
real(wp)              ::   r(1:5)           ! Backward Off-Center Coefficients in xi- 
real(wp)              ::   s(1:5)           ! Backward Off-Center Coefficients in eta -

! Repeated Constants in Poisson Solver
real(wp) :: r6_5
real(wp) :: r3_2
real(wp) :: r2_5
real(wp) :: const

!integer :: it

! Reteated Constants
r6_5 = 6.0_wp/5.0_wp
r3_2 = 3.0_wp/2.0_wp
r2_5 = 2.0_wp/5.0_wp

! ========================================================================================== !
!                                   Memory Allocation                                        !
! ========================================================================================== !

! Find the Specified Grid to Solve on
call GetMesh(TargetGRID,grid_ptr)

! Number of Points in x and y:
 nXi = grid_ptr%nxi
nEta = grid_ptr%neta

! Step Sizes
 r1_dXi = 1.0_wp/grid_ptr%dxi
r1_dEta = 1.0_wp/grid_ptr%deta

! Allocate Temporary Varibales
! Stream Function Arrays
allocate(psiTemp(1:nXi,1:nEta))          ! Stream Function Arrays

#ifdef xGridRefinement
allocate( DxiDxsq(1:nXi ))               ! Transformation Metrics, 1st Derivative
allocate( D2xiDx2(1:nXi ))               ! Transformation Metrics, 2nd Derivative
#endif
#ifdef yGridRefinement
allocate(DetaDysq(1:nEta))               ! Transformation Metrics, 1st Derivative
allocate(D2etaDy2(1:nEta))               ! Transformation Metrics, 2nd Derivative
#endif

! Populate Temporary Variables
 psiTemp(:,:) =  psi(:,:)

#ifdef xGridRefinement
DxiDxsq (:) = grid_ptr%DxiDxsq (:)      ! Transformation Metrics, 1st Derivative
D2xiDx2 (:) = grid_ptr%D2xiDx2 (:)      ! Transformation Metrics, 2nd Derivative
#endif
#ifdef yGridRefinement
DetaDysq(:) = grid_ptr%DetaDysq(:)      ! Transformation Metrics, 1st Derivative
D2etaDy2(:) = grid_ptr%D2etaDy2(:)      ! Transformation Metrics, 2nd Derivative
#endif

    ! Allocate Denomenator
    allocate(den(1:nXi,1:nEta))

    ! Allocate Central Difference Coefficient Arrays
    allocate(a(1:nXi ))
    allocate(b(1:nXi ))
    allocate(c(1:nXi ))
    allocate(d(1:nXi ))
    allocate(e(1:nEta))
    allocate(f(1:nEta))
    allocate(g(1:nEta))
    allocate(h(1:nEta))

!=============================================================================================!
!                                 Poisson Solver SOR Method                                   !
!=============================================================================================!
!
! 2D Computational Domain:
!
!                    _____________________________________________________________ 
! Region            |    |    |                                         |    |    |
! numbers           | BC | BC |                    BC                   | BC | BC |
! correspond        |____|____|_________________________________________|____|____|
! to the            |    |    |                                         |    |    |
! numerical         | BC | 08 |                    01                   | 02 | BC |
! methods           |____|____|_________________________________________|____|____|
! used to           |    |    |                                         |    |    |
! solve at          |    |    |                                         |    |    |
! that point.       |    |    |                                         |    |    |
!                   |    |    |                                         |    |    |
!                   |    |    |                                         |    |    |
!                   |    |    |                                         |    |    |
!                   | BC | 07 |                    09                   | 03 | BC |
!                   |    |    |                                         |    |    |
!                   |    |    |                                         |    |    |
!                   |    |    |                                         |    |    |
!                   |    |    |                                         |    |    |
!                   |    |    |                                         |    |    |
!   Eta             |____|____|_________________________________________|____|____|
!   ^               |    |    |                                         |    |    |
!   |               | BC | 06 |                    05                   | 04 | BC |
!   |               |____|____|_________________________________________|____|____|
!   |               |    |    |                                         |    |    |
!   |               | BC | BC |                    BC                   | BC | BC |
!   *------> Xi     |____|____|_________________________________________|____|____|
!

! Populate Finite Difference Stensil Vectors
#if defined(xGridRefinement).OR.defined(yGridRefinement)
! 1st Derivative
a1 = (/1.0_wp, -8.0_wp,  8.0_wp, -1.0_wp/)/12.0_wp
b1 = (/-3.0_wp, 18.0_wp, -6.0_wp, 1.0_wp/)/12.0_wp
#endif

! 2nd Derivative
a2 = (/-1.0_wp, 16.0_wp, 16.0_wp, -1.0_wp/)/12.0_wp
b2 = (/10.0_wp, -4.0_wp, 14.0_wp, -6.0_wp, 1.0_wp/)/12.0_wp

! Repeated Constants
 r1_dXisq =  r1_dXi*r1_dXi
r1_dEtasq = r1_dEta*r1_dEta

! ===================================================================================== !
!                              Build Coefficient Arrays                                 !
! ===================================================================================== !

! Initialise Coefficients
    a(:) = 0.0_wp
    b(:) = 0.0_wp
    c(:) = 0.0_wp
    d(:) = 0.0_wp
    e(:) = 0.0_wp
    f(:) = 0.0_wp
    g(:) = 0.0_wp
    h(:) = 0.0_wp
    p(:) = 0.0_wp
    q(:) = 0.0_wp
    r(:) = 0.0_wp
    s(:) = 0.0_wp

! xi - Forward Off-Centre Differencing ==========================
#ifdef xGridRefinement
    p(1) = b2(1)*DxiDxsq(2)*r1_dXisq + b1(1)*D2xiDx2(2)*r1_dXi  ! i - 1
    p(2) = b2(2)*DxiDxsq(2)*r1_dXisq + b1(2)*D2xiDx2(2)*r1_dXi  ! i + 1
    p(3) = b2(3)*DxiDxsq(2)*r1_dXisq + b1(3)*D2xiDx2(2)*r1_dXi  ! i + 2
    p(4) = b2(4)*DxiDxsq(2)*r1_dXisq + b1(4)*D2xiDx2(2)*r1_dXi  ! i + 3
    p(5) = b2(5)*DxiDxsq(2)*r1_dXisq                            ! i + 4
#endif
#ifndef xGridRefinement
    p(1) = b2(1)*r1_dXisq   ! i - 1
    p(2) = b2(2)*r1_dXisq   ! i + 1
    p(3) = b2(3)*r1_dXisq   ! i + 2
    p(4) = b2(4)*r1_dXisq   ! i + 3
    p(5) = b2(5)*r1_dXisq   ! i + 4
#endif

! eta - Forward Off-Centre Differencing =========================
#ifdef yGridRefinement
    q(1) = b2(1)*DetaDysq(2)*r1_dEtasq + b1(1)*D2etaDy2(2)*r1_dEta  ! j - 1
    q(2) = b2(2)*DetaDysq(2)*r1_dEtasq + b1(2)*D2etaDy2(2)*r1_dEta  ! j + 1
    q(3) = b2(3)*DetaDysq(2)*r1_dEtasq + b1(3)*D2etaDy2(2)*r1_dEta  ! j + 2
    q(4) = b2(4)*DetaDysq(2)*r1_dEtasq + b1(4)*D2etaDy2(2)*r1_dEta  ! j + 3
    q(5) = b2(5)*DetaDysq(2)*r1_dEtasq                              ! j + 4
#endif
#ifndef yGridRefinement
    q(1) = b2(1)*r1_dEtasq  ! j - 1
    q(2) = b2(2)*r1_dEtasq  ! j + 1
    q(3) = b2(3)*r1_dEtasq  ! j + 2
    q(4) = b2(4)*r1_dEtasq  ! j + 3
    q(5) = b2(5)*r1_dEtasq  ! j + 4
#endif

! xi - Backward Off-Centre Differencing ==========================
#ifdef xGridRefinement
    r(1) = b2(1)*DxiDxsq(nXi-1)*r1_dXisq - b1(1)*D2xiDx2(nXi-1)*r1_dXi ! i + 1
    r(2) = b2(2)*DxiDxsq(nXi-1)*r1_dXisq - b1(2)*D2xiDx2(nXi-1)*r1_dXi ! i - 1
    r(3) = b2(3)*DxiDxsq(nXi-1)*r1_dXisq - b1(3)*D2xiDx2(nXi-1)*r1_dXi ! i - 2
    r(4) = b2(4)*DxiDxsq(nXi-1)*r1_dXisq - b1(4)*D2xiDx2(nXi-1)*r1_dXi ! i - 3
    r(5) = b2(5)*DxiDxsq(nXi-1)*r1_dXisq                           ! i - 4
#endif
#ifndef xGridRefinement
    r(1) = b2(1)*r1_dXisq  ! i + 1
    r(2) = b2(2)*r1_dXisq  ! i - 1
    r(3) = b2(3)*r1_dXisq  ! i - 2
    r(4) = b2(4)*r1_dXisq  ! i - 3
    r(5) = b2(5)*r1_dXisq  ! i - 4
#endif

! eta - Backward Off-Centre Differencing =========================
#ifdef yGridRefinement
    s(1) = b2(1)*DetaDysq(nEta-1)*r1_dEtasq - b1(1)*D2etaDy2(nEta-1)*r1_dEta  ! j + 1
    s(2) = b2(2)*DetaDysq(nEta-1)*r1_dEtasq - b1(2)*D2etaDy2(nEta-1)*r1_dEta  ! j - 1
    s(3) = b2(3)*DetaDysq(nEta-1)*r1_dEtasq - b1(3)*D2etaDy2(nEta-1)*r1_dEta  ! j - 2
    s(4) = b2(4)*DetaDysq(nEta-1)*r1_dEtasq - b1(4)*D2etaDy2(nEta-1)*r1_dEta  ! j - 3
    s(5) = b2(5)*DetaDysq(nEta-1)*r1_dEtasq                              ! j - 4
#endif
#ifndef yGridRefinement
    s(1) = b2(1)*r1_dEtasq  ! j + 1 
    s(2) = b2(2)*r1_dEtasq  ! j - 1
    s(3) = b2(3)*r1_dEtasq  ! j - 2
    s(4) = b2(4)*r1_dEtasq  ! j - 3
    s(5) = b2(5)*r1_dEtasq  ! j - 4
#endif


! xi - Central Differencing =====================================

do i = 3,(nXi-2)
#ifdef xGridRefinement
    a(i) = a2(1)*DxiDxsq(i)*r1_dXisq + a1(1)*D2xiDx2(i)*r1_dXi
    b(i) = a2(2)*DxiDxsq(i)*r1_dXisq + a1(2)*D2xiDx2(i)*r1_dXi
    c(i) = a2(3)*DxiDxsq(i)*r1_dXisq + a1(3)*D2xiDx2(i)*r1_dXi
    d(i) = a2(4)*DxiDxsq(i)*r1_dXisq + a1(4)*D2xiDx2(i)*r1_dXi
#endif
#ifndef xGridRefinement
    a(i) = a2(1)*r1_dXisq
    b(i) = a2(2)*r1_dXisq
    c(i) = a2(3)*r1_dXisq
    d(i) = a2(4)*r1_dXisq
#endif
end do

! eta - Central Differencing ====================================

do j = 3,(nEta-2)
#ifdef yGridRefinement
    e(j) = a2(1)*DetaDysq(j)*r1_dEtasq + a1(1)*D2etaDy2(j)*r1_dEta
    f(j) = a2(2)*DetaDysq(j)*r1_dEtasq + a1(2)*D2etaDy2(j)*r1_dEta
    g(j) = a2(3)*DetaDysq(j)*r1_dEtasq + a1(3)*D2etaDy2(j)*r1_dEta
    h(j) = a2(4)*DetaDysq(j)*r1_dEtasq + a1(4)*D2etaDy2(j)*r1_dEta
#endif
#ifndef yGridRefinement
    e(j) = a2(1)*r1_dEtasq
    f(j) = a2(2)*r1_dEtasq
    g(j) = a2(3)*r1_dEtasq
    h(j) = a2(4)*r1_dEtasq
#endif
end do

! ===================================================================================== !
!                              Build Denominator Array                                  !
! ===================================================================================== !

    den(:,:) = 0.0_wp

! Region (1)

       j = (nEta-1)
    do i = 3,(nXi-2)

#if .NOT.defined(xGridRefinement) && .NOT.defined(yGridRefinement)
     den(i,j) = r6_5/(3.0_wp*r1_dXisq+r3_2*r1_dEtasq)
#endif

#if .NOT.defined(xGridRefinement) && defined(yGridRefinement)
    den(i,j) = r6_5/(3.0_wp*r1_dXisq                            + &
               r3_2*DetaDysq(j)*r1_dEtasq - D2etaDy2(j)*r1_dEta)
#endif

#if defined(xGridRefinement) && .NOT.defined(yGridRefinement)
    den(i,j) = r6_5/(3.0_wp*DxiDxsq(i)*r1_dXisq+r3_2*r1_dEtasq)
#endif

#if defined(xGridRefinement) && defined(yGridRefinement)

    den(i,j) = r6_5/(3.0_wp*DxiDxsq (i)*r1_dXisq               + &
                   r3_2*DetaDysq(j)*r1_dEtasq - D2etaDy2(j)*r1_dEta)
#endif

    end do

! Region (2)

    i = (nXi-1)
    j = (nEta-1)

#if .NOT.defined(xGridRefinement) && .NOT.defined(yGridRefinement)
    den(i,j) = r6_5/(r3_2*r1_dXisq + r3_2*r1_dEtasq)
#endif

#if .NOT.defined(xGridRefinement) && defined(yGridRefinement)
    den(i,j) = r6_5/(r3_2*r1_dXisq + &
        r3_2*DetaDysq(j)*r1_dEtasq-D2etaDy2(j)*r1_dEta)
#endif

#if defined(xGridRefinement) && .NOT.defined(yGridRefinement)
    den(i,j) = r6_5/(r3_2*DxiDxsq(i)*r1_dXisq + &
               r3_2*r1_dEtasq-D2xiDx2(i)*r1_dXi)
#endif

#if defined(xGridRefinement) && defined(yGridRefinement)
    den(i,j) = r6_5/(r3_2*DxiDxsq(i)*r1_dXisq + &
        r3_2*DetaDysq(j)*r1_dEtasq-D2xiDx2(i)*r1_dXi-D2etaDy2(j)*r1_dEta)
#endif

! Region (3)

    i = (nXi-1)
    do j = 3,(nEta-2)

#if .NOT.defined(xGridRefinement) && .NOT.defined(yGridRefinement)
    den(i,j) = r6_5/(r3_2*r1_dXisq + 3.0_wp*r1_dEtasq)
#endif

#if .NOT.defined(xGridRefinement) && defined(yGridRefinement)
    den(i,j) = r6_5/(r3_2*r1_dXisq + 3.0_wp*DetaDysq(j)*r1_dEtasq)
#endif

#if defined(xGridRefinement) && .NOT.defined(yGridRefinement)
    den(i,j) = r6_5/(r3_2*DxiDxsq(i)*r1_dXisq + &
               3.0_wp*r1_dEtasq-D2xiDx2(i)*r1_dXi)
#endif

#if defined(xGridRefinement) && defined(yGridRefinement)
    den(i,j) = r6_5/(r3_2*DxiDxsq(i)*r1_dXisq + &
               3.0_wp*DetaDysq(j)*r1_dEtasq-D2xiDx2(i)*r1_dXi)
#endif

    end do

! Region (4)

    i = (nXi-1)
    j = 2
#if .NOT.defined(xGridRefinement) && .NOT.defined(yGridRefinement)
    den(i,j) = r6_5/(r3_2*r1_dXisq + r3_2*r1_dEtasq)
#endif

#if .NOT.defined(xGridRefinement) && defined(yGridRefinement)
    den(i,j) = r6_5/(r3_2*r1_dXisq + &
               r3_2*DetaDysq(j)*r1_dEtasq+D2etaDy2(j)*r1_dEta)
#endif

#if defined(xGridRefinement) && .NOT.defined(yGridRefinement)
    den(i,j) = r6_5/(r3_2*DxiDxsq(i)*r1_dXisq + &
               r3_2*r1_dEtasq-D2xiDx2(i)*r1_dXi)
#endif

#if defined(xGridRefinement) && defined(yGridRefinement)
    den(i,j) = r6_5/(r3_2*DxiDxsq(i)*r1_dXisq + &
        r3_2*DetaDysq(j)*r1_dEtasq-D2xiDx2(i)*r1_dXi+D2etaDy2(j)*r1_dEta)
#endif

! Region (5)

    j = 2
    do i = 3,(nXi-2)

#if .NOT.defined(xGridRefinement) && .NOT.defined(yGridRefinement)
    den(i,j) = r6_5/(3.0_wp*r1_dXisq + r3_2*r1_dEtasq)
#endif

#if .NOT.defined(xGridRefinement) && defined(yGridRefinement)
    den(i,j) = r6_5/(3.0_wp*r1_dXisq + &
               r3_2*DetaDysq(j)*r1_dEtasq+D2etaDy2(j)*r1_dEta)
#endif

#if defined(xGridRefinement) && .NOT.defined(yGridRefinement)
    den(i,j) = r6_5/(3.0_wp*DxiDxsq(i)*r1_dXisq + r3_2*r1_dEtasq)
#endif

#if defined(xGridRefinement) && defined(yGridRefinement)
    den(i,j) = r6_5/(3.0_wp*DxiDxsq(i)*r1_dXisq + &
               r3_2*DetaDysq(j)*r1_dEtasq+D2etaDy2(j)*r1_dEta)
#endif

    end do

! Region (6)

    i = 2
    j = 2

#if .NOT.defined(xGridRefinement) && .NOT.defined(yGridRefinement)
    den(i,j) = r6_5/(r3_2*r1_dXisq + r3_2*r1_dEtasq)
#endif

#if .NOT.defined(xGridRefinement) && defined(yGridRefinement)
    den(i,j) = r6_5/(r3_2*r1_dXisq + & 
        r3_2*DetaDysq(j)*r1_dEtasq+D2etaDy2(j)*r1_dEta)
#endif

#if defined(xGridRefinement) && .NOT.defined(yGridRefinement)
    den(i,j) = r6_5/(r3_2*DxiDxsq(i)*r1_dXisq + & 
               r3_2*r1_dEtasq+D2xiDx2(i)*r1_dXi)
#endif

#if defined(xGridRefinement) && defined(yGridRefinement)
    den(i,j) = r6_5/(r3_2*DxiDxsq(i)*r1_dXisq + & 
        r3_2*DetaDysq(j)*r1_dEtasq+D2xiDx2(i)*r1_dXi+D2etaDy2(j)*r1_dEta)
#endif

! Region (7)

    i = 2
    do j = 3,(nEta-2)

#if .NOT.defined(xGridRefinement) && .NOT.defined(yGridRefinement)
    den(i,j) = r6_5/(r3_2*r1_dXisq+3.0_wp*r1_dEtasq)
#endif

#if .NOT.defined(xGridRefinement) && defined(yGridRefinement)
    den(i,j) = r6_5/(r3_2*r1_dXisq+3.0_wp*DetaDysq(j)*r1_dEtasq)
#endif

#if defined(xGridRefinement) && .NOT.defined(yGridRefinement)
    den(i,j) = r6_5/(r3_2*DxiDxsq(i)*r1_dXisq + &
               3.0_wp*r1_dEtasq+D2xiDx2(i)*r1_dXi)
#endif

#if defined(xGridRefinement) && defined(yGridRefinement)
    den(i,j) = r6_5/(r3_2*DxiDxsq(i)*r1_dXisq + &
        3.0_wp*DetaDysq(j)*r1_dEtasq+D2xiDx2(i)*r1_dXi)
#endif

    end do

! Region (8)
    i = 2
    j = (nEta-1)

#if .NOT.defined(xGridRefinement) && .NOT.defined(yGridRefinement)
    den(i,j) = r6_5/(r3_2*r1_dXisq + r3_2*r1_dEtasq)
#endif

#if .NOT.defined(xGridRefinement) && defined(yGridRefinement)
    den(i,j) = r6_5/(r3_2*r1_dXisq + &
        r3_2*DetaDysq(j)*r1_dEtasq-D2etaDy2(j)*r1_dEta)
#endif

#if defined(xGridRefinement) && .NOT.defined(yGridRefinement)
    den(i,j) = r6_5/(r3_2*DxiDxsq(i)*r1_dXisq + &
        r3_2*r1_dEtasq+D2xiDx2(i)*r1_dXi)
#endif

#if defined(xGridRefinement) && defined(yGridRefinement)
    den(i,j) = r6_5/(r3_2*DxiDxsq(i)*r1_dXisq + &
        r3_2*DetaDysq(j)*r1_dEtasq+D2xiDx2(i)*r1_dXi-D2etaDy2(j)*r1_dEta)
#endif

! Region (9)

    do j = 3,(nEta-2)
        do i = 3,(nXi-2)

#if .NOT.defined(xGridRefinement) && .NOT.defined(yGridRefinement)
    den(i,j) = r2_5/(r1_dXisq + r1_dEtasq)
#endif

#if .NOT.defined(xGridRefinement) && defined(yGridRefinement)
    den(i,j) = r2_5/(r1_dXisq + DetaDysq(j)*r1_dEtasq)
#endif

#if defined(xGridRefinement) && .NOT.defined(yGridRefinement)
    den(i,j) = r2_5/(DxiDxsq(i)*r1_dXisq + r1_dEtasq)
#endif

#if defined(xGridRefinement) && defined(yGridRefinement)
    den(i,j) = r2_5/(DxiDxsq(i)*r1_dXisq + DetaDysq(j)*r1_dEtasq)
#endif

        end do
    end do

! Along Lid (x,y=end) , u = d psi/dy = 1
#ifdef yGridRefinement
    const = (12.0_wp*grid_ptr%deta)/grid_ptr%DetaDy(nEta)
#endif
#ifndef yGridRefinement
    const = 12.0_wp*grid_ptr%deta
#endif


! Gauss-Siedel Iteration ============================================================

! Relaxation Parameter
    alpha = 1.18_wp
dif1alpha = 1.0_wp-alpha

! Solving Vorticity Stream-function or residual-error?
! Apply Boundary Conditions to Regions Marked BC

do it = 1,nIterations

    ! Apply Boundary Conditions to Regions Marked BC
    
    SELECT CASE(SorR)
        CASE(0)
            !call psiBC(psiTemp,TargetGRID)
            psi(nXi    ,:       ) = 0.0_wp
            psi(1      ,:       ) = 0.0_wp  
            psi(:      ,1       ) = 0.0_wp
            psi(:      ,nEta    ) = const
        CASE(1)
            psiTemp(nXi,:   ) = 0.0_wp
            psiTemp(1  ,:   ) = 0.0_wp
            psiTemp(:  ,1   ) = 0.0_wp
            psiTemp(:  ,nEta) = 0.0_wp
    END SELECT

! Region (1) Central Xi, Backward Off-Centre Eta
       j = (nEta-1)
    do i = 3,(nXi-2)
            
    psiTemp(i,j) = omega(i,j) + &
!
        a(i)*psiTemp(i-2,j)+b(i)*psiTemp(i-1,j)                     + &
        c(i)*psiTemp(i+1,j)+d(i)*psiTemp(i+2,j)                     + &
!
        s(1)*psiTemp(i,j+1)+s(2)*psiTemp(i,j-1)+s(3)*psiTemp(i,j-2) + &
        s(4)*psiTemp(i,j-3)+s(5)*psiTemp(i,j-4)

        psiTemp(i,j) = psiTemp(i,j)*den(i,j)
            
    end do

    ! Region (2) Backward Off-Centre Xi, Backward Off-Centre Eta
    i = (nXi-1)
    j = (nEta-1)
            
    psiTemp(i,j) = omega(i,j) + &
!
        r(1)*psiTemp(i+1,j)+r(2)*psiTemp(i-1,j)+r(3)*psiTemp(i-2,j) + &
        r(4)*psiTemp(i-3,j)+r(5)*psiTemp(i-4,j)                     + &
!
        s(1)*psiTemp(i,j+1)+s(2)*psiTemp(i,j-1)+s(3)*psiTemp(i,j-2) + &
        s(4)*psiTemp(i,j-3)+s(5)*psiTemp(i,j-4)
              
    psiTemp(i,j) = psiTemp(i,j)*den(i,j)

    ! Region (3) Backward Off-Centre Xi, Central Eta
    i = (nXi-1)
 do j = 3,(nEta-2)
            
    psiTemp(i,j) = omega(i,j) + &
!
        r(1)*psiTemp(i+1,j)+r(2)*psiTemp(i-1,j)+r(3)*psiTemp(i-2,j) + &
        r(4)*psiTemp(i-3,j)+r(5)*psiTemp(i-4,j)                     + &
!
        e(j)*psiTemp(i,j-2)+f(j)*psiTemp(i,j-1)                     + &
        g(j)*psiTemp(i,j+1)+h(j)*psiTemp(i,j+2)                               
            
    psiTemp(i,j) = psiTemp(i,j)*den(i,j)

 end do
    
    ! Region (4) Backward Off-Centre Xi, Forward Off-Centre Eta
    i = (nXi-1)
    j = 2
            
    psiTemp(i,j) = omega(i,j) + &
!
        r(1)*psiTemp(i+1,j)+r(2)*psiTemp(i-1,j)+r(3)*psiTemp(i-2,j) + &
        r(4)*psiTemp(i-3,j)+r(5)*psiTemp(i-4,j)                     + &
!
        q(1)*psiTemp(i,j-1)+q(2)*psiTemp(i,j+1)+q(3)*psiTemp(i,j+2) + &
        q(4)*psiTemp(i,j+3)+q(5)*psiTemp(i,j+4)

    psiTemp(i,j) = psiTemp(i,j)*den(i,j)
            
    ! Region (5) Central Xi, Forward Off-Centre Eta
    j = 2
 do i = 3,(nXi-2)
            
    psiTemp(i,j) = omega(i,j) + &
!
        a(i)*psiTemp(i-2,j)+b(i)*psiTemp(i-1,j) + &
        c(i)*psiTemp(i+1,j)+d(i)*psiTemp(i+2,j) + &
!
        q(1)*psiTemp(i,j-1)+q(2)*psiTemp(i,j+1)+q(3)*psiTemp(i,j+2) + &
        q(4)*psiTemp(i,j+3)+q(5)*psiTemp(i,j+4)
            
        psiTemp(i,j) = psiTemp(i,j)*den(i,j)

 end do
    
    ! Region (6) Forward Off-Centre Xi, Forward Off-Centre Eta
    i = 2
    j = 2
            
    psiTemp(i,j) = omega(i,j) + &
!
        p(1)*psiTemp(i-1,j)+p(2)*psiTemp(i+1,j)+p(3)*psiTemp(i+2,j) + &
        p(4)*psiTemp(i+3,j)+p(5)*psiTemp(i+4,j)                     + &
!
        q(1)*psiTemp(i,j-1)+q(2)*psiTemp(i,j+1)+q(3)*psiTemp(i,j+2) + &
        q(4)*psiTemp(i,j+3)+q(5)*psiTemp(i,j+4)
            
    psiTemp(i,j) = psiTemp(i,j)*den(i,j)

    ! Region (7) Forward Off-Centre Xi, Central Eta
    i = 2
 do j = 3,(nEta-2)
            
    psiTemp(i,j) = omega(i,j) + &
!
        p(1)*psiTemp(i-1,j)+p(2)*psiTemp(i+1,j)+p(3)*psiTemp(i+2,j) + &
        p(4)*psiTemp(i+3,j)+p(5)*psiTemp(i+4,j)                     + &
!
        e(j)*psiTemp(i,j-2)+f(j)*psiTemp(i,j-1)                     + &
        g(j)*psiTemp(i,j+1)+h(j)*psiTemp(i,j+2)

    psiTemp(i,j) = psiTemp(i,j)*den(i,j)
            
 end do
    
    ! Region (8) Forward Off-Centre Xi, Backward Off-Centre Eta
    i = 2
    j = (nEta-1)
            
    psiTemp(i,j) = omega(i,j) + &
!
        p(1)*psiTemp(i-1,j)+p(2)*psiTemp(i+1,j)+p(3)*psiTemp(i+2,j) + &
        p(4)*psiTemp(i+3,j)+p(5)*psiTemp(i+4,j)                     + &
!
        s(1)*psiTemp(i,j+1)+s(2)*psiTemp(i,j-1)+s(3)*psiTemp(i,j-2) + &
        s(4)*psiTemp(i,j-3)+s(5)*psiTemp(i,j-4)

    psiTemp(i,j) = psiTemp(i,j)*den(i,j)

    ! Region (9) Central Xi, Central Eta
 do j = 3,(nEta-2)
    do i = 3,(nXi-2)

            
    psiTemp(i,j) = omega(i,j) + &
!
        a(i)*psiTemp(i-2,j)+b(i)*psiTemp(i-1,j) + &
        c(i)*psiTemp(i+1,j)+d(i)*psiTemp(i+2,j) + &
!
        e(j)*psiTemp(i,j-2)+f(j)*psiTemp(i,j-1) + &
        g(j)*psiTemp(i,j+1)+h(j)*psiTemp(i,j+2)

    psiTemp(i,j) = psiTemp(i,j)*den(i,j)
            
    end do
 end do
        
    !Over Relaxation
    psiTemp(:,:) = dif1alpha*psi(:,:)+alpha*psiTemp(:,:)

    psi(:,:) = psiTemp(:,:)


end do

    ! Tidy Up
    grid_ptr    => null()

    deallocate(psiTemp,a,b,c,d,e,f,g,h,den)

#ifdef xGridRefinement
    deallocate( DxiDxsq)               ! Transformation Metrics, 1st Derivative
    deallocate( D2xiDx2)               ! Transformation Metrics, 2nd Derivative
#endif
#ifdef yGridRefinement
    deallocate(DetaDysq)              ! Transformation Metrics, 1st Derivative
    deallocate(D2etaDy2)              ! Transformation Metrics, 2nd Derivative
#endif

!=================================================================================================!
# endif

! 2nd Order Poisson's Equation Solver ============================================================!
# ifdef Order2

! psi_{GS} is a 2nd Order, 5-point Gauss-Seidel Iteration, of the form:
!  
!                    [       omega(i  ,j  )                     + ]
! psi(i,j) =  a(i,j) [    b(i)*psi(i+1,j  ) + c(i)*psi(i-1,j  ) + ]
!                    [    d(j)*psi(i  ,j-1) + e(j)*psi(i  ,j-1)   ] 
!
! where
!
!           1  [xi_{x}(i)^2   eta_{y}(j)^2 ]^{-1}  
! a(i,j) =  -* [----------- + ------------ ]     ; 
!           2  [   dxi^2         deta^2    ]       
!
!               xi_{xx}(i)      xi_{x}(i)^2                     xi_{x}(i)^2      xi_{xx}(i)    
! b(i) =       ----------   +  ------------      ;  c(j)   =    -----------   -  ---------- ; 
!                 2*dxi            dxi^2                          dxi^2            2*dxi       
!
!             eta_{y}(j)^2      eta_{yy}(j)                    eta_{y}(j)^2       eta_{yy}(j) 
!  d(j) =     ------------  +   -----------      ;  e(j)   =   ------------   -   ----------- 
!                deta^2           2*deta                        deta^2              2*deta    
!
! and 
!
!  phi_{z}(k) is the kth element in a generic 1st derivative grid transformation metric
! phi_{zz}(k) is the kth element in a generic 2nd derivative grid transformation metric

!==============================================================================================!
!                              2nd Order Variable Declaration                                  !
!==============================================================================================!

! Coefficient Arrays - See Comments Above for Definition
real(wp), allocatable :: a(:,:)
real(wp), allocatable :: b(:)
real(wp), allocatable :: c(:)
real(wp), allocatable :: d(:)
real(wp), allocatable :: e(:)

! Repeated Constants

! Half The Reciporcal of the Step Size
real(wp) :: r1_2dXi
real(wp) :: r1_2dEta

! Iteration Counter
!integer :: it

! others
integer :: step
integer :: p,q

! ========================================================================================== !
!                                   Memory Allocation                                        !
! ========================================================================================== !

! Find the Specified Grid to Solve on
call GetMesh(TargetGRID,grid_ptr)

! Number of Points in x and y:
 nXi = grid_ptr%nxi
nEta = grid_ptr%neta

! Step Sizes
 r1_dXi = 1.0_wp/grid_ptr%dxi
r1_dEta = 1.0_wp/grid_ptr%deta 

! Allocate Temporary Varibales
! Stream Function Arrays
allocate(psiTemp(1:nXi,1:nEta))          ! Stream Function Arrays

#ifdef xGridRefinement
allocate( DxiDxsq(1:nXi ))               ! Transformation Metrics, 1st Derivative
allocate( D2xiDx2(1:nXi ))               ! Transformation Metrics, 2nd Derivative
#endif
#ifdef yGridRefinement
allocate(DetaDysq(1:nEta))               ! Transformation Metrics, 1st Derivative
allocate(D2etaDy2(1:nEta))               ! Transformation Metrics, 2nd Derivative
#endif

! Allocate Coefficient Arrays
allocate(a(1:nXi,1:nEta))
allocate(b(1:nXi ))
allocate(c(1:nXi ))
allocate(d(1:nEta))
allocate(e(1:nEta))

! Populate Temporary Variables
 psiTemp(:,:) =  psi(:,:)

#ifdef xGridRefinement
DxiDxsq (:) = grid_ptr%DxiDxsq (:)      ! Transformation Metrics, 1st Derivative
D2xiDx2 (:) = grid_ptr%D2xiDx2 (:)      ! Transformation Metrics, 2nd Derivative
#endif
#ifdef yGridRefinement
DetaDysq(:) = grid_ptr%DetaDysq(:)      ! Transformation Metrics, 1st Derivative
D2etaDy2(:) = grid_ptr%D2etaDy2(:)      ! Transformation Metrics, 2nd Derivative
#endif
!=============================================================================================!
!                                 Poisson Solver SOR Method                                   !
!=============================================================================================!
!
! 2D Computational Domain:
!
!                    _____________________________________________________________ 
!                   |    |                                                   |    |
!                   | BC |                         BC                        | BC |
!                   |____|___________________________________________________|____|
!                   |    |                                                   |    |
!                   |    |                                                   |    |
!                   |    |                                                   |    |
!                   |    |                                                   |    |
!                   |    |                                                   |    |
!                   |    |                                                   |    |
!                   |    |                    FLUID DOMAIN                   |    |
!                   |    |                                                   |    |
!                   | BC |                SOLVE IN THIS REGION               | BC |
!                   |    |                                                   |    |
!                   |    |                                                   |    |
!                   |    |                                                   |    |
!   Eta             |    |                                                   |    |
!   ^               |    |                                                   |    |
!   |               |    |                                                   |    |
!   |               |____|___________________________________________________|____|
!   |               |    |                                                   |    |
!   |               | BC |                         BC                        | BC |
!   *------> Xi     |____|___________________________________________________|____|
!

! Repeated Constants
 r1_dXisq =  r1_dXi*r1_dXi
r1_dEtasq = r1_dEta*r1_dEta

 r1_2dXi =  r1_dXi*0.5_wp
r1_2dEta = r1_dEta*0.5_wp

! Optimum alpha for SOR
alpha = 1.7_wp/(1.0_wp+SIN(pi/(REAL(nXi+1,wp))))
!print*, 'alpha = ', alpha
dif1alpha = 1.0_wp-alpha



! ===================================================================================== !
!                            Build Repeated Coefficient Arrays                          !
! ===================================================================================== !

! Initialise arrays a,b,c,d & e amd set to zero
a(:,:) = 0.0_wp
  b(:) = 0.0_wp
  c(:) = 0.0_wp
  d(:) = 0.0_wp
  e(:) = 0.0_wp

do j = 2,(nEta-1)
    do i = 2,(nXi-1)

#if .NOT.defined(xGridRefinement) && .NOT.defined(yGridRefinement)
        a(i,j) = 0.5_wp/(r1_dXisq + r1_dEtasq)
#endif
#if defined(xGridRefinement) && .NOT.defined(yGridRefinement)
        a(i,j) = 0.5_wp/(DxiDxsq(i)*r1_dXisq + r1_dEtasq)
#endif
#if .NOT.defined(xGridRefinement) && defined(yGridRefinement)
        a(i,j) = 0.5_wp/(r1_dXisq + DetaDysq(j)*r1_dEtasq)
#endif
#if defined(xGridRefinement) && defined(yGridRefinement)
        a(i,j) = 0.5_wp/(DxiDxsq(i)*r1_dXisq + DetaDysq(j)*r1_dEtasq)
#endif

    end do
end do

do i = 2,(nXi-1)

#ifndef xGridRefinement
    b(i) = r1_dXisq
    c(i) = r1_dXisq
#endif
#ifdef xGridRefinement
    b(i) = DxiDxsq(i)*r1_dXisq + D2xiDx2(i)*r1_2dXi
    c(i) = DxiDxsq(i)*r1_dXisq - D2xiDx2(i)*r1_2dXi
#endif

end do

do j = 2,(nEta-1)

#ifndef yGridRefinement
    d(j) = r1_dEtasq
    e(j) = r1_dEtasq
#endif
#ifdef yGridRefinement
    d(j) = DetaDysq(j)*r1_dEtasq + D2etaDy2(j)*r1_2dEta
    e(j) = DetaDysq(j)*r1_dEtasq - D2etaDy2(j)*r1_2dEta
#endif

end do


! Gauss-Siedel Iteration ============================================================

do it = 1,nIterations

    ! Apply Boundary Conditions to Regions 10 to 25
    call psiBC(psi,TargetGRID)

    ! SOR Scheme with 2nd Order Gauss-Seidel Iteration

    do j = 2,(nEta-1)
        do i = 2,(nXi-1)

            psiTemp(i,j) =  dif1alpha*psi(i  ,j  ) + &
!
                      alpha*a(i,j)*(omega(i  ,j  ) + &
                             b(i)*psiTemp(i+1,j  ) + c(i)*psiTemp(i-1,j  ) + &
                             d(j)*psiTemp(i  ,j+1) + e(j)*psiTemp(i  ,j-1))

        end do

        ! Apply Line Relaxation
!        psiTemp(:,j) = (1.0_wp-alpha)*psi(:,j)+alpha*psiTemp(:,j)

    end do

    psi(:,:) = psiTemp(:,:)

end do

!=================================================================================================!
# endif

RETURN

end subroutine Relax

subroutine MG(psi,omega)
use MeshGen
implicit none

! MG takes the inputs
real(wp) ::   psi(:,:)      ! 2D Stream Function Array
real(wp) :: omega(:,:)      ! 2D Vorticity Array

#ifndef MultiGrid
! If MultiGrid is not enabled, solve on Grid 1 only
    call Relax(psi,omega,1,tol)
#endif

#ifdef MultiGrid

# ifndef nGrids3

!======================================================================================!
!                                 2 TIER V- MULTIGRID                                  !
!======================================================================================!

! MG solves the 2D Poisson's equation calling Relax on 2 interstitial grids of
! varying resolution using a V-MultiGrid method.
!
! Grid 1 -------o-----o
!                \   /
!                 \ /
! Grid 2 ----------o

!
! Where Grid labels have their usual meaning
!
!    Algorythm: (Moin, P., Fundamentals of Engineering Numerical Analysis)
! 1. Perform a few iterations on the original equation, Aφ = b, on the fine grid
!    with the mesh spacing h. Let the resulting solution be denoted by ψ. Calculate
!    the residual r = b − Aψ on the same grid.

! 2. Transfer the residual to a coarse grid (restriction) of mesh spacing 2h, and on
!    this grid iterate on the error equation A*epsilon = r, with the initial guess 
!    epsilon = 0

! 3. Interpolate (prolongation) the resulting psi to the fine grid.Make a correction
!    on the previous ψ by adding it to epsilon, i.e., ψnew = ψ + epsilon. Use ψnew as the
!    initial guess to iterate on the original problem, Aφ = b.

! 4. Repeat the process.

! Also require l additional temporary arrays of the same size as the l coarser meshes
real(wp) ::  epsilon(grid1%nxi,grid1%neta),  rho(grid1%nxi,grid1%neta)
real(wp) :: epsilon2(grid2%nxi,grid2%neta), rho2(grid2%nxi,grid2%neta)

! Convergence Criterion (Max Residual)
real(wp) :: ERR

! Initialise Convergence Criterion
ERR = 1.0_wp

! Initialse epsilon 2 as zeros
epsilon2(:,:) = 0.0_wp

do while (ERR > tol)    
    
    ! STEP 1 =============================================================================!
    
    !  Relax psi n times on Grid 1
    call Relax(psi,omega,1,10,0)

    ! Compute the Residual
    rho = Residual(omega,psi,1)
    
    ! STEP 2 =============================================================================!
                                        
    ! Restrict the residual onto grid 2
    call Restrict(1,rho,rho2)
    
    ! Relax epsilon2 on Grid 2
    call Relax(epsilon2,rho2,2,100,1)
    
    ! STEP 3 =============================================================================!
    
    ! Interpolate the error onto grid 1
    epsilon = Interp2(2,epsilon2)

    ! Calculate the new stream function
    psi = psi + epsilon
    
    ! Relax on Grid 1
    call Relax(psi,omega,1,10,0)

    ! STEP 4 =============================================================================!

    ! Compute the residual
    rho = Residual(omega,psi,1)

    ERR = MAXVAL(ABS(rho))

    !print*, ERR
end do

# endif

#  ifdef nGrids3
#   ifndef nGrids4

!======================================================================================!
!                               3 TIER FULL MULTIGRID                                  !
!======================================================================================!

! MG solves the 2D Poisson's equation calling Relax on 3 interstitial grids of
! varying resolution using a Full MultiGrid method.
!
! Grid 1 ------------------o-----------o
!                         / \         /
!                        /   \       /
! Grid 2 ---------o-----o-----o-----o
!                / \   /       \   /
!               /   \ /         \ /
! Grid 3 ------o-----o-----------o

!
! Where Grid labels have their usual meaning
!
!    Algorythm:

!  1. Restrict (by sampling) omega onto grid 3 and converge to tol
!  2. Interpolate psi onto grid 2 and converge to tol
!  3. Richardson Extrapolation on grids 2->3 and converge to tol
!  4. Update grid 2 and converge to tol
!  5. Interpolate onto grid 1 and converge to tol
!  6. Richardson Extrapolation on grids 1->2 and converge to tol
!  7. Richardson Extrapolation on grids 2->3 and converge to tol
!  8. Update grid 2 and converge to tol
!  9. Update grid 1 and converge to tol

#if .NOT.defined(xGridRefinement) .OR. .NOT.defined(yGridRefinement)
    ! Also require l additional temporary arrays of the same size as the l coarser meshes
    real(wp) :: psi2(grid2%nxi,grid2%neta), omega2(grid2%nxi,grid2%neta)
    real(wp) :: psi3(grid3%nxi,grid3%neta), omega3(grid3%nxi,grid3%neta)
#endif

! error and residual arrays
real(wp) ::  epsilon(grid1%nxi,grid1%neta),  rho(grid1%nxi,grid1%neta)
real(wp) :: epsilon2(grid2%nxi,grid2%neta), rho2(grid2%nxi,grid2%neta)
real(wp) :: epsilon3(grid3%nxi,grid3%neta), rho3(grid3%nxi,grid3%neta)

! Residual of residual array
real(wp) :: rhorho2(grid2%nxi,grid2%neta)

! Need Loop indices for Restriction
integer :: i,j,p,q

! Convergence Criterion (Max Residual)
real(wp) :: ERR

! Initialise Convergence Criterion
ERR = 1.0_wp

! Initialise epsilon 2 as zeros
epsilon2(:,:) = 0.0_wp
epsilon3(:,:) = 0.0_wp

#if .NOT.defined(xGridRefinement) .OR. .NOT.defined(yGridRefinement)
    ! Initialise psi3 as zeros
    psi3(:,:) = 0.0_wp
    
    ! Sample omega onto grids 2
    p = 0
    q = 0
    do j = 1,grid1%neta,2
        q = q + 1
        do i = 1,grid1%nxi,2
                  p = p + 1
        omega2(p,q) = omega(i,j)            ! Restrict Vorticity Field
        end do
    p = 0
    end do
    
    ! Sample omega onto grids 2
    p = 0
    q = 0
    do j = 1,grid1%neta,4
        q = q + 1
        do i = 1,grid1%nxi,4
                  p = p + 1
        omega3(p,q) = omega(i,j)            ! Restrict Vorticity Field
        end do
    p = 0
    end do
      
    
    ! STEP 1 =============================================================================!
    ! 1. Restrict psi and omega onto Grid 3 by sampling
    
    p = 0
    q = 0
    do j = 1,grid1%neta,4
        q = q + 1
        do i = 1,grid1%nxi,4
                  p = p + 1
          psi3(p,q) = psi(i,j)    ! Restrict Stream Funciton Field 
        end do
    p = 0
    end do
    
    ! Converge Relax on Grid 3
    call Relax(psi3,omega3,3,300,0)
    
    ! STEP 2 =============================================================================!
                                        
    ! 2. Interpolate onto Grid 2
    
    psi2 = Interp2(3,psi3)
    
    ! Converge Relax on Grid 2
    call Relax(psi2,omega2,2,50,0)
    
    ! Compute the Residual
    rho2 = Residual(omega2,psi2,2)
    
    ! STEP 3 =============================================================================!
    
    ! 3. Restrict the residual onto grid 3
    call Restrict(2,rho2,rho3)
    
    ! Converge Relax on Grid 3
    call Relax(epsilon3,rho3,3,100,1)
    
    ! STEP 4 =============================================================================!
    
    ! 4. Interpolate the error onto grid 2
    epsilon2 = Interp2(3,epsilon3)
    
    ! Calculate new stream function array
    psi2 = psi2 + epsilon2
    ! Converge Relax on Grid 2
    call Relax(psi2,omega2,2,50,0)
    
    ! STEP 5 =============================================================================!
    
    ! 5. Interpolate onto grid 1 and converge to tol
    
    ! Interpolate psi2 onto Grid 1
    psi = Interp2(2,psi2)
#endif

! Converge Relax on Grid 1
call Relax(psi,omega,1,5,0)

do while (ERR > tol)  

    rho = Residual(omega,psi,1)
    
    call Restrict(1,rho,rho2)
    
    ! STEP 6 =============================================================================!
    
    ! 6. Relax epsilon2 on grid 2
    epsilon2(:,:) = 0.0_wp
    
    ! Converge Relax on Grid 2
    call Relax(epsilon2,rho2,2,25,1)
    
    ! Compute the residual
    rhorho2 = Residual(rho2,epsilon2,2)
    
    ! Restrict onto grid 3
    call Restrict(2,rhorho2,rho3)
    
    ! STEP 7 =============================================================================!
    
    epsilon3(:,:) = 0.0_wp
    
    ! 7. Relax epsilon3 on grid 3
    call Relax(epsilon3,rho3,3,700,1)
    
    ! Interpolate epsilon3 onto grid 2
    epsilon2 = epsilon2+Interp2(3,epsilon3)
    
    ! STEP 8 =============================================================================!
    
    ! 8. Relax epsilon on grid 2
    
    call Relax(epsilon2,rho2,2,25,1)
    
    ! Interpolate epsilon onto grid 1
    epsilon = Interp2(2,epsilon2)
    
    ! STEP 9 =============================================================================!
    
    ! 9. Update Grid 1 and Converge
    psi = psi + epsilon
    
    ! Converge Relax on Grid 1
    call Relax(psi,omega,1,5,0)

    ! Compute the residual
    rho = Residual(omega,psi,1)

    ERR = MAXVAL(ABS(rho))

    !print*, ERR
end do

#   endif
#  endif

#  ifdef nGrids3 
#   ifdef nGrids4
#    ifndef nGrids5

!======================================================================================!
!                               4 TIER FULL MULTIGRID                                  !
!======================================================================================!

! MG solves the 2D Poisson's equation calling Relax on 4 interstitial grids of
! varying resolution using a Full MultiGrid method.
!
! Grid 1 -----------------------------o-----------------o
!                                    / \               /
!                                   /   \             /
! Grid 2 --------------o-----------o-----o-----------o
!                     / \         /       \         /
!                    /   \       /         \       /
! Grid 3 -----o-----o-----o-----o-----------o-----o
!            / \   /       \   /             \   /
!           /   \ /         \ /               \ /
! Grid 4 --o-----o-----------o-----------------o
!
! Where Grid labels have their usual meaning
!
#if .NOT.defined(xGridRefinement) .OR. .NOT.defined(yGridRefinement)
! Also require l additional temporary arrays of the same size as the l coarser meshes
real(wp) :: psi2(grid2%nxi,grid2%neta), omega2(grid2%nxi,grid2%neta)
real(wp) :: psi3(grid3%nxi,grid3%neta), omega3(grid3%nxi,grid3%neta)
real(wp) :: psi4(grid4%nxi,grid4%neta), omega4(grid4%nxi,grid4%neta)
#endif

! error and residual arrays
real(wp) ::  epsilon(grid1%nxi,grid1%neta),  rho(grid1%nxi,grid1%neta)
real(wp) :: epsilon2(grid2%nxi,grid2%neta), rho2(grid2%nxi,grid2%neta)
real(wp) :: epsilon3(grid3%nxi,grid3%neta), rho3(grid3%nxi,grid3%neta)
real(wp) :: epsilon4(grid4%nxi,grid4%neta), rho4(grid4%nxi,grid4%neta)

! residual of residual arrays 
real(wp) :: rhorho2(grid2%nxi,grid2%neta)
real(wp) :: rhorho3(grid3%nxi,grid3%neta)

! Need Loop indices for Restriction
integer :: i,j,p,q

! Convergence Criterion (Max Residual)
real(wp) :: ERR

! Initialise Convergence Criterion
ERR = 1.0_wp

! Initialise epsilon 2 as zeros
epsilon2(:,:) = 0.0_wp
epsilon3(:,:) = 0.0_wp
epsilon4(:,:) = 0.0_wp

#if .NOT.defined(xGridRefinement) .OR. .NOT.defined(yGridRefinement)
    ! Initialise psi3 as zeros
    psi4(:,:) = 0.0_wp
    
    !! Sample omega onto grids 2
    p = 0
    q = 0
    do j = 1,grid1%neta,2
        q = q + 1
        do i = 1,grid1%nxi,2
                  p = p + 1
        omega2(p,q) = omega(i,j)            ! Restrict Vorticity Field
        end do
    p = 0
    end do
    
    ! Sample omega onto grid 3
    p = 0
    q = 0
    do j = 1,grid1%neta,4
        q = q + 1
        do i = 1,grid1%nxi,4
                  p = p + 1
        omega3(p,q) = omega(i,j)            ! Restrict Vorticity Field
        end do
    p = 0
    end do
      
    ! Sample omega onto grid 4
    p = 0
    q = 0
    do j = 1,grid1%neta,8
        q = q + 1
        do i = 1,grid1%nxi,8
                  p = p + 1
        omega4(p,q) = omega(i,j)            ! Restrict Vorticity Field
        end do
    p = 0
    end do
    
    ! STEP 1 =============================================================================!
    ! 1. Restrict psi and omega onto Grid 4 by sampling
    
    p = 0
    q = 0
    do j = 1,grid1%neta,8
        q = q + 1
        do i = 1,grid1%nxi,8
                  p = p + 1
          psi4(p,q) = psi(i,j)    ! Restrict Stream Funciton Field 
        end do
    p = 0
    end do
    
    ! Converge Relax on Grid 4
    call Relax(psi4,omega4,4,200,0)
    
    ! STEP 2 =============================================================================!
                                        
    ! 2. Interpolate onto Grid 3
    psi3 = Interp2(4,psi4)
    
    ! Converge Relax on Grid 3
    call Relax(psi3,omega3,3,20,0)
    
    ! Compute the Residual
    rho3 = Residual(omega3,psi3,3)
    
    ! STEP 3 =============================================================================!
    
    ! 3. Restrict the residual onto grid 4
    call Restrict(3,rho3,rho4)
    
    ! Converge Relax on Grid 4
    call Relax(epsilon4,rho4,4,50,1)
    
    ! STEP 4 =============================================================================!
    
    ! 4. Interpolate the error onto grid 3
    epsilon3 = Interp2(4,epsilon4)
    
    ! Calculate new stream function array
    psi3 = psi3 + epsilon3
    
    ! Converge Relax on Grid 3
    call Relax(psi3,omega3,3,50,0)
    
    ! STEP 5 =============================================================================!
    
    ! Interpolate psi2 onto Grid 2
    psi2 = Interp2(3,psi3)
    
    ! Relax on Grid 2
    call Relax(psi2,omega2,2,5,0)
    
    ! Compute Residual
    rho2 = Residual(omega2,psi2,2)
    
    ! STEP 6 =============================================================================!
    
    ! Restrict the residual onto grid 3
    call Restrict(2,rho2,rho3)
    
    epsilon3(:,:) = 0.0_wp
    ! Relax the error
    call Relax(epsilon3,rho3,3,20,1)
    
    ! Find the residual of the residual
    rhorho3 = Residual(rho3,epsilon3,3)
    
    ! STEP 7 =============================================================================!
    
    ! Restrict onto grid 4
    call Restrict(3,rhorho3,rho4)
    
    ! Relax
    call Relax(epsilon4,rho4,4,100,1)
    
    ! STEP 8 =============================================================================!
    
    ! Interpolate onto grid 3
    epsilon3 = epsilon3 + Interp2(4,epsilon4)
    
    ! Relax
    call Relax(epsilon3,rho3,3,40,1)
    
    ! STEP 9 =============================================================================!
    
    ! Interpolate onto grid 2
    epsilon2 = Interp2(3,epsilon3)
    
    ! Update stream function array
    psi2 = psi2 + epsilon2
    
    ! Relax on Grid 2
    call Relax(psi2,omega2,2,5,0)
    ! STEP 10 =============================================================================!
    ! Interpolate onto grid 1
    psi = Interp2(2,psi2)
#endif
    
! Relax on grid 1
call Relax(psi,omega,1,5,0)

! V-MultiGrid ============================================================================!

do while (ERR > tol)  

    rho = Residual(omega,psi,1)
    
    call Restrict(1,rho,rho2)
    
    ! STEP 11 ============================================================================!
    
    epsilon2(:,:) = 0.0_wp
    
    ! Converge Relax on Grid 2
    call Relax(epsilon2,rho2,2,10,1)
    
    ! Compute the residual
    rhorho2 = Residual(rho2,epsilon2,2)
    
    ! Restrict onto grid 3
    call Restrict(2,rhorho2,rho3)
    
    ! STEP 12 ============================================================================!
    
    epsilon3(:,:) = 0.0_wp
    
    ! 7. Relax epsilon3 on grid 3
    call Relax(epsilon3,rho3,3,20,1)
    
    ! Compute residual
    rhorho3 = Residual(rho3,epsilon3,3)

    
    ! STEP 13 ============================================================================!

    ! Restrict onto grid 4
    call Restrict(3,rhorho3,rho4)

    epsilon4(:,:) = 0.0_wp

    ! Relax
    call Relax(epsilon4,rho4,4,200,1)

    ! STEP 14 ============================================================================!

    ! Interpolate epsilon4 onto grid 3
    epsilon3 = epsilon3+Interp2(4,epsilon4)
    
    ! Relax epsilon on grid 3
    call Relax(epsilon3,rho3,3,40,1)
    
    ! STEP 15  ===========================================================================!
    
    ! 9. Update Grid 2 and Converge
    epsilon2 = epsilon2 + Interp2(3,epsilon3)
    
    ! Converge Relax on Grid 2
    call Relax(epsilon2,rho2,2,20,1)

    ! STEP 16  ===========================================================================!

    ! Interpolate epsilon onto grid 1
    epsilon = Interp2(2,epsilon2)
    
    ! 9. Update Grid 1 and Converge
    psi = psi + epsilon
    
    ! Converge Relax on Grid 1
    call Relax(psi,omega,1,10,0)

    ! Compute the residual
    rho = Residual(omega,psi,1)

    ERR = MAXVAL(ABS(rho))

    !print*, ERR
end do

#    endif
#   endif
#  endif

#  ifdef nGrids3 
#   ifdef nGrids4
#    ifdef nGrids5
#     ifndef nGrids6

!======================================================================================!
!                               5 TIER FULL MULTIGRID                                  !
!======================================================================================!

! MG solves the 2D Poisson's equation calling Relax on 5 interstitial grids of
! varying resolution using a Full MultiGrid method.
!
! Grid 1 --------------------------------------------------o-----------------------o 
!                                                         / \                     /
!                                                        /   \                   /
! Grid 2 -----------------------------o-----------------o-----o-----------------o
!                                    / \               /       \               /
!                                   /   \             /         \             /
! Grid 3 --------------o-----------o-----o-----------o-----------o-----------o 
!                     / \         /       \         /             \         /    
!                    /   \       /         \       /               \       /         
! Grid 4 -----o-----o-----o-----o-----------o-----o-----------------o-----o              
!            / \   /       \   /             \   /                   \   /          
!           /   \ /         \ /               \ /                     \ /              
! Grid 5 --o-----o-----------o-----------------o-----------------------o           
!
! Where Grid labels have their usual meaning

#if .NOT.defined(xGridRefinement) .OR. .NOT.defined(yGridRefinement)
    ! Also require l additional temporary arrays of the same size as the l coarser meshes
    real(wp) :: psi2(grid2%nxi,grid2%neta), omega2(grid2%nxi,grid2%neta)
    real(wp) :: psi3(grid3%nxi,grid3%neta), omega3(grid3%nxi,grid3%neta)
    real(wp) :: psi4(grid4%nxi,grid4%neta), omega4(grid4%nxi,grid4%neta)
    !real(wp) :: psi5(grid5%nxi,grid5%neta), omega5(grid5%nxi,grid5%neta)
#endif

! error and residual arrays
real(wp) ::  epsilon(grid1%nxi,grid1%neta),  rho(grid1%nxi,grid1%neta)
real(wp) :: epsilon2(grid2%nxi,grid2%neta), rho2(grid2%nxi,grid2%neta)
real(wp) :: epsilon3(grid3%nxi,grid3%neta), rho3(grid3%nxi,grid3%neta)
real(wp) :: epsilon4(grid4%nxi,grid4%neta), rho4(grid4%nxi,grid4%neta)
real(wp) :: epsilon5(grid5%nxi,grid5%neta), rho5(grid5%nxi,grid5%neta)

! residual of residual arrays 
real(wp) :: rhorho2(grid2%nxi,grid2%neta)
real(wp) :: rhorho3(grid3%nxi,grid3%neta)
real(wp) :: rhorho4(grid4%nxi,grid4%neta)
! Need Loop indices for Restriction
integer :: i,j,p,q

! Timing Variables
!real(wp) :: start, finish

! Convergence Criterion (Max Residual)
real(wp) :: ERR

! Initialise Convergence Criterion
ERR = 1.0_wp

! Initialise epsilon 2 as zeros
epsilon2(:,:) = 0.0_wp
epsilon3(:,:) = 0.0_wp
epsilon4(:,:) = 0.0_wp

#if .NOT.defined(xGridRefinement) .OR. .NOT.defined(yGridRefinement)
     ! Initialise psi3 as zeros
    psi4(:,:) = 0.0_wp
    
    !! Sample omega onto grids 2
    p = 0
    q = 0
    do j = 1,grid1%neta,2
        q = q + 1
        do i = 1,grid1%nxi,2
                  p = p + 1
        omega2(p,q) = omega(i,j)            ! Restrict Vorticity Field
        end do
    p = 0
    end do
    
    ! Sample omega onto grid 3
    p = 0
    q = 0
    do j = 1,grid1%neta,4
        q = q + 1
        do i = 1,grid1%nxi,4
                  p = p + 1
        omega3(p,q) = omega(i,j)            ! Restrict Vorticity Field
        end do
    p = 0
    end do
      
    ! Sample omega onto grid 4
    p = 0
    q = 0
    do j = 1,grid1%neta,8
        q = q + 1
        do i = 1,grid1%nxi,8
                  p = p + 1
        omega4(p,q) = omega(i,j)            ! Restrict Vorticity Field
        end do
    p = 0
    end do
    
    ! STEP 1 =============================================================================!
    ! 1. Restrict psi and omega onto Grid 4 by sampling
    
    p = 0
    q = 0
    do j = 1,grid1%neta,8
        q = q + 1
        do i = 1,grid1%nxi,8
                  p = p + 1
          psi4(p,q) = psi(i,j)    ! Restrict Stream Funciton Field 
        end do
    p = 0
    end do
    
    ! Converge Relax on Grid 4
    call Relax(psi4,omega4,4,200,0)
    
    ! STEP 2 =============================================================================!
                                        
    ! 2. Interpolate onto Grid 3
    psi3 = Interp2(4,psi4)
    
    ! Converge Relax on Grid 3
    call Relax(psi3,omega3,3,20,0)
    
    ! Compute the Residual
    rho3 = Residual(omega3,psi3,3)
    
    ! STEP 3 =============================================================================!
    
    ! 3. Restrict the residual onto grid 4
    call Restrict(3,rho3,rho4)
    
    ! Converge Relax on Grid 4
    call Relax(epsilon4,rho4,4,50,1)
    
    ! STEP 4 =============================================================================!
    
    ! 4. Interpolate the error onto grid 3
    epsilon3 = Interp2(4,epsilon4)
    
    ! Calculate new stream function array
    psi3 = psi3 + epsilon3
    
    ! Converge Relax on Grid 3
    call Relax(psi3,omega3,3,50,0)
    
    ! STEP 5 =============================================================================!
    
    ! Interpolate psi2 onto Grid 2
    psi2 = Interp2(3,psi3)

    ! Relax on Grid 2
    call Relax(psi2,omega2,2,5,0)
    
    ! Compute Residual
    rho2 = Residual(omega2,psi2,2)
    
    ! STEP 6 =============================================================================!
    
    ! Restrict the residual onto grid 3
    call Restrict(2,rho2,rho3)
    
    epsilon3(:,:) = 0.0_wp
    ! Relax the error
    call Relax(epsilon3,rho3,3,20,1)
    
    ! Find the residual of the residual
    rhorho3 = Residual(rho3,epsilon3,3)
    
    ! STEP 7 =============================================================================!
    
    ! Restrict onto grid 4
    call Restrict(3,rhorho3,rho4)
    
    ! Relax
    call Relax(epsilon4,rho4,4,100,1)
    
    ! STEP 8 =============================================================================!
    
    ! Interpolate onto grid 3
    epsilon3 = epsilon3 + Interp2(4,epsilon4)
    
    ! Relax
    call Relax(epsilon3,rho3,3,40,1)
    
    ! STEP 9 =============================================================================!
    
    ! Interpolate onto grid 2
    epsilon2 = Interp2(3,epsilon3)
    
    ! Update stream function array
    psi2 = psi2 + epsilon2
    
    ! Relax on Grid 2
    call Relax(psi2,omega2,2,5,1)
    ! STEP 10 =============================================================================!
    ! Interpolate onto grid 1
    psi = Interp2(2,psi2)
#endif

! Relax on grid 1
call Relax(psi,omega,1,5,0)

! V-MultiGrid ============================================================================!

do while (ERR > tol)  

    rho = Residual(omega,psi,1)
    
    call Restrict(1,rho,rho2)
    
    ! STEP 28 ============================================================================!
    
    epsilon2(:,:) = 0.0_wp
    
    ! Converge Relax on Grid 2
    call Relax(epsilon2,rho2,2,5,1)
    
    ! Compute the residual
    rhorho2 = Residual(rho2,epsilon2,2)
    
    ! Restrict onto grid 3
    call Restrict(2,rhorho2,rho3)
    
    ! STEP 29 ============================================================================!
    
    epsilon3(:,:) = 0.0_wp
    
    ! 7. Relax epsilon3 on grid 3
    call Relax(epsilon3,rho3,3,10,1)
    
    ! Compute residual
    rhorho3 = Residual(rho3,epsilon3,3)

    
    ! STEP 20 ============================================================================!

    ! Restrict onto grid 4
    call Restrict(3,rhorho3,rho4)

    epsilon4(:,:) = 0.0_wp

    ! Relax
    call Relax(epsilon4,rho4,4,30,1)

    ! Compute residual
    rhorho4 = Residual(rho4,epsilon4,4)


    ! STEP 21 ============================================================================!

    ! Restrict onto grid 5
    call Restrict(4,rhorho4,rho5)

    epsilon5(:,:) = 0.0_wp

    ! Relax
    call Relax(epsilon5,rho5,5,200,1)

    ! STEP 22 ============================================================================!

    epsilon4 = epsilon4 + Interp2(5,epsilon5)

    ! Relax
    call Relax(epsilon4,rho4,4,40,1)

    ! STEP 23 ============================================================================!

    ! Interpolate epsilon4 onto grid 3
    epsilon3 = epsilon3+Interp2(4,epsilon4)
    
    ! Relax epsilon on grid 3
    call Relax(epsilon3,rho3,3,30,1)
    
    ! STEP 24  ===========================================================================!
    
    ! 9. Update Grid 2 and Converge
    epsilon2 = epsilon2 + Interp2(3,epsilon3)
    
    ! Converge Relax on Grid 2
    call Relax(epsilon2,rho2,2,20,1)

    ! STEP 25  ===========================================================================!

    ! Interpolate epsilon onto grid 1
    epsilon = Interp2(2,epsilon2)
    
    ! 9. Update Grid 1 and Converge
    psi = psi + epsilon
    
    ! Converge Relax on Grid 1
    call Relax(psi,omega,1,10,0)

    ! Compute the residual
    rho = Residual(omega,psi,1)

    ERR = MAXVAL(ABS(rho))

    !print*, ERR


end do

#     endif
#    endif
#   endif
#  endif

#endif

RETURN
end subroutine MG

function Interp2(srcGrid, A)
! Interp2 Bilinearly Interpolates array A from the grid specified by integer srcGrid
! onto a finer grid with size(2**(n+1)+1), where n is an integer exponent of the srcGrid.
use MeshGen
implicit none

! Variable Declaration
real(wp), allocatable :: Interp2(:,:)  ! Output

integer :: srcGrid
! Specifies the Computational Grid of A
! Can Take Values :     2 -> Grid 2
!                       3 -> Grid 3
!                       4 -> Grid 4

real(wp) :: A(:,:)                         ! Array to be interpolated onto target grid

! Declare Grid Pointers
type(GRID), pointer :: src_ptr => null()              ! Source Grid Pointer
type(GRID), pointer :: tgt_ptr => null()              ! Target Grid Pointer

integer :: i, j                                       ! Loop Index Variables
integer :: p, q                                       ! Index Variables


! Find the Meshes to interpolate between
    call GetMesh(srcGrid  , src_ptr)                  ! Get the source grid
    call GetMesh(srcGrid-1, tgt_ptr)                  ! Get the target mesh

allocate(Interp2(1:tgt_ptr%nxi,1:tgt_ptr%neta))

!=======================================================================================================!
!                                     Bilinear Interpolation                                            !
!=======================================================================================================!

! Cell Arrangement

!
!       0----*----0             KEY:
!       |    |    |             0 = Point Shared By Coarse and Fine Grid
!       |    |    |             * = Point Belonging to Grid l Only.
!  eta  *----*----*
!   ^   |    |    |             Fine Grid Step Size = 1/2 Coarse Grid Step Size
!   |   |    |    |
!   |   0----*----0             Fine Grid has 3 Points for every 2 Coarse Points
!   |
!   *-----> xi
!

! Mapping from coarser grid to finer

!           i --> 2i - 1

p = 0
q = 0
! Direct Mapping Between Grids - Values at equal (xi,eta) Co-ordinates are the same
do j = 1,tgt_ptr%neta,2
q = q + 1                           ! Increment Coarser Grid xi Index
    do i = 1,tgt_ptr%nxi,2

    p = p + 1                       ! Increment Coarser Grid eta Index
    Interp2(i,j) = A(p,q)           ! Store Every Coarser Grid Value into Evey Other Finer Grid Entry

    end do
p = 0                               ! Reset Coarser Grid eta Index
end do

! Linear Interpolation Along Lines of Constant y
do j = 1,tgt_ptr%neta,2
    do i = 2,(tgt_ptr%nxi-1),2

    Interp2(i,j) = 0.5_wp*(Interp2(i+1,j) + Interp2(i-1,j))

    end do
end do


! Linear Interpolation Along Lines of Constant x
do j = 2,(tgt_ptr%neta-1),2
    do i = 1,tgt_ptr%nxi,2

    Interp2(i,j) = 0.5_wp*(Interp2(i,j+1) + Interp2(i,j-1))

    end do
end do

! Bilinear Interpolation (x=even,y=even)
do j = 2,(tgt_ptr%neta-1),2
    do i = 2,(tgt_ptr%nxi-1),2

    Interp2(i,j) = 0.25_wp*(Interp2(i+1,j+1) + Interp2(i+1,j-1) + &
                            Interp2(i-1,j+1) + Interp2(i-1,j-1))

    end do
end do

! Nullify Pointer
src_ptr => null()
tgt_ptr => null()


RETURN
end function Interp2

subroutine CopyValues(srcGrid,u_c,u_f)
! CopyValues copies elemests from coarse array u_c from the grid specified by integer srcGrid
! onto every other element on a finer grid array u_f
use MeshGen
implicit none

! Variable Declaration
integer :: srcGrid
! The srcGrid input is an integer and defines the grid on which Relax will solve:
!
!   1 --> Grid 1
!   2 --> Grid 2
!   3 --> Grid 3
!   4 --> Grid 4
!   5 --> Grid 5
!   6 --> Grid 6         

real(wp) :: u_c(:,:)                         ! Array to be copied onto target grid
real(wp) :: u_f(:,:)                         ! Finer Array to recieve Coped Values

! Declare Grid Pointers
type(GRID), pointer :: tgt_ptr => null()              ! Target Grid Pointer

integer :: i, j                                       ! Loop Index Variables
integer :: p, q                                       ! Index Variables


! Find the Meshes to copy between
    call GetMesh(srcGrid-1, tgt_ptr)                  ! Get the target mesh

! Cell Arrangement

!
!       0----*----0             KEY:
!       |    |    |             0 = Point Shared By Coarse and Fine Grid
!       |    |    |             * = Point Belonging to Grid l Only.
!  eta  *----*----*
!   ^   |    |    |             Fine Grid Step Size = 1/2 Coarse Grid Step Size
!   |   |    |    |
!   |   0----*----0             Fine Grid has 3 Points for every 2 Coarse Points
!   |
!   *-----> xi
!

! Mapping from coarser grid to finer

!           i --> 2i - 1

p = 0
q = 0

! Direct Mapping Between Grids - Values at equal (xi,eta) Co-ordinates are the same
do j = 1,tgt_ptr%neta,2
q = q + 1                         ! Increment Coarser Grid xi Index
    do i = 1,tgt_ptr%nxi,2

    p = p + 1                     ! Increment Coarser Grid eta Index
    u_f(i,j) = u_c(p,q)             ! Store Every Coarser Grid Value into Evey Other Finer Grid Entry
    
    end do
p = 0                             ! Reset Coarser Grid eta Index
end do

! Nullify Pointer
tgt_ptr => null()

RETURN
end subroutine CopyValues

subroutine Restrict(srcGrid,u_f,u_c)
! Restriction Operation from a fine grid array u_f on gris srcGrid to a coarser grid u_c
! The srcGRID input is an integer and defines the grid on which Relax will solve:
!
!   1 --> Grid 1
!   2 --> Grid 2
!   3 --> Grid 3
!   4 --> Grid 4
!   5 --> Grid 5
!   6 --> Grid 6

! Assuming the grid on which u_c is based is is 2x coarser than that of u_f, he restriction
! operation is a follows:

!                    /    u_f(2i-1,2j-1) + u_f(2i+1,2j-1) + u_f(2i-1,2j+1) + u_f(2i+1,2j+1)   + \
! u_c(i,j) = r1_16 * | 2*[u_f(2i-1,2j  ) + u_f(2i+1,2j  ) + u_f(2i  ,2j-1) + u_f(2i  ,2j+1)]  + |
!                    \  4*u_f(2i  ,2j  )                                                        /

! where r1_16 = 1/16

! This subroutine is used to restrict the residual in the MG method and so it is assumed that residuals lying 
! on a boundary are zero, thus we resrict from 2,nx-1 and 2,ny-1 in the coarsest grid. Note restrictions occur on
! computational mesh.

use MeshGen
implicit none

! Inputs
integer  :: srcGrid         ! Grid which is being restricted FROM
real(wp) :: u_f(:,:)        ! Fine grid variables to be restricted
real(wp) :: u_c(:,:)        ! Coarse Grid variabled to be determined

! We require a grid pointer to find which grid we a restricting to and from
type(GRID), pointer :: ptr => null()

! And loop indices
integer :: i,j,p,q

! The restriction operation uses the following repeated constant
real(wp) :: r1_16

! Point at the finer grid
call GetMesh(srcGrid,ptr)

! Define 1/16
r1_16 = 1.0_wp/16.0_wp

p = 1
q = 1
do j = 2,ptr%neta-1,2
    q = q + 1
    do i = 2,ptr%nxi-1,2
        p = p + 1
    
        u_c(p,q) =         u_f(i-1,j-1) + u_f(i+1,j-1) + u_f(i-1,j+1) + u_f(i+1,j+1)   + &
                   2.0_wp*(u_f(i-1,j  ) + u_f(i+1,j  ) + u_f(i  ,j-1) + u_f(i  ,j+1))  + &
                    4.0_wp*u_f(i  ,j  )
 

        u_c(p,q) = r1_16*u_c(p,q) 


    end do
p = 1
end do

! Again, we assume that the residual at the boundaries is zero
call GetMesh(srcGrid+1,ptr)

u_c(1      ,:       ) = 0.0_wp
u_c(ptr%nxi,:       ) = 0.0_wp
u_c(:      ,1       ) = 0.0_wp
u_c(:      ,ptr%neta) = 0.0_wp

! Nullift pointer
ptr => null()


RETURN
end subroutine Restrict

function Residual(omega,psi,srcGrid)
! Residual computes the residual, r, of Poisson's equation after the nth relaxation:
! Poisson's eqution for stream funtion and vorticity is

! omega = -(d psi/dx + d psi/dy)

! And the residual is defined as

! r = omega + A*psi
!
! where A is the 4th order finite difference approximation to the laplacian.

! The srcGRID input is an integer and defines the grid on which Relax will solve:
!
!   1 --> Grid 1
!   2 --> Grid 2
!   3 --> Grid 3
!   4 --> Grid 4
!   5 --> Grid 5
!   6 --> Grid 6

use MeshGen
use Derivatives
implicit none

! Inputs
real(wp) ::   omega(:,:)
real(wp) ::     psi(:,:)
integer  :: srcGrid

! Outputs
real(wp), allocatable :: Residual(:,:)

! Need a poiter to point to the grid on which to calculate the residual:
type(GRID), pointer :: ptr => null()

! And we use the following temporary variables for the derivatives of
! psi w.r.t. x and y
real(wp), allocatable :: psi_xx(:,:)
real(wp), allocatable :: psi_yy(:,:)

#ifdef xGridRefinement
    ! Need loop index to apply x(xi) metrics of transformation
    integer :: i
    ! And we must also compute the 1st derivative
    real(wp), allocatable :: psi_x(:,:)
#endif
#ifdef yGridRefinement
    ! Need loop index to apply y(eta) metrics of transformation
    integer :: j
    ! And we must also compute the 1st derivative
    real(wp), allocatable :: psi_y(:,:)
#endif

! Point to the grid on which we calculate the residual
call GetMesh(srcGrid,ptr)

! Allocate memory for output and temporary variables
allocate(Residual(1:ptr%nxi,1:ptr%neta))
allocate(  psi_xx(1:ptr%nxi,1:ptr%neta))
allocate(  psi_yy(1:ptr%nxi,1:ptr%neta))

! Compute the 2nd derivatives
psi_xx = D2f(psi,1,srcGrid)     ! Calcualte the 2nd derivative w.r.t. x
psi_yy = D2f(psi,2,srcGrid)     ! Calcualte the 2nd derivative w.r.t. y

#ifdef xGridRefinement
    ! Must compute 1st derivatives too..
    ! Allocate Temporary Variables
    allocate(psi_x(1:ptr%nxi,1:ptr%neta))
    ! Compute 1st derivative
    psi_x = Df(psi,1,srcGrid)

    ! Apply finite difference stencils to 1st and 2nd derivatives
    do i = 1,ptr%nxi
        psi_x (i,:) = psi_x (i,:)*ptr%D2xiDx2(i)
        psi_xx(i,:) = psi_xx(i,:)*ptr%DxiDxsq(i)
    end do

    ! Combine 1st and 2nd derivatives

    psi_xx = psi_xx + psi_x

#endif

#ifdef yGridRefinement
    ! Must compute 1st derivatives too..
    ! Allocate Temporary Variables
    allocate(psi_y(1:ptr%nxi,1:ptr%neta))
    ! Compute 1st derivative
    psi_y = Df(psi,2,srcGrid)

    ! Apply finite difference stencils to 1st and 2nd derivatives
    do j = 1,ptr%neta
        psi_y (:,j) = psi_y (:,j)*ptr%D2etaDy2(j)
        psi_yy(:,j) = psi_yy(:,j)*ptr%DetaDysq(j)
    end do

    ! Combine 1st and 2nd derivatives

    psi_yy = psi_yy + psi_y

#endif

! Compute residual

Residual = omega + psi_xx + psi_yy

! Ensure that the residual at the boundaries is zero
Residual(1      ,:       ) = 0.0_wp
Residual(ptr%nxi,:       ) = 0.0_wp
Residual(:      ,1       ) = 0.0_wp
Residual(:      ,ptr%neta) = 0.0_wp

ptr => null()

RETURN
end function Residual


end module Poissonmod