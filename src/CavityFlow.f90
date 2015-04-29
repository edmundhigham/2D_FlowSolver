! Module for Vortex Rebound Simulation
! (a.k.a. Normal Dipole No-Slip Boundary Collision)
! E. Higham Mar 2014
! Final Year Project - 2D Laminar Flow Solver
module CavityFlow
#include <src/PreProcVars.f90>
use SetPrecision
implicit none

! Vortex Rebound Problem Methodology From:

! Clercx, H., Breneau, C., "The Normal and Oblique Collision of a Dipole with
! a No-Slip Boundary Condition"
! Journal of Computers and Fluids, (2006) 245-279

! Noting The Co-ordinate Convention Used in the Following

!       Physical                          Computational
!
!       y                                 eta
!       ^             ==> eta(y) ==>      ^
!       |                                 |
!       |                                 |
!       |             ==>  xi(x) ==>      |
!       *------> x                        *------> xi


contains

subroutine Initialise(omega)
implicit none
real(wp) :: omega(:,:)
! In the Lid-Driven Cavity-Flow simulation, there are no defined initial
! conditions - set initial vorticity field to zero.
omega(:,:) = 0.0_wp

RETURN
end subroutine Initialise

subroutine psiBC(psi, TargetGrid)
! Apply Boundary Conditions to the Stream Function field
! The TargetGRID input is an integer and defines the grid on which PoissonSolver will solve:
!
!   1 --> Fine Grid
!   2 --> Upper Intermediate Grid
!   3 --> Lower Intermediate Grid
!   4 --> Coarse Grid
use MeshGen
implicit none

! Inputs
real(wp) :: psi(:,:)
integer  :: TargetGRID

! Pointer to Point to the computational mesh to which psi belongs
type(GRID), pointer :: ptr => null()

! Constant Coeffient to enforce 1st derivative at the lid is 1
real(wp) :: const


call GetMesh(TargetGRID,ptr)

! Apply Boundary Conditions to psi

psi(ptr%nxi,:       ) = 0.0_wp
psi(1      ,:       ) = 0.0_wp  
psi(:      ,1       ) = 0.0_wp

! Along Lid (x,y=end) , u = d psi/dy = 1
#ifdef yGridRefinement
    const = (12.0_wp*ptr%deta)/ptr%DetaDy(ptr%neta)
#endif
#ifndef yGridRefinement
    const = 12.0_wp*ptr%deta
#endif

psi(:      ,ptr%neta) = const


! Nullift Pointer
ptr => null()

RETURN
end subroutine psiBC

subroutine psi2u(psi,u,v)
! Calculate the Velocity Field from the Stream Function by
!
!            d                  d
!       u = --(psi),   v = (-) --(psi) 
!           dy                 dx
!
! And Apply Boundary Conditions
use MeshGen
use Derivatives
implicit none

! Declare Inputs
real(wp) :: psi(1:nx,1:ny)
real(wp) ::   u(1:nx,1:ny)
real(wp) ::   v(1:nx,1:ny)

#if defined(xGridRefinement) .OR. defined(yGridRefinement)
    ! Pointer to Point to Fine Grid for Transformation Metrics
    type(GRID), pointer :: ptr => null()
    
    ! Need Loop Indices
    integer :: i,j
    
    call GetMesh(1,ptr)     ! Point to the Fine Computational Grid

#endif

u =  Df(psi,2,1)
v = -Df(psi,1,1)

#ifdef xGridRefinement
    ! Apply Metrics of Transformation to Detivatives in Computational Domain
    do j = 1,ny
        u(:,j) = u(:,j)*ptr%DetaDy(j)
    end do
#endif

#ifdef yGridRefinement
    ! Apply Metrics of Transformation to Detivatives in Computational Domain
    do i = 1,nx
        v(i,:) = v(i,:)*ptr%DxiDx(i)
    end do
#endif


! Apply No-Slip Boundary Conditions To all 3 walls

u(nx,: ) = 0.0_wp
u(1 ,: ) = 0.0_wp
u(: ,ny) = 1.0_wp       ! Except at the lid were u=1
u(: ,1 ) = 0.0_wp

v(nx,: ) = 0.0_wp
v(1 ,: ) = 0.0_wp
v(: ,ny) = 0.0_wp
v(: ,1 ) = 0.0_wp

#if defined(xGridRefinement) .OR. defined(yGridRefinement)
    ! Nullift Pointer
    ptr => null()
#endif

RETURN
end subroutine psi2u

function u2omega(u,v)
! u2omega computes the vorticity, omega, from the curl of velocity.
! In 2D, it is:
! omega = dv/dx - du/dy
use MeshGen
use Derivatives
implicit none
! Inputs
real(wp) ::       u(1:nx,1:ny)
real(wp) ::       v(1:nx,1:ny)

! Returns
real(wp) :: u2omega(1:nx,1:ny)

! If Grid Refinement is used, we calculate derivatives in computational domain and 
! transform back to physical domain using the fine grid transformation metrics
#ifdef xGridRefinement
    real(wp) :: v_xi (1:nx,1:ny)     ! Derivatives in Computaional Domain
    integer :: i
#endif
#ifdef yGridRefinement
    real(wp) :: u_eta(1:nx,1:ny)     ! Derivatives in Computaional Domain
    integer :: j
#endif

! Pointer to Point to Fine Grid for Transformation Metrics
type(GRID), pointer :: ptr => null()

! Define Recurring Coefficients
real(wp) :: r1_12dxi
real(wp) :: r1_12deta

call GetMesh(1,ptr)

! Constants in BCs
 r1_12dxi = 1.0_wp/(12.0_wp*ptr%dxi)
r1_12deta = 1.0_wp/(12.0_wp*ptr%deta)

! Calculate Vorticity Field
#if .NOT.defined(xGridRefinement) && .NOT.defined(yGridRefinement)
    ! Physical Domain is uniform and no transformation is required
    
    u2omega = Df(v,1,1)-Df(u,2,1)

#elif defined(xGridRefinement) && .NOT.defined(yGridRefinement)
    ! Must transform derivatives w.r.t. xi
    
    ! Compute Derivatives in Computational Space
    v_xi = Df(v,1,1)
    
    ! Map back to Physical space
    do i = 1,nx
        v_xi(i,:) = v_xi(i,:)*ptr%DxiDx(i)
    end do
    
    ! Contruct Vorticity Field
    u2omega = v_xi-Df(u,2,1)

#elif .NOT.defined(xGridRefinement) && defined(yGridRefinement)
    ! Must transform derivatives w.r.t. eta
    
    ! Compute Derivatives in Computational Space
    u_eta = Df(u,2,1)
    
    ! Map back to Physical space
    do j = 1,ny
        u_eta(:,j) = u_eta(:,j)*ptr%Detady(j)
    end do
    
    ! Contruct Vorticity Field
    u2omega = Df(v,1,1)-u_eta
    
#elif defined(xGridRefinement) && defined(yGridRefinement)
    ! Must transform derivatives w.r.t. xi & eta
    
    ! Compute Derivatives in Computational Space
    v_xi = Df(v,1,1)
    u_eta = Df(u,2,1)
    
    ! Map back to Physical space
    do i = 1,nx
        v_xi(i,:) = v_xi(i,:)*ptr%DxiDx(i)
    end do
    do j = 1,ny
        u_eta(:,j) = u_eta(:,j)*ptr%Detady(j)
    end do
    
    ! Contruct Vorticity Field
    u2omega = v_xi-u_eta

#endif


! Boundary Conditions =============================================

! At no slip boundary (x=1,y), omega = dv/dx
u2omega(nx,:) =   (3.0_wp*v(nx-4,:)-16.0_wp*v(nx-3,:) + &
                  36.0_wp*v(nx-2,:)-48.0_wp*v(nx-1,:))*r1_12dxi

#ifdef xGridRefinement
! Transform the derivative into physical space
    u2omega(nx,:) = u2omega(nx,:)*ptr%DxiDx(nx)
#endif

! At no slip boundary (x=-1,y), omega = dv/dx
u2omega(1,:)  =   (48.0_wp*v(2,:)-36.0_wp*v(3,:) + &
                   16.0_wp*v(4,:) -3.0_wp*v(5,:))*r1_12dxi

#ifdef xGridRefinement
! Transform the derivative into physical space
    u2omega(1,:) = u2omega(1,:)*ptr%DxiDx(1)
#endif

! At no slip boundary (x, y=1), omega = -du/dy
u2omega(:,ny) = - (3.0_wp*u(:,ny-4)-16.0_wp*u(:,ny-3) + &
                  36.0_wp*u(:,ny-2)-48.0_wp*u(:,ny-1))*r1_12deta

#ifdef yGridRefinement
! Transform the derivative into physical space
    u2omega(:,ny) = u2omega(:,ny)*ptr%DetaDy(ny)
#endif

! At no slip boundary (x, y=-1), omega = -du/dy
u2omega(:,1)  = - (48.0_wp*u(:,2)-36.0_wp*u(:,3) + &
                   16.0_wp*u(:,4)- 3.0_wp*u(:,5))*r1_12deta

#ifdef yGridRefinement
! Transform the derivative into physical space
    u2omega(:,1) = u2omega(:,1)*ptr%DetaDy(1)
#endif

! At Bottom 2 courners , omega = 0.0
u2omega(1 ,1 ) = 0.0_wp
u2omega(nx,1 ) = 0.0_wp

! Nullift Pointer
ptr => null()

RETURN
end function u2omega

end module CavityFlow
