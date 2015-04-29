! Module for Marching Vorticity Equation in Time
! using Runge-Kutta 4
! E. Higham Feb 2014
! Final Year Project - 2D Laminar Flow Solver
module RKMod
#include <src/PreProcVars.f90>
use SetPrecision
use MeshGen
implicit none

!===============================================================================!
!                                TIME DOMAIN                                    !
!===============================================================================!
! Declare Variables
real(wp) :: tLim(1:2)               ! Start and End Limits of Time Integration
integer  :: nt                      ! Number of Discrete Time Points
real(wp) :: dt                      ! Time Step
real(wp) :: t                       ! Time Variable

contains

subroutine time(tLim,nt)
implicit none
! time establishes the discrete time domain in which the Vorticity equation will
! be solved.
! tLim is a real(wp) vector of limits for the time integration scheme of the form:
!
!       tLim(/t0, tEnd/), where t0      is the initial time and 
!                               tEnd    is the end time of the simulation
! time takes the imput variables

real(wp) :: tLim(1:2)               ! Start and End Limits of Time Integration
integer  :: nt                      ! Number of Discrete Time Points

dt = (tLim(2)-tLim(1))/(REAL(nt-1,wp))
t = tlim(1)

RETURN
end subroutine time


function RHSVorticity(omega,u,v,Re)
! RHSVorticity evaluates the RHS of the 2D Incompressible Vorticity equation:
!
!   d               d             d              -1 [  d2            d2        ]
!  --(omega) = - u --(omega) - v --(omega) + (Re)   [ ---(omega) +  ---(omega) ]
!  dt              dx            dy                 [ dx2           dy2        ]
!
! For use with Runge-Kutta scheme to solve in time. Noting the Cartesian
! Co-ordinate System Convention:

!       Physical                          Computational
!
!       y                                 eta
!       ^             ==> eta(y) ==>      ^
!       |                                 |
!       |                                 |
!       |             ==>  xi(x) ==>      |
!       *------> x                        *------> xi 

! Spacial Deriviative are approximated using 4th orther finite difference
! schemes given in module Derivatives
use Derivatives
implicit none

! RHSVorticity takes the Input variables
real(wp) :: omega(1:nx,1:ny)          ! 2D Vorticity Field
real(wp) ::     u(1:nx,1:ny)          ! 2D x-Component of Velocity Array
real(wp) ::     v(1:nx,1:ny)          ! 2D y-Component of Velocity Array
real(wp) ::    Re                     ! Reynolds Number

! And outputs the time derivative evaluated at that instant
real(wp) :: RHSVorticity(1:nx,1:ny)

! Need to evaluate and temporarily store spacial derivatives
! 1st Derivative of omega...
real(wp) ::  omega_x(1:nx,1:ny)  ! w.r.t. x
real(wp) ::  omega_y(1:nx,1:ny)  ! w.r.t. y
! 2nd Derivative of omega...
real(wp) :: omega_xx(1:nx,1:ny)  ! w.r.t. x
real(wp) :: omega_yy(1:nx,1:ny)  ! w.r.t. y

! Need loop indices
integer  :: i,j

! For efficiency, define the reciprocal of the Re
real(wp) :: r1_Re


#if defined(xGridRefinement).OR.defined(yGridRefinement)
! Need a pointer to point to fine grid transformation metrics
type(GRID), pointer :: ptr

call GetMesh(1,ptr)

#endif

! Evaluate 1st Derivatives
omega_x  =  Df(omega,1,1)
omega_y  =  Df(omega,2,1)

! Evaluate 2nd Derivatives
omega_xx = D2f(omega,1,1)
omega_yy = D2f(omega,2,1)

! Reciprocal of Re
r1_Re = 1.0_wp/Re

! Combine Terms

do j = 1,ny
    do i = 1,nx

#if .NOT.defined(xGridRefinement) && .NOT.defined(yGridRefinement)

    RHSVorticity(i,j) =  -u(i,j)*omega_x(i,j) - v(i,j)*omega_y(i,j) + &
                         r1_Re*(omega_xx(i,j) +       omega_yy(i,j))

#elif defined(xGridRefinement) && .NOT.defined(yGridRefinement)

    RHSVorticity(i,j) =  -u(i,j)*omega_x (i,j)*ptr%DxiDx  (i) - &
                          v(i,j)*omega_y (i,j)                + &
                          r1_Re*(omega_x (i,j)*ptr%D2xiDx2(i) + &
                                 omega_xx(i,j)*ptr%DxiDxsq(i) + &
                                 omega_yy(i,j))

#elif .NOT.defined(xGridRefinement) && defined(yGridRefinement)

    RHSVorticity(i,j) =  -u(i,j)*omega_x (i,j)                 - &
                          v(i,j)*omega_y (i,j)*ptr%DetaDy  (j) + &
                          r1_Re*(omega_xx(i,j)                 + &
                                 omega_y (i,j)*ptr%D2etaDy2(j) + &
                                 omega_yy(i,j)*ptr%DetaDysq(j))

#elif defined(xGridRefinement) && defined(yGridRefinement)

    RHSVorticity(i,j) =  -u(i,j)*omega_x (i,j)*ptr%DxiDx   (i) - &
                          v(i,j)*omega_y (i,j)*ptr%DetaDy  (j) + &
                          r1_Re*(omega_x (i,j)*ptr%D2xiDx2 (i) + &
                                 omega_xx(i,j)*ptr%DxiDxsq (i) + &
                                 omega_y (i,j)*ptr%D2etaDy2(j) + &
                                 omega_yy(i,j)*ptr%DetaDysq(j))

#endif

    end do
end do

RETURN 
end function RHSVorticity

subroutine RKsolver(omega,u,v,psi,Re)
! RKsolver marches the 2D incompressible virticity equation in time by
! step dt using a 4th order Runge-Kutta numerical method.
! The interstitial steps of RK4 are:
! 
! omega(t+dt) = omega(t) + dt*(k1 + 2k2 + 2k3 + k4)/6

! where     k1 = RHSVorticity(omega        ,u,v,Re,t     )
!           k2 = RHSVorticity(omega+dt*k1/2,u,v,Re,t+dt/2) 
!           k3 = RHSVorticity(omega+dt*k2/2,u,v,Re,t+dt/2)
!           k4 = RHSVorticity(omega+dt*k3  ,u,v,Re,t+dt  )

! Note u,v must be updated at each interstitial step
use SimulationVars
use Poissonmod
implicit none

! RKsolver takes the Input variables
real(wp) :: omega(1:nx,1:ny)          ! 2D Vorticity Field
real(wp) ::     u(1:nx,1:ny)          ! 2D x-Component of Velocity Array
real(wp) ::     v(1:nx,1:ny)          ! 2D y-Component of Velocity Array
real(wp) ::   psi(1:nx,1:ny)          ! 2D Steam Function Array
real(wp) ::    Re                     ! Reynolds Number

! And requires the following iterstitial steps
real(wp) :: k1(1:nx,1:ny)
real(wp) :: k2(1:nx,1:ny)
real(wp) :: k3(1:nx,1:ny)
real(wp) :: k4(1:nx,1:ny)

! With a Temporary Vorticity Variable at Each Step
real(wp) :: omegaTemp(1:nx,1:ny)

! Repeated Coefficient
real(wp) :: r1_6 

r1_6 = 1.0_wp/6.0_wp
!==============================================================================!
!                              RK4 Method                                      !
!==============================================================================!

!******************************************************************************!
! Interstitial Step k1

            k1 = RHSVorticity(omega,u,v,Re)
omegaTemp(:,:) = omega(:,:) + 0.5_wp*dt*k1(:,:)

!*****************************************************************************!
! Interstitial Step k2
call MG(psi,omegaTemp)
call psi2u(psi,u,v)

     omegaTemp = u2omega(u,v)
            k2 = RHSVorticity(omegaTemp,u,v,Re)
omegaTemp(:,:) = omega(:,:) + 0.5_wp*dt*k2(:,:)

!*****************************************************************************!
! Interstitial Step k3
call MG(psi,omegaTemp)
call psi2u(psi,u,v)

     omegaTemp = u2omega(u,v)
            k3 = RHSVorticity(omegaTemp,u,v,Re)
omegaTemp(:,:) = omega(:,:) + dt*k3(:,:)

!*****************************************************************************!
! Interstitial Step k4
call MG(psi,omegaTemp)
call psi2u(psi,u,v)

     omegaTemp = u2omega(u,v)
            k4 = RHSVorticity(omegaTemp,u,v,Re)

!*****************************************************************************!
! Evaluate omega(t+dt)
omegaTemp(:,:) = r1_6*dt*(k1(:,:) + 2.0_wp*k2(:,:) + 2.0_wp*k3(:,:) + k4(:,:))
omega(:,:) = omega(:,:) + omegaTemp(:,:)

call MG(psi,omega)
call psi2u(psi,u,v)

omega = u2omega(u,v)

RETURN
end subroutine RKsolver

subroutine RK4_CFL(u,v)
! RK_CFL performs CFL analyis given max streamwise and vertical velocities
! uMax and vMax.
! CFL Analysis from Moin, P., "Fundamentals of Engineering Numerical Analysis"

! Declare Inputs
real(wp) ::   u(1:nx,1:ny)
real(wp) ::   v(1:nx,1:ny)
real(wp) :: CFL(1:nx,1:ny)

! Working Variables
real(wp) :: uMax, vMax
real(wp) :: CFLMax, CFL_RK4

CFL_RK4 = 2.83_wp

! Caclulate the CFL number

CFL = u(:,:)*dt/grid1%dxi + v(:,:)*dt/grid1%deta
CFLMax = MAXVAL(ABS(CFL))

! Ensure that CFL analytis is satisfied
if (CFLMax>CFL_RK4) then
!    dt = 0.9_wp*CFLMax/(uMax/grid1%dxi+vMax/grid1%deta)
!    write(*,*) 'CFL NOT SATISFIED, dt = ', dt
!    CFL = uMax*dt/grid1%dxi + vMax*dt/grid1%deta
    print*, CFLMax
    write(*,*) 'CFL NOT SATISFIED, ABORT'
    stop
end if

RETURN
end subroutine RK4_CFL

end module RKMod