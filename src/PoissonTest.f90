
program test
use SetPrecision
use MeshGen
use Poissonmod
implicit none

type(GRID), pointer :: ptr
real(wp), allocatable :: omega(:,:), psi(:,:), psiA(:,:), psiui(:,:),omegaui(:,:),psili(:,:),omegali(:,:), psic(:,:), omegac(:,:)
real(wp) :: omega0, A, B, U
integer :: i,j, nrec
real(wp) :: err

xLim = (/-1.0_wp, 1.0_wp/)
yLim = (/-1.0_wp, 1.0_wp/)

call Mesh_MultiGrid(xLim,yLim,0,0)

call GetMesh(1,ptr)


allocate(psi(1:ptr%nxi,1:ptr%neta),omega(1:ptr%nxi,1:ptr%neta),psiA(1:ptr%nxi,1:ptr%neta))
omega0 = 1.0_wp
A = 1.0_wp
B = 1.0_wp
U = 0.1_wp

do i = 1,ptr%nxi
    do j = 1,ptr%neta

        omega(i,j) = omega0

        psiA(i,j) = U*y(j)-omega0*((A*x(i)**2.0_wp + &
             B*y(j)**2.0_wp)/(A+B))/2.0_wp

        psi(i,j) = 0.0_wp

    end do
end do

!print*, ptr%xi
call GetMesh(2,ptr)

allocate(psiui(1:ptr%nxi,1:ptr%neta),omegaui(1:ptr%nxi,1:ptr%neta))


do i = 1,ptr%nxi
    do j = 1,ptr%neta

        omegaui(i,j) = omega0

        psiui(i,j) = 0.0_wp

    end do
end do

call GetMesh(3,ptr)

allocate(psili(1:ptr%nxi,1:ptr%neta),omegali(1:ptr%nxi,1:ptr%neta))

do i = 1,ptr%nxi
    do j = 1,ptr%neta

        omegali(i,j) = omega0

        psili(i,j) = 0.0_wp

    end do
end do

call GetMesh(4,ptr)

allocate(psic(1:ptr%nxi,1:ptr%neta),omegac(1:ptr%nxi,1:ptr%neta))

do i = 1,ptr%nxi
    do j = 1,ptr%neta

        omegac(i,j) = omega0

        psic(i,j) = 0.0_wp

    end do
end do

alpha = 1.4_wp
tol = 0.000000001_wp

call PoissonSolver(psic,omegac,4,tol,alpha)

psili = Interp2(4,psic)

call PoissonSolver(psili,omegali,3,tol,alpha)

psiui = Interp2(3,psili)

call PoissonSolver(psiui,omegaui,2,tol,alpha)

psi = Interp2(2,psiui)

call PoissonSolver(psi,omega,1,tol,alpha)

err = MAXVAL(ABS(psi-psiA))

print*, err

!do i = 1,ptr%nxi
!print*, psiA(i,1:ptr%neta)
!print*, psi(i,1:ptr%neta)
!end do

!do i = 1,ptr%nxi
!print*, psi(i,1:ptr%neta)
!end do

call GetMesh(1,ptr)

nrec = (ptr%nxi*ptr%neta)*8

OPEN(UNIT=12, FILE="psi4.dat", ACTION="write", STATUS="replace")
 
do j = 1,ptr%neta
!print*, psi(1:ptr%nxi,j)-psiA(1:ptr%nxi,j)
!print*, psiA(1:ptr%nxi,j)
!print*, '___________________'
  do i=1,ptr%nxi
    WRITE(12,*) x(i), y(j), psi(i,j)-psiA(i,j)
  end do
end do

CLOSE(unit=12)
end program test
