! Module for Numerical Differentiation
! E. Higham Mar 2014
! Final Year Project - 2D Laminar Flow Solver
module Derivatives
use MeshGen
use SetPrecision
implicit none

contains
function Df(f,dim,tgtGrid)
    ! Df approximates the 1st derivative w.r.t. dim of an uniformly spaced 2D
    ! array to 4th order accuracy on grid tgtGrid. tgtGrid is optional and by default
    ! is the fine grid, grid 1

!   Dimension dim must be an integer and can take the values:

!   dim =   1 --> Calculate the x-derivative
!           2 --> Calculate the y-derivative

! ! The tgtGrid input is an integer and defines the grid on which the derivative will be calculated:
!
!   1 --> Grid 1
!   2 --> Grid 2
!   3 --> Grid 3
!   4 --> Grid 4
!   5 --> Grid 5
!   6 --> Grid 6
        
    implicit none

    real(wp), allocatable :: Df(:,:)
    real(wp)              ::  f(:,:)       ! array to be differentiated numerically
    real(wp)              :: a(4), b(5), c(5)    ! Central, Off-Center and Single Sided Difference Coefficients
    integer               :: dim
    integer               :: i,j                 ! Loop Indicies 
    integer               :: tgtGrid

    type(GRID), pointer :: ptr                   ! Pointer to point to grid type

    call GetMesh(tgtGrid,ptr)

    allocate(Df(1:ptr%nxi,1:ptr%neta))

    ! Finite Difference Coeffcients - to double precision
    a = (/1.0_wp,-8.0_wp,8.0_wp,-1.0_wp/)
    b = (/-3.0_wp,-10.0_wp,18.0_wp,-6.0_wp,1.0_wp/)
    c = (/-25.0_wp, 48.0_wp, -36.0_wp, 16.0_wp, -3.0_wp/)

    SELECT CASE (dim)

    CASE (1)            ! Calculate the x-derivative

        ! Central Difference Scheme
        do i = 3,ptr%nxi-2
            Df(i,:) =  a(1)*f(i-2,:)+a(2)*f(i-1,:)+a(3)*f(i+1,:)+a(4)*f(i+2,:)
        end do

        i = ptr%nxi
        ! Off-Centre Scheme
        Df(2  ,:) =   b(1)*f(1,:)+b(2)*f(2  ,:)+b(3)*f(3  ,:)+b(4)*f(4  ,:)+b(5)*f(5  ,:)
        Df(i-1,:) = -(b(1)*f(i,:)+b(2)*f(i-1,:)+b(3)*f(i-2,:)+b(4)*f(i-3,:)+b(5)*f(i-4,:))

        ! Single Sided Scheme
        Df(1,:) =   c(1)*f(1,:)+c(2)*f(2  ,:)+c(3)*f(3  ,:)+c(4)*f(4  ,:)+c(5)*f(5  ,:)
        Df(i,:) = -(c(1)*f(i,:)+c(2)*f(i-1,:)+c(3)*f(i-2,:)+c(4)*f(i-3,:)+c(5)*f(i-4,:))

    Df(:,:) = Df(:,:)/(12.0_wp*ptr%dxi)

    CASE (2)       ! Caluclate the y-derivative

        ! Central Difference Scheme
        do j = 3,ptr%neta-2
            Df(:,j) =  a(1)*f(:,j-2)+a(2)*f(:,j-1)+a(3)*f(:,j+1)+a(4)*f(:,j+2)
        end do

        j = ptr%neta
        ! Off-Centre Scheme
        Df(:,2  ) =   b(1)*f(:,1)+b(2)*f(:,2  )+b(3)*f(:,3  )+b(4)*f(:,4  )+b(5)*f(:,5  )
        Df(:,j-1) = -(b(1)*f(:,j)+b(2)*f(:,j-1)+b(3)*f(:,j-2)+b(4)*f(:,j-3)+b(5)*f(:,j-4))

        ! Single Sided Scheme
        Df(:,1  ) =   c(1)*f(:,1)+c(2)*f(:,2  )+c(3)*f(:,3  )+c(4)*f(:,4  )+c(5)*f(:,5  )
        Df(:,j  ) = -(c(1)*f(:,j)+c(2)*f(:,j-1)+c(3)*f(:,j-2)+c(4)*f(:,j-3)+c(5)*f(:,j-4))

    Df(:,:) = Df(:,:)/(12.0_wp*ptr%deta)

    END SELECT

    ! Nullify Pointer
    ptr => null() 

    RETURN
end function Df

function D2f(f,dim,tgtGrid)
!   D2f approximates the 2nd derivative w.r.t. dim of an uniformly spaced 2D
!   array to 4th order accuracy on grid tgtGrid.

!   Dimension dim must be an integer and can take the values:

!   dim =   1 --> Calculate the x-derivative
!           2 --> Calculate the y-derivative

! ! The tgtGrid input is an integer and defines the grid on which the derivative will be calculated:
!
!   1 --> Grid 1
!   2 --> Grid 2
!   3 --> Grid 3
!   4 --> Grid 4
!   5 --> Grid 5
!   6 --> Grid 6

    implicit none

    real(wp), allocatable :: D2f(:,:)
    real(wp)              ::   f(:,:)            ! array to be differentiated numerically
    real(wp)              :: a(5), b(6), c(6)    ! Central, Off-Center and Single Sided Difference Coefficients
    integer               :: i,j                 ! Array Indicies
    integer               :: dim
    integer               :: tgtGrid

    type(GRID), pointer   :: ptr                 ! Pointer to point to grid type

    ! Finite Difference Coeffcients - to double precision
    a = (/-1.0_wp,   16.0_wp, -30.0_wp,   16.0_wp, -1.0_wp          /)
    b = (/10.0_wp,  -15.0_wp,  -4.0_wp,   14.0_wp, -6.0_wp,   1.0_wp/)
    c = (/45.0_wp, -154.0_wp, 214.0_wp, -156.0_wp, 61.0_wp, -10.0_wp/)

    ! Get the grid to caluclate the derivative on
    call GetMesh(tgtGrid,ptr)

    allocate(D2f(1:ptr%nxi,1:ptr%neta))

    SELECT CASE (dim)

    CASE (1)            ! Calculate the x-derivative

            ! Central Difference Scheme
        do i = 3,ptr%nxi-2
            D2f(i,:) =  a(1)*f(i-2,:)+a(2)*f(i-1,:)+a(3)*f(i,:)+a(4)*f(i+1,:)+a(5)*f(i+2,:)
        end do

            i = ptr%nxi
            ! Off-Centered Scheme
        D2f(2  ,:) =  b(1)*f(1  ,:)+b(2)*f(2  ,:)+b(3)*f(3  ,:) + &
                      b(4)*f(4  ,:)+b(5)*f(5  ,:)+b(6)*f(6  ,:)
        D2f(i-1,:) =  b(1)*f(i  ,:)+b(2)*f(i-1,:)+b(3)*f(i-2,:) + &
                      b(4)*f(i-3,:)+b(5)*f(i-4,:)+b(6)*f(i-5,:)

        ! Single-Sided Scheme
        D2f(1  ,:) = c(1)*f(1  ,:)+c(2)*f(2  ,:)+c(3)*f(3  ,:) + &
                     c(4)*f(4  ,:)+c(5)*f(5  ,:)+c(6)*f(6  ,:)
        D2f(i  ,:) = c(1)*f(i  ,:)+c(2)*f(i-1,:)+c(3)*f(i-2,:) + &
                     c(4)*f(i-3,:)+c(5)*f(i-4,:)+c(6)*f(i-5,:)


        D2f(:,:) = D2f(:,:)/((12.0_wp)*(ptr%dxi*ptr%dxi))

    CASE (2)       ! Calculate the y-derivative

        ! Central Difference Scheme
        do j = 3,ptr%neta-2
            D2f(:,j) =  a(1)*f(:,j-2)+a(2)*f(:,j-1)+a(3)*f(:,j)+a(4)*f(:,j+1)+a(5)*f(:,j+2)
        end do

            j = ptr%neta
        ! Off-Centered Scheme
        D2f(:,2  ) =   b(1)*f(:,1  )+b(2)*f(:,2  )+b(3)*f(:,3  ) + &
                       b(4)*f(:,4  )+b(5)*f(:,5  )+b(6)*f(:,6  )
        D2f(:,j-1) =   b(1)*f(:,j  )+b(2)*f(:,j-1)+b(3)*f(:,j-2) + &
                       b(4)*f(:,j-3)+b(5)*f(:,j-4)+b(6)*f(:,j-5)

        ! Single-Sided Scheme
        D2f(:,1  ) =  c(1)*f(:,1  )+c(2)*f(:,2  )+c(3)*f(:,3  ) + &
                      c(4)*f(:,4  )+c(5)*f(:,5  )+c(6)*f(:,6  )
        D2f(:  ,j) =  c(1)*f(:,j  )+c(2)*f(:,j-1)+c(3)*f(:,j-2) + &
                      c(4)*f(:,j-3)+c(5)*f(:,j-4)+c(6)*f(:,j-5)

        D2f(:,:) = D2f(:,:)/((12.0_wp)*(ptr%deta*ptr%deta))

    END SELECT
    
    ! Nullify Pointer
    ptr => null() 

    RETURN
end function D2f

end module Derivatives