module SetPrecision
    implicit none

    ! Define Precision Suffixes
integer, parameter :: sp = kind(0.0E0)          ! Single Precision
integer, parameter :: dp = kind(0.0D0)          ! Double Precision
!integer, parameter :: qp = kind(0.0Q0)          ! Quadruple Precision

! Working Precision
integer, parameter :: wp = dp   

end module SetPrecision