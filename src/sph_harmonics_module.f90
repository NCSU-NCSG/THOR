module sph_harmonics_module
!*********************************************************************
!
! Module spherical harmonics constructs real spherical harmonic
! expansion coefficient
!
!*********************************************************************
! User derived-type modules

  use types
  use parameter_types
  use vector_types
  use angle_types
  use multindex_types

! Use modules that pertain setting up problem

  use termination_module

  implicit none

contains

  subroutine spherical_harmonics(sord,nang,quad,max_p,Y)
  !*********************************************************************
  !
  ! Subroutine spherical harmonics constructs real spherical harmonic
  ! expansion coefficient
  !
  !*********************************************************************
  ! Pass input parameters

    integer(kind=li), intent(in) :: sord, nang

  ! Pass angular derived type used globally

    type(ordinate), dimension(nang), intent(in) :: &
         quad

  ! Define temporary variables

    integer(kind=li) :: q, octant,l,m,indx
    integer(kind=li), intent(inout) :: max_p
    real(kind=d_t), dimension(nang,8,max_p) :: Y
    real(kind=d_t) :: mu,eta,xi
    type(vector) :: omega

  ! Do even contributions  
    do l=0,sord
      do m=0,l 
        indx=1_li + m + (l+1_li)*l/2_li
        do octant=1,8
           do q=1, nang
             if(octant == 1)then
               mu  = quad(q)%mu%x1 
               eta = quad(q)%mu%x2
               xi  = quad(q)%mu%x3
             elseif(octant == 2)then
               mu  =-quad(q)%mu%x1 
               eta = quad(q)%mu%x2
               xi  = quad(q)%mu%x3
             elseif(octant == 3)then
               mu  =-quad(q)%mu%x1 
               eta =-quad(q)%mu%x2
               xi  = quad(q)%mu%x3
             elseif(octant == 4)then
               mu  = quad(q)%mu%x1 
               eta =-quad(q)%mu%x2
               xi  = quad(q)%mu%x3
             elseif(octant == 5)then
               mu  = quad(q)%mu%x1 
               eta = quad(q)%mu%x2
               xi  = -quad(q)%mu%x3
             elseif(octant == 6)then
               mu  =-quad(q)%mu%x1 
               eta = quad(q)%mu%x2
               xi  = -quad(q)%mu%x3
             elseif(octant == 7)then
               mu  =-quad(q)%mu%x1 
               eta =-quad(q)%mu%x2
               xi  = -quad(q)%mu%x3
             else
               mu  = quad(q)%mu%x1 
               eta =-quad(q)%mu%x2
               xi  = -quad(q)%mu%x3
             end if
             Y(q,octant,indx)=y_e(l,m,mu,eta,xi)
           end do
        end do
      end do
    end do
  ! Odd contributions
    do l=1,sord
      do m=1,l 
        indx= neven + m  + (l-1_li)*l/2_li
        do octant=1,8
           do q=1, nang
             if(octant == 1)then
               mu  = quad(q)%mu%x1 
               eta = quad(q)%mu%x2
               xi  = quad(q)%mu%x3
             elseif(octant == 2)then
               mu  =-quad(q)%mu%x1 
               eta = quad(q)%mu%x2
               xi  = quad(q)%mu%x3
             elseif(octant == 3)then
               mu  =-quad(q)%mu%x1 
               eta =-quad(q)%mu%x2
               xi  = quad(q)%mu%x3
             elseif(octant == 4)then
               mu  = quad(q)%mu%x1 
               eta =-quad(q)%mu%x2
               xi  = quad(q)%mu%x3
             elseif(octant == 5)then
               mu  = quad(q)%mu%x1 
               eta = quad(q)%mu%x2
               xi  = -quad(q)%mu%x3
             elseif(octant == 6)then
               mu  =-quad(q)%mu%x1 
               eta = quad(q)%mu%x2
               xi  = -quad(q)%mu%x3
             elseif(octant == 7)then
               mu  =-quad(q)%mu%x1 
               eta =-quad(q)%mu%x2
               xi  = -quad(q)%mu%x3
             else
               mu  = quad(q)%mu%x1 
               eta =-quad(q)%mu%x2
               xi  = -quad(q)%mu%x3
             end if
             Y(q,octant,indx)=y_o(l,m,mu,eta,xi)
           end do
        end do
      end do
    end do

  end subroutine spherical_harmonics

  function y_e(l,m,mu,eta,xi)
     integer(kind=li), intent(in) :: l,m
     real(kind=d_t), intent(in) :: mu,eta,xi
     real(kind=d_t) :: y_e
     real(kind=d_t) :: phi

!     phi = acos(eta/sqrt(1.0_d_t-mu**2))
     phi = atan2(eta,mu)

     y_e = sqrt(real(2*l+1,d_t)*factorial_d_t(l-m)/factorial_d_t(l+m))*aLegendreP(l,m,xi)*&
           cos(real(m,d_t)*phi)*(-1.0_d_t)**m

  end function

  function y_o(l,m,mu,eta,xi)
     integer(kind=li), intent(in) :: l,m
     real(kind=d_t), intent(in) :: mu,eta,xi
     real(kind=d_t) :: y_o
     real(kind=d_t) :: phi

!     phi = acos(eta/sqrt(1.0_d_t-mu**2))
     phi = atan2(eta,mu)

     y_o = sqrt(real(2*l+1,d_t)*factorial_d_t(l-m)/factorial_d_t(l+m))*aLegendreP(l,m,xi)*&
           sin(real(m,d_t)*phi)*(-1.0_d_t)**m

  end function

  function aLegendreP(l,m,x)
     integer(kind=li), intent(in) :: l,m
     real(kind=d_t), intent(in) :: x
     real(kind=d_t) :: aLegendreP

     if(l == 0_li) then
        if(m == 0_li) then
           aLegendreP=1.0_d_t
        else 
           call stop_thor(32_li)
        end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     else if(l == 1_li) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(m == -1_li) then
           aLegendreP=0.5_d_t*sqrt(1.0_d_t-x*x) 
        else if (m == 0_li) then 
           aLegendreP=x
        else if (m == 1_li) then 
           aLegendreP=-sqrt(1.0_d_t-x*x) 
        else
           call stop_thor(32_li)  
        end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     else if(l == 2_li) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(m == -2_li) then
           aLegendreP=0.125_d_t*(1.0_d_t-x*x)
        else if(m == -1_li) then
           aLegendreP=0.5_d_t*x*sqrt(1.0_d_t-x*x)
        else if(m == 0_li) then
           aLegendreP=0.5_d_t*(3.0_d_t*x*x-1.0_d_t)
        else if(m == 1_li) then
           aLegendreP=-3.0_d_t*x*sqrt(1.0_d_t-x*x)
        else if(m == 2_li) then
           aLegendreP=-3.0_d_t*(x*x-1.0_d_t) 
        else
           call stop_thor(32_li)
        end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     else if(l == 3_li) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(m == -3_li) then
           aLegendreP=0.02083333333333333_d_t*sqrt((1.0_d_t-x*x)**3)   
        else if(m == -2_li) then
           aLegendreP=-0.125_d_t*x*(x*x-1.0_d_t)
        else if(m == -1_li) then
           aLegendreP=0.125_d_t*sqrt(1.0_d_t-x*x)*(5.0_d_t*x*x-1.0_d_t)
        else if(m == 0_li) then
           aLegendreP=0.5_d_t*(-3.0_d_t*x+5.0_d_t*x*x*x)
        else if(m == 1_li) then
           aLegendreP=-1.5_d_t*sqrt(1.0_d_t-x*x)*(5.0_d_t*x*x-1.0_d_t)
        else if(m == 2_li) then
           aLegendreP=-15.0_d_t*x*(x*x-1.0_d_t)
        else if(m == 3_li) then
           aLegendreP=-15.0_d_t**sqrt((1.0_d_t-x*x)**3)
        else
           call stop_thor(32_li)
        end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     else if(l == 4_li) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(m == -4_li) then
           aLegendreP=0.002604166666666667_d_t*(-1.0_d_t+x*x)**2
        else if(m == -3_li) then
           aLegendreP=0.02083333333333333_d_t*x*sqrt((1.0_d_t-x*x)**3)
        else if(m == -2_li) then
           aLegendreP=-0.02083333333333333_d_t*(-1.0_d_t+x*x)*(-1.0_d_t+7.0_d_t*x*x)
        else if(m == -1_li) then
           aLegendreP=0.125_d_t*sqrt(1.0_d_t-x*x)*(-3.0_d_t*x+7.0_d_t*x*x*x) 
        else if(m == 0_li) then
           aLegendreP=0.125_d_t*(3.0_d_t-30.0_d_t*x*x+35.0_d_t*x*x*x*x)
        else if(m == 1_li) then
           aLegendreP=-2.5_d_t*sqrt(1.0_d_t-x*x)*(-3.0_d_t*x+7.0_d_t*x*x*x)
        else if(m == 2_li) then
           aLegendreP=-7.5_d_t*(-1.0_d_t+x*x)*(-1.0_d_t+7.0_d_t*x*x)
        else if(m == 3_li) then
           aLegendreP=-105.0_d_t*x*sqrt((1.0_d_t-x*x)**3)
        else if(m == 4_li) then
           aLegendreP=105.0_d_t*(-1.0_d_t+x*x)**2
        else
           call stop_thor(32_li)
        end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     else if(l == 5_li) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(m == -5_li) then
           aLegendreP=0.0002604166666666667_d_t*sqrt((1.0_d_t-x*x)**5)
        else if(m == -4_li) then
           aLegendreP=0.002604166666666667_d_t*x*(-1.0_d_t+x*x)**2
        else if(m == -3_li) then
           aLegendreP=-0.002604166666666667_d_t*sqrt(1.0_d_t-x*x)*(-1.0_d_t+x*x)*(-1.0_d_t+9.0_d_t*x*x)
        else if(m == -2_li) then
           aLegendreP=-0.0625_d_t*(-1.0_d_t+x*x)*(-x+3.0_d_t*x*x*x)
        else if(m == -1_li) then
           aLegendreP=0.0625_d_t*sqrt(1.0_d_t-x*x)*(1.0_d_t-14.0_d_t*x*x+21.0_d_t*x*x*x*x)
        else if(m == 0_li) then
           aLegendreP=0.125_d_t*(15.0_d_t*x-70.0_d_t*x*x*x+63.0_d_t*x*x*x*x*x)
        else if(m == 1_li) then
           aLegendreP=-1.875_d_t*sqrt(1.0_d_t-x*x)*(1.0_d_t-14.0_d_t*x*x+21.0_d_t*x*x*x*x)
        else if(m == 2_li) then
           aLegendreP=-52.5_d_t*(-1.0_d_t+x*x)*(-x+3.0_d_t*x*x*x)
        else if(m == 3_li) then
           aLegendreP=52.5_d_t*sqrt(1.0_d_t-x*x)*(-1.0_d_t+x*x)*(-1.0_d_t+9.0_d_t*x*x)
        else if(m == 4_li) then
           aLegendreP=945.0_d_t*x*(-1.0_d_t+x*x)**2
        else if(m == 5_li) then
           aLegendreP=-945.0_d_t*sqrt((1.0_d_t-x*x)**5) 
        else
           call stop_thor(32_li)
        end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     else  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call stop_thor(33_li)
     end if
 
  end function

  function LegendreP(l,x)
    integer(kind=li), intent(in) :: l
    real(kind=d_t), intent(in) :: x
    real(kind=d_t) :: LegendreP
    
    if(l == 0_li)then
       LegendreP=1
    elseif(l == 1_li)then
       LegendreP=x
    elseif(l == 2_li)then
       LegendreP=0.5_d_t*(3.0_d_t*x*x-1)
    elseif(l == 3_li)then
       LegendreP=0.5_d_t*(5.0_d_t*x*x*x-3.0_d_t*x)
    elseif(l == 4_li)then
       LegendreP=0.125_d_t*(35.0_d_t*x*x*x*x-30.0_d_t*x*x+3.0_d_t)
    elseif(l == 5_li)then
       LegendreP=0.125_d_t*(63.0_d_t*x*x*x*x*x-70.0_d_t*x*x*x+&
            15.0_d_t*x)
    elseif(l == 6_li)then
       LegendreP=0.0625_d_t*(231.0_d_t*x*x*x*x*x*x-315.0_d_t*x*x*x*x+&
            105.0_d_t*x*x-5.0_d_t)
    elseif(l == 7_li)then
       LegendreP=0.0625_d_t*(429.0_d_t*x*x*x*x*x*x*x-&
            693.0_d_t*x*x*x*x*x+315.0_d_t*x*x*x-35.0_d_t*x)
    else
      call stop_thor(23_li) 
    end if
    
  end function LegendreP

end module sph_harmonics_module

