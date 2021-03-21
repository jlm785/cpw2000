!------------------------------------------------------------!
! This file is distributed as part of the cpw2000 code and   !
! under the terms of the GNU General Public License. See the !
! file `LICENSE' in the root directory of the cpw2000        !
! distribution, or http://www.gnu.org/copyleft/gpl.txt       !
!                                                            !
! The webpage of the cpw2000 code is not yet written         !
!                                                            !
! The cpw2000 code is hosted on GitHub:                      !
!                                                            !
! https://github.com/jlm785/cpw2000                          !
!------------------------------------------------------------!

!>     subroutine calculates the real spherical harmonics (multiplied by r**l and a normalization)
!>     and the corresponding derivatives

!>     zylm(x,y,z,l,m) = r**l sqrt(4 pi / 2l+1)   Y_l,0(theta,phi)

!>     rylm(x,y,z,l,m) = r**l sqrt(4 pi / 2l+1) (1/sqrt(2)) [ (-1)**m Y_l,m(theta,phi) +  Y_l,-m(theta,phi) ]  IF  m> 0
!>     rylm(x,y,z,l,0) = r**l sqrt(4 pi / 2l+1)   Y_l,0(theta,phi)   IF  m=0
!>     rylm(x,y,z,l,m) = r**l sqrt(4 pi / 2l+1) (1/sqrt(2)) (i) [ Y_l,m(theta,phi) - (-1)**m Y_l,-m(theta,phi) ]  IF  m< 0

       subroutine ylm_rc(rc,rylm,drylmdrc,d2rylmdrc2,lc16,               &
     &                  zylm,dzylmdrc,d2zylmdrc2,lmax,nder,ifail)

!      Written 28 February 2014. JLM
!      Modified, real, 1 April 2014. JLM
!      Modified, complex and real, 6 April 2014. JLM
!      Modified documentation August 2019.  JLM

!      copyright INESC-MN/Jose Luis Martins

!      version 4.94

       implicit none
       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       real(REAL64), intent(in)           ::  rc(3)                      !<  r in cartesian coordinates 
       integer, intent(in)                ::  lmax                       !<  maximum angular momentum (lmax <= 3)
       integer, intent(in)                ::  nder                       !<  order of derivatives (nder <= 2)
       logical, intent(in)                ::  lc16                       !<  calculates also the complex version (usual spherical harmonics)

!      output

       real(REAL64), intent(out)          ::  rylm(0:3,-3:3)             !<  r**l sqrt(4 pi / 2l+1) Y_lm(theta,phi), REAL VERSION
       real(REAL64), intent(out)          ::  drylmdrc(3,0:3,-3:3)       !<  d rylm / d rc
       real(REAL64), intent(out)          ::  d2rylmdrc2(3,3,0:3,-3:3)   !<  d^2 rylm / d rc^2

       complex(REAL64), intent(out)       ::  zylm(0:3,-3:3)             !<  r**l sqrt(4 pi / 2l+1) Y_lm(theta,phi)
       complex(REAL64), intent(out)       ::  dzylmdrc(3,0:3,-3:3)       !<  d zylm / d rc
       complex(REAL64), intent(out)       ::  d2zylmdrc2(3,3,0:3,-3:3)   !<  d^2 zylm / d rc^2

       integer, intent(out)               ::  ifail                      !<  wrong (not implemented) input lmax or nder

!      local variables

       real(REAL64)    ::  x, y, z                !  r cos(phi)sin(theta), r sin(phi)sin(theta), r cos(theta)
       real(REAL64)    ::  c

!      constants

       real(REAL64), parameter :: ZERO = 0.0_REAL64, UM = 1.0_REAL64

!      counters

       integer         ::  l, m, j, k

!      check input

       ifail = 0
       if(lmax < 0 .or. lmax > 3) ifail = 1
       if(nder < 0 .or. nder > 2) ifail = 1

       x = rc(1)
       y = rc(2)
       z = rc(3)

!      angular functions

       rylm(0, 0) = UM

       if(lmax > 0) then

         rylm(1,-1) = y
         rylm(1, 0) = z
         rylm(1, 1) = x

         if(lmax > 1) then

           rylm(2,-2) = sqrt(3*UM)*x*y
           rylm(2,-1) = sqrt(3*UM)*z*y
           rylm(2, 0) = (UM/2)*(2*z*z-x*x-y*y)
           rylm(2, 1) = sqrt(3*UM)*z*x
           rylm(2, 2) = sqrt((3*UM)/4)*(x*x-y*y)

           if(lmax > 2) then

             rylm(3,-3) = sqrt((5*UM)/8)*(3*x*x-y*y)*y
             rylm(3,-2) = sqrt(15*UM)*z*x*y
             rylm(3,-1) = sqrt((3*UM)/8)*(4*z*z-x*x-y*y)*y
             rylm(3, 0) = (UM/2)*(2*z*z-3*x*x-3*y*y)*z
             rylm(3, 1) = sqrt((3*UM)/8)*(4*z*z-x*x-y*y)*x
             rylm(3, 2) = sqrt((15*UM)/4)*(x*x-y*y)*z
             rylm(3, 3) =-sqrt((5*UM)/8)*(3*y*y-x*x)*x

           endif

         endif

       endif

       if(nder > 0) then

         drylmdrc(1,0,0) = ZERO
         drylmdrc(2,0,0) = ZERO
         drylmdrc(3,0,0) = ZERO

         if(lmax > 0) then

           drylmdrc(1,1,-1) = ZERO
           drylmdrc(2,1,-1) = UM
           drylmdrc(3,1,-1) = ZERO     

           drylmdrc(1,1, 0) = ZERO
           drylmdrc(2,1, 0) = ZERO
           drylmdrc(3,1, 0) = UM

           drylmdrc(1,1, 1) = UM
           drylmdrc(2,1, 1) = ZERO
           drylmdrc(3,1, 1) = ZERO     

           if(lmax > 1) then

             drylmdrc(1,2,-2) = sqrt(3*UM)*y
             drylmdrc(2,2,-2) = sqrt(3*UM)*x
             drylmdrc(3,2,-2) = ZERO

             drylmdrc(1,2,-1) = ZERO
             drylmdrc(2,2,-1) = sqrt(3*UM)*z
             drylmdrc(3,2,-1) = sqrt(3*UM)*y

             drylmdrc(1,2, 0) = -x
             drylmdrc(2,2, 0) = -y
             drylmdrc(3,2, 0) = 2*z

             drylmdrc(1,2, 1) = sqrt(3*UM)*z
             drylmdrc(2,2, 1) = ZERO
             drylmdrc(3,2, 1) = sqrt(3*UM)*x

             drylmdrc(1,2, 2) = sqrt(3*UM)*x
             drylmdrc(2,2, 2) =-sqrt(3*UM)*y
             drylmdrc(3,2, 2) = ZERO

             if(lmax > 2) then

               drylmdrc(1,3,-3) = sqrt((45*UM)/2)*x*y
               drylmdrc(2,3,-3) = sqrt((45*UM)/8)*(x*x-y*y)
               drylmdrc(3,3,-3) = ZERO

               drylmdrc(1,3,-2) = sqrt(15*UM)*z*y
               drylmdrc(2,3,-2) = sqrt(15*UM)*z*x
               drylmdrc(3,3,-2) = sqrt(15*UM)*x*y

               drylmdrc(1,3,-1) =-sqrt((3*UM)/2)*x*y
               drylmdrc(2,3,-1) = sqrt((3*UM)/8)*(4*z*z-x*x-3*y*y)
               drylmdrc(3,3,-1) = sqrt(24*UM)*z*y
               
               drylmdrc(1,3, 0) =-3*x*z
               drylmdrc(2,3, 0) =-3*y*z
               drylmdrc(3,3, 0) = 3*(UM/2)*(2*z*z-x*x-y*y)
 
               drylmdrc(1,3, 1) = sqrt((3*UM)/8)*(4*z*z-3*x*x-y*y)
               drylmdrc(2,3, 1) =-sqrt((3*UM)/2)*(x*y)
               drylmdrc(3,3, 1) = sqrt(24*UM)*z*x

               drylmdrc(1,3, 2) = sqrt(15*UM)*z*x
               drylmdrc(2,3, 2) =-sqrt(15*UM)*z*y
               drylmdrc(3,3, 2) = sqrt((15*UM)/4)*(x*x-y*y)

               drylmdrc(1,3, 3) = sqrt((45*UM)/8)*(x*x-y*y)
               drylmdrc(2,3, 3) =-sqrt((45*UM)/2)*x*y
               drylmdrc(3,3, 3) = ZERO

            endif

           endif

         endif

         if(nder > 1) then

           do j=1,3
           do k=1,3
             d2rylmdrc2(k,j,0,0) = ZERO
           enddo
           enddo

           if(lmax > 0) then

             do j=1,3
             do k=1,3
               d2rylmdrc2(k,j,1,-1) = ZERO
               d2rylmdrc2(k,j,1, 0) = ZERO
               d2rylmdrc2(k,j,1, 1) = ZERO
             enddo
             enddo

             if(lmax > 1) then

               d2rylmdrc2(1,1,2,-2) = ZERO
               d2rylmdrc2(2,1,2,-2) = sqrt(3*UM)
               d2rylmdrc2(3,1,2,-2) = ZERO
               d2rylmdrc2(2,2,2,-2) = ZERO
               d2rylmdrc2(3,2,2,-2) = ZERO
               d2rylmdrc2(3,3,2,-2) = ZERO

               d2rylmdrc2(1,1,2,-1) = ZERO
               d2rylmdrc2(2,1,2,-1) = ZERO
               d2rylmdrc2(3,1,2,-1) = ZERO
               d2rylmdrc2(2,2,2,-1) = ZERO
               d2rylmdrc2(3,2,2,-1) = sqrt(3*UM)
               d2rylmdrc2(3,3,2,-1) = ZERO
               
               d2rylmdrc2(1,1,2, 0) =-UM
               d2rylmdrc2(2,1,2, 0) = ZERO
               d2rylmdrc2(3,1,2, 0) = ZERO
               d2rylmdrc2(2,2,2, 0) =-UM
               d2rylmdrc2(3,2,2, 0) = ZERO
               d2rylmdrc2(3,3,2, 0) = 2*UM

               d2rylmdrc2(1,1,2, 1) = ZERO
               d2rylmdrc2(2,1,2, 1) = ZERO
               d2rylmdrc2(3,1,2, 1) = sqrt(3*UM)
               d2rylmdrc2(2,2,2, 1) = ZERO
               d2rylmdrc2(3,2,2, 1) = ZERO
               d2rylmdrc2(3,3,2, 1) = ZERO

               d2rylmdrc2(1,1,2, 2) = sqrt(3*UM)
               d2rylmdrc2(2,1,2, 2) = ZERO
               d2rylmdrc2(3,1,2, 2) = ZERO
               d2rylmdrc2(2,2,2, 2) =-sqrt(3*UM)
               d2rylmdrc2(3,2,2, 2) = ZERO
               d2rylmdrc2(3,3,2, 2) = ZERO

               do m =-2,2
               do j = 2,3
               do k =1,j-1
                 d2rylmdrc2(k,j,2,m) = d2rylmdrc2(j,k,2,m)
               enddo
               enddo
               enddo 

               if(lmax > 2) then
 
                 d2rylmdrc2(1,1,3,-3) = sqrt((45*UM)/2)*y
                 d2rylmdrc2(2,1,3,-3) = sqrt((45*UM)/2)*x
                 d2rylmdrc2(3,1,3,-3) = ZERO
                 d2rylmdrc2(2,2,3,-3) =-sqrt((45*UM)/2)*y
                 d2rylmdrc2(3,2,3,-3) = ZERO
                 d2rylmdrc2(3,3,3,-3) = ZERO

                 d2rylmdrc2(1,1,3,-2) = ZERO
                 d2rylmdrc2(2,1,3,-2) = sqrt(15*UM)*z
                 d2rylmdrc2(3,1,3,-2) = sqrt(15*UM)*y
                 d2rylmdrc2(2,2,3,-2) = ZERO
                 d2rylmdrc2(3,2,3,-2) = sqrt(15*UM)*x
                 d2rylmdrc2(3,3,3,-2) = ZERO
 
                 d2rylmdrc2(1,1,3,-1) =-sqrt((3*UM)/2)*y
                 d2rylmdrc2(2,1,3,-1) =-sqrt((3*UM)/2)*x
                 d2rylmdrc2(3,1,3,-1) = ZERO
                 d2rylmdrc2(2,2,3,-1) =-sqrt((27*UM)/2)*y
                 d2rylmdrc2(3,2,3,-1) = sqrt(24*UM)*z
                 d2rylmdrc2(3,3,3,-1) = sqrt(24*UM)*y
 
                 d2rylmdrc2(1,1,3, 0) =-3*z
                 d2rylmdrc2(2,1,3, 0) = ZERO
                 d2rylmdrc2(3,1,3, 0) =-3*x
                 d2rylmdrc2(2,2,3, 0) =-3*z
                 d2rylmdrc2(3,2,3, 0) =-3*y
                 d2rylmdrc2(3,3,3, 0) = 6*z
 
                 d2rylmdrc2(1,1,3, 1) =-sqrt((27*UM)/2)*x
                 d2rylmdrc2(2,1,3, 1) =-sqrt((3*UM)/2)*y
                 d2rylmdrc2(3,1,3, 1) = sqrt(24*UM)*z
                 d2rylmdrc2(2,2,3, 1) =-sqrt((3*UM)/2)*x
                 d2rylmdrc2(3,2,3, 1) = ZERO
                 d2rylmdrc2(3,3,3, 1) = sqrt(24*UM)*x
 
                 d2rylmdrc2(1,1,3, 2) = sqrt(15*UM)*z
                 d2rylmdrc2(2,1,3, 2) = ZERO
                 d2rylmdrc2(3,1,3, 2) = sqrt(15*UM)*x
                 d2rylmdrc2(2,2,3, 2) =-sqrt(15*UM)*z
                 d2rylmdrc2(3,2,3, 2) =-sqrt(15*UM)*y
                 d2rylmdrc2(3,3,3, 2) = ZERO
 
                 d2rylmdrc2(1,1,3, 3) = sqrt((45*UM)/2)*x
                 d2rylmdrc2(2,1,3, 3) =-sqrt((45*UM)/2)*y
                 d2rylmdrc2(3,1,3, 3) = ZERO
                 d2rylmdrc2(2,2,3, 3) =-sqrt((45*UM)/2)*x
                 d2rylmdrc2(3,2,3, 3) = ZERO
                 d2rylmdrc2(3,3,3, 3) = ZERO
 
                 do m =-3,3
                 do j = 2,3
                 do k =1,j-1
                   d2rylmdrc2(k,j,3,m) = d2rylmdrc2(j,k,3,m)
                 enddo
                 enddo
                 enddo 

               endif

             endif

           endif

         endif

       endif
       
       if(lc16) then
       
         do l = 0,lmax

           zylm(l,0) = cmplx(rylm(l,0),ZERO,REAL64)

         enddo

         if(lmax > 0) then
           do l = 1,lmax
             c = UM
             do m = 1,l
               c = -c
               zylm(l,-m) = sqrt(UM/2)*cmplx(rylm(l,m),-rylm(l,-m),REAL64)
               zylm(l,m) = c*conjg(zylm(l,-m))

             enddo
           enddo
         endif
         
         if(nder > 0) then
       
           do l = 0,lmax
             do j = 1,3
               dzylmdrc(j,l,0) = cmplx(drylmdrc(j,l,0),ZERO,REAL64)
             enddo
           enddo

           if(lmax > 0) then
             do l = 1,lmax
               c = UM
               do m = 1,l
                 c = -c
                 do j = 1,3
                   dzylmdrc(j,l,-m) = sqrt(UM/2)*                        &
     &                cmplx(drylmdrc(j,l,m),-drylmdrc(j,l,-m),REAL64)
                   dzylmdrc(j,l,m) = c*conjg(dzylmdrc(j,l,-m))
                 enddo
               enddo
             enddo
           endif
         
         endif
         
         if(nder > 1) then
       
           do l = 0,lmax
             do j = 1,3
             do k = 1,3
               d2zylmdrc2(k,j,l,0) = cmplx(d2rylmdrc2(k,j,l,0),          &
     &                                                ZERO,REAL64)
             enddo
             enddo
           enddo

           if(lmax > 0) then
             do l = 1,lmax
               c = UM
               do m = 1,l
                 c = -c
                 do j = 1,3
                 do k = 1,3
                   d2zylmdrc2(k,j,l,-m) = sqrt(UM/2)*cmplx               &
     &             (d2rylmdrc2(k,j,l,m),-d2rylmdrc2(k,j,l,-m),REAL64)
                   d2zylmdrc2(k,j,l,m) = c*conjg(d2zylmdrc2(k,j,l,-m))
                 enddo
                 enddo
               enddo
             enddo
           endif
         
         endif

       endif                                                             !   lc16 == true


       return

       end subroutine ylm_rc
