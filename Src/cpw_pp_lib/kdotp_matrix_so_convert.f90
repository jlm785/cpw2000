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

!>     Converts the matrices with spin-orbit from the perturbation
!>     ordering (4 spin blocks neig*neig) to the spinor representation

       subroutine kdotp_matrix_so_convert(neig,hso,dhsodrk,d2hsodrk2,    &
     &      nder,                                                        &
     &      mxdbnd)

!      written 29 February 2020. JLM
!      copyright  Jose Luis Martins/INESC-MN

!      version 4.95


       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdbnd                     !<  array dimension for number of bands

       integer, intent(in)                ::  neig                       !<  number of wavefunctions

       integer, intent(in)                ::  nder                       !<  order of derivatives to be calculated. kdotp nder =2, oscillator nder = 1.

!      input and output

       complex(REAL64), intent(inout)     ::  hso(2*mxdbnd,2*mxdbnd)     !<  <Psi|H|Psi> without spin-orbit
       complex(REAL64), intent(inout)     :: dhsodrk(2*mxdbnd,2*mxdbnd,3)    !<  d <Psi|H|Psi> d k
       complex(REAL64), intent(inout)   :: d2hsodrk2(2*mxdbnd,2*mxdbnd,3,3)  !<  d^2 <Psi|H|Psi> d k^2

!      local allocatable arrays

       complex(REAL64), allocatable       ::  temp(:,:)

!      counters

       integer    ::  i, j, n, m


       allocate(temp(2*neig,2*neig))


       do i=1,2*neig
       do j=1,2*neig
         temp(j,i) = hso(j,i)
       enddo
       enddo

       do i=1,neig
       do j=1,neig
         hso(2*j-1,2*i-1) = temp(j     ,i     )
         hso(2*j-1,2*i  ) = temp(j     ,i+neig)
         hso(2*j  ,2*i-1) = temp(j+neig,i     )
         hso(2*j  ,2*i  ) = temp(j+neig,i+neig)
       enddo
       enddo


       do n = 1,3
         do i=1,2*neig
         do j=1,2*neig
           temp(j,i) = dhsodrk(j,i,n)
         enddo
         enddo
         do i=1,neig
         do j=1,neig
           dhsodrk(2*j-1,2*i-1,n) = temp(j     ,i     )
           dhsodrk(2*j-1,2*i  ,n) = temp(j     ,i+neig)
           dhsodrk(2*j  ,2*i-1,n) = temp(j+neig,i     )
           dhsodrk(2*j  ,2*i  ,n) = temp(j+neig,i+neig)
         enddo
         enddo
       enddo

       if(nder > 1) then

         do m = 1,3
         do n = 1,3
           do i=1,2*neig
           do j=1,2*neig
             temp(j,i) = d2hsodrk2(j,i,n,m)
           enddo
           enddo
           do i=1,neig
           do j=1,neig
             d2hsodrk2(2*j-1,2*i-1,n,m) = temp(j     ,i     )
             d2hsodrk2(2*j-1,2*i  ,n,m) = temp(j     ,i+neig)
             d2hsodrk2(2*j  ,2*i-1,n,m) = temp(j+neig,i     )
             d2hsodrk2(2*j  ,2*i  ,n,m) = temp(j+neig,i+neig)
           enddo
           enddo
         enddo
         enddo

       endif

       deallocate(temp)

       return

       end subroutine kdotp_matrix_so_convert
