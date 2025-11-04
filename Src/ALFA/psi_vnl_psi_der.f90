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

!>     calculates the matrix <Psi|V_NL|Psi> for a separable non-local
!>     pseudopotential V_NL for neig wavevectors, and their derivatives
!>     with respect to k-vector. complex version

       subroutine psi_vnl_psi_der(mtxd,neig,nanl,psi,nder,               &
     & vnl,dvnldrk,d2vnldrk2,                                            &
     & anlga,xnlkb,danlgadrk,d2anlgadrk2,                                &
     & mxddim,mxdbnd,mxdanl)

!      written June 2012. jlm
!      Modified 7 January 2014, style. jlm
!      Modified (ladd) 6 February 2014. JLM
!      Modified, documentation, 19 February 2020. JLM
!      copyright INESC-MN/Jose Luis Martins

!      version 4.95

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxddim                     !<  array dimension of plane-waves
       integer, intent(in)                ::  mxdanl                     !<  array dimension of number of projectors
       integer, intent(in)                ::  mxdbnd                     !<  array dimension for number of bands

       integer, intent(in)                ::  mtxd                       !<  wavefunction dimension
       integer, intent(in)                ::  neig                       !<  wavefunction dimension

       integer, intent(in)                ::  nder                       !<  order of derivatives to be calculated. kdotp nder =2, oscillator nder = 1.

       integer, intent(in)                ::  nanl                       !<  number of projectors
       complex(REAL64), intent(in)        ::  anlga(mxddim,mxdanl)       !<  Kleinman-Bylander projectors
       complex(REAL64), intent(in)        ::  danlgadrk(mxddim,mxdanl,3) !<  d anlga / d rkpt
       complex(REAL64), intent(in)    ::  d2anlgadrk2(mxddim,mxdanl,3,3) !<  d^2 anlga / d rkpt^2
       real(REAL64), intent(in)           ::  xnlkb(mxdanl)              !<  Kleinman-Bylander normalization

       complex(REAL64), intent(in)        ::  psi(mxddim,mxdbnd)         !<  wavevectors

!      output

       complex(REAL64), intent(out)       ::  vnl(mxdbnd,mxdbnd)         !<  <Psi|V_NL|Psi>
       complex(REAL64), intent(out)       ::  dvnldrk(mxdbnd,mxdbnd,3)   !<  d <Psi|V_NL|Psi> d k
       complex(REAL64), intent(out)    ::  d2vnldrk2(mxdbnd,mxdbnd,3,3)  !<  d^2 <Psi|V_NL|Psi> d k^2

!      local variables

       complex(REAL64), allocatable       ::  xdhd(:,:)
       complex(REAL64), allocatable       ::  dhd(:,:)
       complex(REAL64), allocatable       ::  ddhddrk(:,:,:)
       complex(REAL64), allocatable       ::  xddhddrk(:,:,:)
       complex(REAL64), allocatable       ::  d2dhddrk2(:,:,:,:)

!      constants

       real(REAL64), parameter :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
       complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
       complex(REAL64), parameter  ::  C_UM = cmplx(UM,ZERO,REAL64)

!      counters

       integer    ::  i,j,m,n

       if(nanl > 0) then

         allocate(dhd(nanl,neig))
         allocate(xdhd(nanl,neig))
         allocate(ddhddrk(nanl,neig,3))

         if(nder > 1) then
           allocate(xddhddrk(nanl,neig,3))
           allocate(d2dhddrk2(nanl,neig,3,3))
         endif

!        dhd = < anl | psi >

         call zgemm('c','n',nanl,neig,mtxd,C_UM,anlga,mxddim,psi,        &
     &                mxddim,C_ZERO,dhd,nanl)

         do j = 1,3

           call zgemm('c','n',nanl,neig,mtxd,C_UM,danlgadrk(1,1,j),      &
     &                mxddim,psi,mxddim,C_ZERO,ddhddrk(1,1,j),nanl)

         enddo

         if(nder > 1) then
           do i = 1,3
           do j = 1,3

             call zgemm('c','n',nanl,neig,mtxd,C_UM,d2anlgadrk2(1,1,j,i),  &
     &                mxddim,psi,mxddim,C_ZERO,d2dhddrk2(1,1,j,i),nanl)

           enddo
           enddo
         endif

!        xdhd := Diag(xnl) dhd

         do n = 1,neig
           do m = 1,nanl
             xdhd(m,n) = xnlkb(m)*dhd(m,n)
           enddo
         enddo

         if(nder > 1) then
           do j = 1,3
             do n = 1,neig
               do m = 1,nanl
                 xddhddrk(m,n,j) = xnlkb(m)*ddhddrk(m,n,j)
               enddo
             enddo
           enddo
         endif

!        vnl := < psi | anl > Diag(xnl) < anl | psi > 

         call zgemm('c','n',neig,neig,nanl,C_UM,dhd,nanl,xdhd,           &
     &                 nanl,C_ZERO,vnl,mxdbnd)

         do j = 1,3

           call zgemm('c','n',neig,neig,nanl,C_UM,ddhddrk(1,1,j),        &
     &                 nanl,xdhd,nanl,C_ZERO,dvnldrk(1,1,j),mxdbnd)

           call zgemm('c','n',neig,neig,nanl,C_UM,xdhd,nanl,             &
     &                 ddhddrk(1,1,j),nanl,C_UM,dvnldrk(1,1,j),mxdbnd)

         enddo

         if(nder > 1) then
           do i = 1,3
           do j = 1,3

             call zgemm('c','n',neig,neig,nanl,C_UM,d2dhddrk2(1,1,j,i),  &
     &                  nanl,xdhd,nanl,C_ZERO,                           &
     &                  d2vnldrk2(1,1,j,i),mxdbnd)

             call zgemm('c','n',neig,neig,nanl,C_UM,xdhd,                &
     &                  nanl,d2dhddrk2(1,1,j,i),nanl,C_UM,               &
     &                  d2vnldrk2(1,1,j,i),mxdbnd)

             call zgemm('c','n',neig,neig,nanl,2*C_UM,xddhddrk(1,1,j),   &
     &                  nanl,ddhddrk(1,1,i),nanl,C_UM,                   &
     &                  d2vnldrk2(1,1,j,i),mxdbnd)
           enddo
           enddo
         endif


         deallocate(dhd)
         deallocate(xdhd)
         deallocate(ddhddrk)

         if(nder > 1) then
           deallocate(xddhddrk)
           deallocate(d2dhddrk2)
         endif

       else

         do n = 1,neig
         do m = 1,neig
           vnl(m,n) = C_ZERO
         enddo
         enddo

         do j = 1,3
           do n = 1,neig
           do m = 1,neig
             dvnldrk(m,n,j) = C_ZERO
           enddo
           enddo
         enddo

         do i = 1,3
         do j = 1,3
           do n = 1,neig
           do m = 1,neig
             d2vnldrk2(m,n,j,i) = C_ZERO
           enddo
           enddo
         enddo
         enddo

       endif

       return

       end subroutine psi_vnl_psi_der
