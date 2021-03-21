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
!>     spin wave-functions are the same
!>     appropriate for spin-perturbation

       subroutine psi_vnl_sp_psi_der(mtxd,neig,nanlso,psi,nder,          &
     &   vnl,dvnldrk,d2vnldrk2,                                          &
     &   anlp,anlm,xnlkb,danlpdrk,danlmdrk,d2anlpdrk2,d2anlmdrk2,        &
     &   mxddim,mxdbnd,mxdaso)

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
       integer, intent(in)                ::  mxdaso                     !<  array dimension of number of projectors
       integer, intent(in)                ::  mxdbnd                     !<  array dimension for number of bands

       integer, intent(in)                ::  mtxd                       !<  wavefunction dimension
       integer, intent(in)                ::  neig                       !<  wavefunction dimension

       complex(REAL64), intent(in)        ::  psi(mxddim,mxdbnd)         !<  wavevectors

       integer, intent(in)                ::  nder                       !<  order of derivatives to be calculated. kdotp nder =2, oscillator nder = 1.

       integer, intent(in)                ::  nanlso                     !<  number of projectors
       complex(REAL64), intent(in)        ::  anlp(mxddim,mxdaso)        !<  Kleinman-Bylander projectors
       complex(REAL64), intent(in)        ::  anlm(mxddim,mxdaso)        !<  Kleinman-Bylander projectors
       complex(REAL64), intent(in)        ::  danlpdrk(mxddim,mxdaso,3)  !<  d anlga / d rkpt
       complex(REAL64), intent(in)        ::  danlmdrk(mxddim,mxdaso,3)  !<  d anlga / d rkpt
       complex(REAL64), intent(in)    ::  d2anlpdrk2(mxddim,mxdaso,3,3)  !<  d^2 anlga / d rkpt^2
       complex(REAL64), intent(in)    ::  d2anlmdrk2(mxddim,mxdaso,3,3)  !<  d^2 anlga / d rkpt^2
       real(REAL64), intent(in)           ::  xnlkb(mxdaso)              !<  Kleinman-Bylander normalization

!      output

       complex(REAL64), intent(out)       ::  vnl(2*mxdbnd,2*mxdbnd)     !<  <Psi|V_NL|Psi>
       complex(REAL64), intent(out)    ::  dvnldrk(2*mxdbnd,2*mxdbnd,3)  !<  d <Psi|V_NL|Psi> d k
       complex(REAL64), intent(out)  :: d2vnldrk2(2*mxdbnd,2*mxdbnd,3,3) !<  d^2 <Psi|V_NL|Psi> d k^2

!      local variables

       complex(REAL64), allocatable       ::  dhdp(:,:)
       complex(REAL64), allocatable       ::  xdhdp(:,:)
       complex(REAL64), allocatable       ::  dhdm(:,:)
       complex(REAL64), allocatable       ::  xdhdm(:,:)

       complex(REAL64), allocatable       ::  ddhdpdrk(:,:,:)
       complex(REAL64), allocatable       ::  ddhdmdrk(:,:,:)
       complex(REAL64), allocatable       ::  xddhdpdrk(:,:,:)
       complex(REAL64), allocatable       ::  xddhdmdrk(:,:,:)

       complex(REAL64), allocatable       ::  d2dhdpdrk2(:,:,:,:)
       complex(REAL64), allocatable       ::  d2dhdmdrk2(:,:,:,:)

       complex(REAL64),allocatable        ::  vnlhalf(:,:)

!      constants

       real(REAL64), parameter :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
       complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
       complex(REAL64), parameter  ::  C_UM = cmplx(UM,ZERO,REAL64)

!      counters

       integer    ::  i, j, m, n

       if(nanlso > 0) then

         allocate(dhdp(nanlso,neig))
         allocate(xdhdp(nanlso,neig))
         allocate(dhdm(nanlso,neig))
         allocate(xdhdm(nanlso,neig))
         
         allocate(ddhdpdrk(nanlso,neig,3))
         allocate(ddhdmdrk(nanlso,neig,3))

         if(nder > 1) then
           allocate(xddhdpdrk(nanlso,neig,3))
           allocate(xddhdmdrk(nanlso,neig,3))
           allocate(d2dhdpdrk2(nanlso,neig,3,3))
           allocate(d2dhdmdrk2(nanlso,neig,3,3))
         else
           allocate(xddhdpdrk(1,1,3))
           allocate(xddhdmdrk(1,1,3))
           allocate(d2dhdpdrk2(1,1,3,3))
           allocate(d2dhdmdrk2(1,1,3,3))
         endif

         allocate(vnlhalf(neig,neig))

!        dhd = < anl | psi >

         call zgemm('c','n',nanlso,neig,mtxd,C_UM,anlp,mxddim,psi,       &
     &                mxddim,C_ZERO,dhdp,nanlso)

         call zgemm('c','n',nanlso,neig,mtxd,C_UM,anlm,mxddim,psi,       &
     &                mxddim,C_ZERO,dhdm,nanlso)

         do j = 1,3

           call zgemm('c','n',nanlso,neig,mtxd,C_UM,danlpdrk(1,1,j),     &
     &                mxddim,psi,mxddim,C_ZERO,ddhdpdrk(1,1,j),nanlso)

           call zgemm('c','n',nanlso,neig,mtxd,C_UM,danlmdrk(1,1,j),     &
     &                mxddim,psi,mxddim,C_ZERO,ddhdmdrk(1,1,j),nanlso)

         enddo

         if(nder > 1) then
           do i = 1,3
           do j = 1,3

             call zgemm('c','n',nanlso,neig,mtxd,C_UM,                   &
     &                                d2anlpdrk2(1,1,j,i),mxddim,        &
     &              psi,mxddim,C_ZERO,d2dhdpdrk2(1,1,j,i),nanlso)

             call zgemm('c','n',nanlso,neig,mtxd,C_UM,                   &
     &                                d2anlmdrk2(1,1,j,i),mxddim,        &
     &              psi,mxddim,C_ZERO,d2dhdmdrk2(1,1,j,i),nanlso)

           enddo
           enddo
         endif

!        xdhd := Diag(xnl) dhd

         do n = 1,neig
           do m = 1,nanlso
             xdhdp(m,n) = xnlkb(m)*dhdp(m,n)
           enddo
         enddo

         do n = 1,neig
           do m = 1,nanlso
             xdhdm(m,n) = xnlkb(m)*dhdm(m,n)
           enddo
         enddo

         if(nder > 1) then

           do j = 1,3
             do n = 1,neig
               do m = 1,nanlso
                 xddhdpdrk(m,n,j) = xnlkb(m)*ddhdpdrk(m,n,j)
               enddo
             enddo
           enddo

           do j = 1,3
             do n = 1,neig
               do m = 1,nanlso
                 xddhdmdrk(m,n,j) = xnlkb(m)*ddhdmdrk(m,n,j)
               enddo
             enddo
           enddo

         endif

!        vnl := < psi | anl > Diag(xnl) < anl | psi > 

         call zgemm('c','n',neig,neig,nanlso,C_UM,dhdp,nanlso,xdhdp,     &
     &                 nanlso,C_ZERO,vnlhalf,mxdbnd)

         do m = 1,neig
         do n = 1,neig
           vnl(n,m) = vnlhalf(n,m)
         enddo
         enddo

         call zgemm('c','n',neig,neig,nanlso,C_UM,dhdm,nanlso,xdhdm,     &
     &                nanlso,C_ZERO,vnlhalf,mxdbnd)

         do m = 1,neig
         do n = 1,neig
           vnl(n+neig,m+neig) = vnlhalf(n,m)
         enddo
         enddo

         call zgemm('c','n',neig,neig,nanlso,C_UM,dhdp,nanlso,xdhdm,     &
     &                nanlso,C_ZERO,vnlhalf,mxdbnd)

         do m = 1,neig
         do n = 1,neig
           vnl(n,m+neig) = vnlhalf(n,m)
         enddo
         enddo

         do m = 1,neig
         do n = 1,neig
           vnl(n+neig,m) = conjg(vnl(m,n+neig))
         enddo
         enddo

         do j = 1,3

           call zgemm('c','n',neig,neig,nanlso,C_UM,ddhdpdrk(1,1,j),     &
     &              nanlso,xdhdp,nanlso,C_ZERO,vnlhalf,mxdbnd)

           call zgemm('c','n',neig,neig,nanlso,C_UM,xdhdp,nanlso,        &
     &                 ddhdpdrk(1,1,j),nanlso,C_UM,vnlhalf,mxdbnd)

           do m = 1,neig
           do n = 1,neig
             dvnldrk(n,m,j) = vnlhalf(n,m)
           enddo
           enddo

           call zgemm('c','n',neig,neig,nanlso,C_UM,ddhdmdrk(1,1,j),     &
     &              nanlso,xdhdm,nanlso,C_ZERO,vnlhalf,mxdbnd)

           call zgemm('c','n',neig,neig,nanlso,C_UM,xdhdm,nanlso,        &
     &                 ddhdmdrk(1,1,j),nanlso,C_UM,vnlhalf,mxdbnd)

           do m = 1,neig
           do n = 1,neig
             dvnldrk(n+neig,m+neig,j) = vnlhalf(n,m)
           enddo
           enddo

           call zgemm('c','n',neig,neig,nanlso,C_UM,ddhdpdrk(1,1,j),     &
     &              nanlso,xdhdm,nanlso,C_ZERO,vnlhalf,mxdbnd)

           call zgemm('c','n',neig,neig,nanlso,C_UM,xdhdp,nanlso,        &
     &                 ddhdmdrk(1,1,j),nanlso,C_UM,vnlhalf,mxdbnd)

           do m = 1,neig
           do n = 1,neig
             dvnldrk(n,m+neig,j) = vnlhalf(n,m)
           enddo
           enddo

           do m = 1,neig
           do n = 1,neig
             dvnldrk(n+neig,m,j) = conjg(dvnldrk(m,n+neig,j))
           enddo
           enddo

         enddo

         if(nder > 1) then

           do i = 1,3
           do j = 1,3

             call zgemm('c','n',neig,neig,nanlso,C_UM,                   &
     &                                        d2dhdpdrk2(1,1,j,i),       &
     &                  nanlso,xdhdp,nanlso,C_ZERO,vnlhalf,mxdbnd)

             call zgemm('c','n',neig,neig,nanlso,C_UM,xdhdp,             &
     &                  nanlso,d2dhdpdrk2(1,1,j,i),nanlso,C_UM,          &
     &                  vnlhalf,mxdbnd)

             call zgemm('c','n',neig,neig,nanlso,2*C_UM,                 &
     &                  xddhdpdrk(1,1,j),                                &
     &                  nanlso,ddhdpdrk(1,1,i),nanlso,C_UM,              &
     &                  vnlhalf,mxdbnd)

             do m = 1,neig
             do n = 1,neig
               d2vnldrk2(n,m,j,i) = vnlhalf(n,m)
             enddo
             enddo

             call zgemm('c','n',neig,neig,nanlso,C_UM,                   &
     &                                        d2dhdmdrk2(1,1,j,i),       &
     &                  nanlso,xdhdm,nanlso,C_ZERO,vnlhalf,mxdbnd)

             call zgemm('c','n',neig,neig,nanlso,C_UM,xdhdm,             &
     &                  nanlso,d2dhdmdrk2(1,1,j,i),nanlso,C_UM,          &
     &                  vnlhalf,mxdbnd)

             call zgemm('c','n',neig,neig,nanlso,2*C_UM,                 &
     &                  xddhdmdrk(1,1,j),                                &
     &                  nanlso,ddhdmdrk(1,1,i),nanlso,C_UM,              &
     &                  vnlhalf,mxdbnd)

             do m = 1,neig
             do n = 1,neig
               d2vnldrk2(n+neig,m+neig,j,i) = vnlhalf(n,m)
             enddo
             enddo

             call zgemm('c','n',neig,neig,nanlso,C_UM,                   &
     &                                        d2dhdpdrk2(1,1,j,i),       &
     &                  nanlso,xdhdm,nanlso,C_ZERO,vnlhalf,mxdbnd)

             call zgemm('c','n',neig,neig,nanlso,C_UM,xdhdp,             &
     &                  nanlso,d2dhdmdrk2(1,1,j,i),nanlso,C_UM,          &
     &                  vnlhalf,mxdbnd)

             call zgemm('c','n',neig,neig,nanlso,2*C_UM,                 &
     &                  xddhdpdrk(1,1,j),                                &
     &                  nanlso,ddhdmdrk(1,1,i),nanlso,C_UM,              &
     &                  vnlhalf,mxdbnd)

             do m = 1,neig
             do n = 1,neig
               d2vnldrk2(n,m+neig,j,i) = vnlhalf(n,m)
             enddo
             enddo

           enddo
           enddo

           do i = 1,3
           do j = 1,3

             do m = 1,neig
             do n = 1,neig
               d2vnldrk2(n+neig,m,j,i) = conjg(d2vnldrk2(m,n+neig,i,j))
             enddo
             enddo

           enddo
           enddo

         endif


         deallocate(dhdp,dhdm)
         deallocate(xdhdp,xdhdm)
         deallocate(ddhdpdrk,ddhdmdrk)

         deallocate(xddhdpdrk,xddhdmdrk)
         deallocate(d2dhdpdrk2,d2dhdmdrk2)

       else

         do n = 1,2*neig
         do m = 1,2*neig
           vnl(m,n) = C_ZERO
         enddo
         enddo

         do j = 1,3
           do n = 1,2*neig
           do m = 1,2*neig
             dvnldrk(m,n,j) = C_ZERO
           enddo
           enddo
         enddo

         if(nder > 1) then
           do i = 1,3
           do j = 1,3
             do n = 1,2*neig
             do m = 1,2*neig
               d2vnldrk2(m,n,j,i) = C_ZERO
             enddo
             enddo
           enddo
           enddo
         endif

       endif

       return

       end subroutine psi_vnl_sp_psi_der
