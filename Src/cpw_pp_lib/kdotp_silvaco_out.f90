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

!>     Writes to a file the data to interface with silvaco software
!>     Silvaco never finished it... but can be used later.

       subroutine kdotp_silvaco_out(filename,iotape,neig,adot,rk0,       &
     &     h0,dh0drk,d2h0drk2,natot,neltot,emidgap,lso,                  &
     &     mxdbnd)

!      Written April 12, 2014. jlm
!      modified documentation 19 February 2020. JLM
!      copyright  Jose Luis Martins/INESC-MN

!      version 4.95

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdbnd                     !<  array dimension for the number of bands

       character(len=*), intent(in)       ::  filename                   !<  file to be written
       integer, intent(in)                ::  iotape                     !<  tape number 

       integer, intent(in)                ::  neig                       !<  number of bands
       real(REAL64), intent(in)           ::  adot(3,3)                  !<  metric in direct space
       real(REAL64), intent(in)           ::  rk0(3)                     !<  reference k-point (in lattice coordinates)
       complex(REAL64), intent(in)        ::  h0(mxdbnd,mxdbnd)          !<  <Psi|H|Psi> 
       complex(REAL64), intent(in)        :: dh0drk(mxdbnd,mxdbnd,3)     !<  d <Psi|H|Psi> d k
       complex(REAL64), intent(in)        :: d2h0drk2(mxdbnd,mxdbnd,3,3) !<  d^2 <Psi|H|Psi> d k^2

       integer, intent(in)                ::  natot                      !<  total number of atoms
       integer, intent(in)                ::  neltot                     !<  total number of electrons
       real(REAL64), intent(in)           ::  emidgap                    !<  rough estimate of the mid-gap
       logical, intent(in)                ::  lso                        !<  true if spin-orbit is included

!      local variables

       real(REAL64)      ::  avec(3,3)        !  a(j,i) is the j-th component of primitive lattice vector i 
       real(REAL64)      ::  bvec(3,3)        !  b(j,i) is the j-th component of reciprocal lattice vector i 


!      counters

       integer               ::  i, j, n, m

!      constants

       real(REAL64), parameter :: HARTREE = 27.21138386_REAL64
       real(REAL64), parameter :: BOHR = 0.05291772109



       open(unit=iotape,file=filename,form='formatted')

       write(iotape,'(2x,i6,"   # Hamiltonian dimension ")') neig

       write(iotape,'(2x,i6,"   # number of atoms in unit cell")') natot

       write(iotape,'(2x,i6,3x,f14.6,"   # number of electrons in unit", &
     &   " cell, midgap level ESTIMATE (eV)")') neltot,emidgap*HARTREE

       if(lso) then
         write(iotape,'("  1   # band degeneracy (H with spin-orbit)")') 
       else
         write(iotape,'("  2   # band degeneracy (no spin)")') 
       endif

       call adot_to_avec_sym(adot,avec,bvec)
       
       do j=1,3
         write(iotape,'(2x,3(1x,f14.6),"  # Reciprocal vector",i3,       &
     &     "  in nm^-1 ")') (bvec(i,j)/BOHR,i=1,3),j
       enddo

       write(iotape,'(2x,3f14.8,"  # Reference k-point in reciprocal ",  &
     &     "lattice coordinates")') (rk0(j),j=1,3)
     
       write(iotape,'("#",5x,"index",11x,"H_0",36x,"d H_0 / d k",94x,     &
     &         "d2 H_0 / d k2")')
     
       do i=1,neig
       do j=1,neig
         write(iotape,'(2(1x,i5),6x,2(1x,f14.8),6x,3(3x,2(1x,f14.8)),    &
     &         6x,9(3x,2(1x,f14.8)))') i,j,                              &
     &              h0(i,j)*HARTREE,(dh0drk(i,j,n)*HARTREE,n=1,3),       &
     &             ((d2h0drk2(i,j,n,m)*HARTREE,n=1,3),m=1,3)
       enddo
       enddo

       close(unit=iotape)

       return

       end subroutine kdotp_silvaco_out
