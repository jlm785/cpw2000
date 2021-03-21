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

!>     Reads the file with the data to interface with kdotp calculations

       subroutine kdotp_in(filename,iotape,neig,adot,rk0,                &
     &     h0,dh0drk,d2h0drk2,                                           &
     &     mxdbnd)

!      Written May 3, 2014, from the writing subroutine. jlm
!      Modified, documentation, name, 15 June 2020. JLM
!      copyright  Jose Luis Martins/INESC-MN

!      version 4.98

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdbnd                     !<  array dimension for the number of bands

       character(len=*), intent(in)       ::  filename                   !<  file to be written
       integer, intent(in)                ::  iotape                     !<  tape number 

!      output

       integer, intent(out)               ::  neig                       !<  number of bands
       real(REAL64), intent(out)          ::  adot(3,3)                  !<  metric in direct space
       real(REAL64), intent(out)          ::  rk0(3)                     !<  reference k-point (in lattice coordinates)
       complex(REAL64), intent(out)       ::  h0(mxdbnd,mxdbnd)          !<  <Psi|H|Psi> 
       complex(REAL64), intent(out)       :: dh0drk(mxdbnd,mxdbnd,3)     !<  d <Psi|H|Psi> d k
       complex(REAL64), intent(out)       :: d2h0drk2(mxdbnd,mxdbnd,3,3) !<  d^2 <Psi|H|Psi> d k^2

!      local variables

       real(REAL64)      ::  bvec(3,3)        !  b(j,i) is the j-th component of reciprocal lattice vector i 
       real(REAL64)      ::  vcell, bdot(3,3)

!      unused variables

       integer                ::  natot                      !  total number of atoms
       integer                ::  neltot                     !  total number of electrons
       real(REAL64)           ::  emidgap                    !  rough estimate of the mid-gap
       integer                ::  ldeg                       !  band degeneracy

!      counters

       integer               ::  i, j, n, m

!      constants

       real(REAL64), parameter :: ZERO = 0.0_REAL64
       real(REAL64), parameter :: HARTREE = 27.21138386_REAL64
       real(REAL64), parameter :: BOHR = 0.05291772109_REAL64


       open(unit=iotape,file=filename,form='formatted',status='old')

       read(iotape,*) neig
       read(iotape,*) natot
       read(iotape,*) neltot,emidgap
       emidgap = emidgap / HARTREE
       read(iotape,*) ldeg
       
       do j=1,3
         read(iotape,'(2x,3(1x,f14.6))') (bvec(i,j),i=1,3)
         do i=1,3
           bvec(i,j) = bvec(i,j)*BOHR
         enddo
       enddo

       do j=1,3
       do i=1,3
         bdot(i,j) = ZERO
         do n=1,3
           bdot(i,j) = bdot(i,j) + bvec(i,n)*bvec(j,n)
         enddo
       enddo
       enddo
       
       call adot_to_bdot(bdot,vcell,adot)

       read(iotape,'(2x,3f14.8)') (rk0(j),j=1,3)

       read(iotape,*)

       do i=1,neig
       do j=1,neig
         read(iotape,'(18x,2(1x,f14.8),6x,3(3x,2(1x,f14.8)),6x,          &
     &        9(3x,2(1x,f14.8)))')  h0(i,j),(dh0drk(i,j,n),n=1,3),       &
     &             ((d2h0drk2(i,j,n,m),n=1,3),m=1,3)
         h0(i,j) = h0(i,j) / HARTREE
         do n=1,3
           dh0drk(i,j,n) = dh0drk(i,j,n) / HARTREE
           do m=1,3
             d2h0drk2(i,j,n,m) = d2h0drk2(i,j,n,m) / HARTREE
           enddo
         enddo
       enddo
       enddo

       close(unit=iotape)

       return

       end subroutine kdotp_in
