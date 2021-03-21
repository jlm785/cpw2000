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

!>     computes, symmetrizes, and adds the charge density
!>     from the eigenvectors.
!>     Should be compiled with openmp

       subroutine charge_by_fft(mtxd,neig,occp,isort,psi,denk,           &
     & ng,kgv,phase,conj,ns,inds,kmax,mstar,                             &
     & mxddim,mxdbnd,mxdgve,mxdnst)

!      written august 14 1987.jlm
!      modified december 9 1993. jlm
!      modified 22 march 1999. jlm
!      modified (chd) 26 july 2002. jlm
!      modified for complex version and f90, 21 September 2015. JLM
!      modified den_sym, mesh_fold, 30 September 2015. JLM
!      Modified 13, December 2019.  denu allocation.  JLM
!      Modified, documentation, January 2020. JLM
!      copyright INESC-MN/Jose Luis Martins

!      version 4.94


       implicit none
       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxddim                     !<  array dimension of plane-waves
       integer, intent(in)                ::  mxdbnd                     !<  array dimension for number of bands
       integer, intent(in)                ::  mxdgve                     !<  array dimension for g-space vectors
       integer, intent(in)                ::  mxdnst                     !<  array dimension for g-space stars

       integer, intent(in)                ::  mtxd                       !<  dimension of the hamiltonian
       integer, intent(in)                ::  neig                       !<  number of eigenvectors
       real(REAL64), intent(in)           ::  occp(mxdbnd)               !<  fractional ocupation of level j
       integer, intent(in)                ::  isort(mxddim)              !<  g-vector associated with row/column i of hamiltonian

       complex(REAL64), intent(in)        ::  psi(mxddim,mxdbnd)         !<  |psi> 

       integer, intent(in)                ::  ng                         !<  size of g-space
       integer, intent(in)                ::  kgv(3,mxdgve)              !<  G-vectors in reciprocal lattice coordinates 
       complex(REAL64), intent(in)        ::  phase(mxdgve)              !<  phase factor of G-vector n
       real(REAL64), intent(in)           ::  conj(mxdgve)               !<  is -1 if one must take the complex conjugate of x*phase
       integer, intent(in)                ::  ns                         !<  number os stars with length less than gmax
       integer, intent(in)                ::  inds(mxdgve)               !<  star to which g-vector n belongs
       integer, intent(in)                ::  kmax(3)                    !<  max value of |kgv(i,n)|
       integer, intent(in)                ::  mstar(mxdnst)              !<  number of g-vectors in the j-th star

!      output

       complex(REAL64), intent(out)       ::  denk(mxdnst)               !<  symmetrized density (in stars of G)

!      local allocatable arrays

       complex(REAL64), allocatable       ::  denu(:)                     !  unsymmetrized density (in G)
       real(REAL64), allocatable          ::  rhomsh(:)
       real(REAL64), allocatable          ::  wrkfft(:)
       complex(REAL64), allocatable       ::  chd(:)
       integer, allocatable               ::  ipoint(:)

!      local variables

       integer         ::  mxdfft,mxdwrk
       integer         ::  jmin, jmax
       integer         ::  nsfft(3)
       integer         ::  n1, n2, n3, id
       integer         ::  nn1, nn2, nn3
       integer         ::  kd1, kd2, kd3
       integer         ::  k1, k2, k3
       integer         ::  ntot, it

!      counters

       integer         ::  i, j
 
!      parameters

       real(REAL64), parameter :: SMALL = 0.000000000001_REAL64
       real(REAL64), parameter :: ZERO = 0.0_REAL64
       complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)


       if (neig > 0) then

!        find min and max band with nonzero occupation

         jmin = 0
         jmax = 0
         do i=1,neig
           if (abs(occp(i)) > SMALL .and. jmin == 0) jmin = i
           if (abs(occp(i)) > SMALL .and. jmin /= 0) jmax = i
         enddo

         if (jmin > 0) then

!          find n for fast fourier transform
!          ni is the number of points used in direction i.

           call size_fft(kmax,nsfft,mxdfft,mxdwrk)
!          call sizfft(kmax,id,n1,n2,n3,mxdfft,mxdwrk)

           allocate(chd(mxdfft))
           allocate(wrkfft(mxdwrk))

           n1 = nsfft(1)
           n2 = nsfft(2)
           n3 = nsfft(3)
           id = n1
           ntot = id * n2 * n3
           nn1 = (n1-1) / 2
           nn2 = (n2-1) / 2
           nn3 = (n3-1) / 2
           kd1 = 0
           kd2 = 0
           kd3 = 0

!          write(6,*) '  n1,n2,n3 = ',n1,n2,n3

!          fills array ipoint (gather-scatter index)

           allocate(ipoint(mtxd))

           do i = 1,mtxd
             it = isort(i)
             k1 = kgv(1,it)
             if (iabs(k1) > kd1) kd1 = iabs(k1)
             if (k1 < 0) k1 = n1 + k1
             k2 = kgv(2,it)
             if (iabs(k2) > kd2) kd2 = iabs(k2)
             if (k2 < 0) k2 = n2 + k2
             k3 = kgv(3,it)
             if (iabs(k3) > kd3) kd3 = iabs(k3)
             if (k3 < 0) k3 = n3 + k3
             ipoint(i) = (k3*n2 + k2)*id + k1 + 1
           enddo
           if (kd1 > nn1 .or. kd2 > nn2 .or. kd3 > nn3) then
             write(6,*)
             write(6,'("     STOPPED in charge_by_fft      ",            &
     &        "size of matrix ",i7," fft mesh ",3i5)') mtxd,n1,n2,n3

             stop

           endif
       
!          initialize rhomsh

           allocate(rhomsh(mxdfft))

           do i = 1,ntot
             rhomsh(i) = ZERO
           enddo

!          start loop over eigenvectors

           do j = jmin,jmax

!$omp parallel do default(shared) private(i)
             do i=1,ntot
               chd(i) = C_ZERO
             enddo
!$omp end parallel do

!$omp parallel do default(shared) private(i)
             do i=1,mtxd
               chd(ipoint(i)) = psi(i,j)
             enddo
!$omp end parallel do

!            fourier transform to real space

             call cfft_wf_c16(chd,id,n1,n2,n3,kd1,kd2,kd3,-1,            &
     &          wrkfft,mxdwrk)

!            calculates the square of chd (wave-function -> charge)

!$omp parallel do default(shared) private(i)
             do i=1,ntot
               rhomsh(i) = rhomsh(i) + occp(j)*                          &
     &               real(chd(i)*conjg(chd(i)),REAL64)
             enddo
!$omp end parallel do

           enddo                  !  end loop over eigenvectors

!          fourier transform to momentum space

           do i=1,ntot
             chd(i) = cmplx(rhomsh(i),ZERO,REAL64)
           enddo
           
           deallocate(rhomsh)
           deallocate(ipoint)

           call rfft_c16(chd,id,n1,n2,n3,1,wrkfft,mxdwrk)

!          initialize denu

           allocate(denu(mxdgve))

           call mesh_fold(denu,chd,id,n1,n2,n3,                          &
     &     ng,kgv,                                                       &
     &     mxdgve,mxdfft)

           deallocate(chd)
           deallocate(wrkfft)


!          CONJ SHOULD BE CONVERTED TO INTEGER OR LOGICAL

           call star_of_g_fold(denk,denu,.FALSE.,                        &
     &     ng,phase,conj,ns,inds,mstar,                                  &
     &     mxdgve,mxdnst)

           deallocate(denu)

         endif
       
       endif

       return
       end subroutine charge_by_fft
