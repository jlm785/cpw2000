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

!>     Tries to match the bands at a given k-point with the prvious
!>     k-point using functions overlap
!>     To be precise it should use symmetry...

       subroutine out_band_match(irk,nlines,ljump,nkstep,                   &
     &    neig,icmplx,xsc,                                                  &
     &    mtxd,psi,isort,psi_so,ei,ei_so,                                   &
     &    e_of_k,e_of_k_so,nrk2,                                            &
     &    mxddim,mxdbnd)

!      extracted from earlier version 4.95 of out_band.f90. JLM
!      Removed icmplx from calls, 10 June 2020. JLM
!      copyright  Jose Luis Martins/INESC-MN

!      version 4.98

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)


!      input

       integer, intent(in)                ::  mxddim                     !<  array dimension of plane-waves
       integer, intent(in)                ::  mxdbnd                     !<  array dimension for number of bands

       integer, intent(in)                ::  nlines                     !<  number of lines in reciprocal space
       integer, intent(in)                ::  nrk2                       !<  total number of k-points

       integer, intent(in)                ::  irk                        !<  index of k-point
       logical, intent(in)                ::  ljump(nlines)              !<  indicates if the new line contains a jump from the preceeding
       integer, intent(in)                ::  nkstep(nlines)             !<  number of steps in line

       integer, intent(in)                ::  neig                       !<  number of eigenvectors
       integer, intent(in)                ::  icmplx                     !<  indicates if the structure factor is complex (not used in present implementation)

       real(REAL64), intent(in)           ::  xsc                        !<  criteria for reduced vector size

       real(REAL64), intent(in)           ::  ei(mxdbnd)                 !<  eigenvalue no. i. (hartree)
       real(REAL64), intent(in)           ::  ei_so(2*mxdbnd)            !<  eigenvalue no. i. (hartree)

       integer, intent(in)                ::  mtxd                       !<  number of coefficients of psi
       complex(REAL64), intent(in)        ::  psi(mxddim, mxdbnd)        !<  component j of eigenvector i
       integer, intent(in)                ::  isort(mxddim)              !<  g-vector associated with row/column i of hamiltonian
       complex(REAL64), intent(in)        ::  psi_so(2*mxddim,2*mxdbnd)  !<  component j of eigenvector i

!      input and output

       real(REAL64), intent(inout)        ::  e_of_k(neig,nrk2)          !<  band energies of k-point in plot
       real(REAL64), intent(inout)        ::  e_of_k_so(2*neig,nrk2)     !<  spin-orbit band energies of k-point in plot

!      local allocatable arrays

       integer, allocatable, save         ::  iperm(:)                   !  permutation of the eigenvalues
       integer, allocatable, save         ::  ipermold(:)                !  permutation of the eigenvalues

       integer, allocatable, save         ::  iperm_so(:)                !  permutation of the eigenvalues
       integer, allocatable, save         ::  ipermold_so(:)             !  permutation of the eigenvalues

       integer, allocatable, save         ::  imatch(:)                  !  points to the old state that is most similar to the new state
       integer, allocatable, save         ::  imatch_so(:)               !  points to the old state that is most similar to the new state

       integer, allocatable, save         ::  isold(:)                   !  g-vector associated with row/column i of old hamiltonian

       complex(REAL64), allocatable, save ::  psiold(:,:)                !  old eigenvectors
       complex(REAL64), allocatable, save ::  psiold_so(:,:)             !  old eigenvectors spin-orbit

!      local variables

       integer, save                      ::  nextline                   !  indicates next vertical line on band circuit
       integer, save                      ::  nl                         !  keeps track of next lines

       integer, save                      ::  mxdold                     !  psiold mxddim
       integer, save                      ::  mtxdold                    !  number of coefficients in old vectors

       real(REAL64)                       ::  xsum
 
!      constants
 
       real(REAL64), parameter  :: ZERO = 0.0_REAL64

!      counters

       integer    ::  j, n


       if(irk == 1) then

         nextline = 1
         nl = 1

         xsum = ZERO
         do j = 1,mtxd
           xsum = xsum + real(psi(j,neig)*conjg(psi(j,neig)),REAL64)
           mxdold = j
           if(xsum > xsc) exit
         enddo
         mtxdold = mxdold

         allocate(iperm(mxdbnd),ipermold(mxdbnd))
         allocate(iperm_so(2*mxdbnd),ipermold_so(2*mxdbnd))
         allocate(imatch(mxdbnd))
         allocate(imatch_so(2*mxdbnd))
         allocate(isold(mxddim))

         allocate(psiold(mxdold,mxdbnd))
         allocate(psiold_so(2*mxdold,2*mxdbnd))

       endif

       if(irk /= nextline .or. .not. ljump(nl)) then

         call out_band_match_state(imatch,neig,                        &
     &      mtxdold,psiold,isold,mtxd,psi,isort,                       &
     &      mxddim,mxdbnd,mxdold)

         call out_band_match_state_so(imatch_so,neig,                  &
     &      mtxdold,psiold_so,isold,mtxd,psi_so,isort,                 &
     &      mxddim,mxdbnd,mxdold)

       endif

!      Defines mtxdold such that it
!      has most of the spectral weight for highest eigenvalue

       xsum = ZERO
       do j = 1,min(mtxd,mxdold)
         xsum = xsum + real(psi(j,neig)*conjg(psi(j,neig)),REAL64)
         mtxdold = j
         if(xsum > xsc) exit
       enddo


       do n=1,neig
       do j=1,mtxdold
         psiold(j,n) = psi(j,n)
       enddo
       enddo
       
       do n=1,2*neig
       do j=1,2*mtxdold
         psiold_so(j,n) = psi_so(j,n)
       enddo
       enddo

       do j=1,mtxd
         isold(j) = isort(j)
       enddo

       if(irk == nextline .and. ljump(nl)) then

         do j=1,neig
           iperm(j) = j
         enddo
         do j=1,2*neig
           iperm_so(j) = j
         enddo

       else

         do j=1,neig
           ipermold(j) = iperm(j)
         enddo
         do j=1,neig
           iperm(j) = imatch(ipermold(j))
         enddo
         do j=1,2*neig
           ipermold_so(j) = iperm_so(j)
         enddo
         do j=1,2*neig
           iperm_so(j) = imatch_so(ipermold_so(j))
         enddo

       endif

       if(irk == nextline) then
         if(ljump(nl)) then
           nextline = nextline + nkstep(nl) + 1
         else
           nextline = nextline + nkstep(nl)
         endif
         if(nl < nlines) nl = nl + 1
       endif

       do j=1,neig
         e_of_k(j,irk) = ei(iperm(j))
       enddo
       
       do j=1,2*neig
         e_of_k_so(j,irk) = ei_so(iperm_so(j))
       enddo

       if(irk == nrk2) then

         deallocate(iperm,ipermold)
         deallocate(iperm_so,ipermold_so)
         deallocate(imatch)
         deallocate(imatch_so)
         deallocate(isold)

         deallocate(psiold)
         deallocate(psiold_so)

       endif

       return
       end subroutine out_band_match
