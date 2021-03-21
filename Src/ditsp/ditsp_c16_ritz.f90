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

!>     Ritz step of the iterative diagonalization

       subroutine ditsp_c16_ritz(lsafe,ipr,ihkp,epsa,iseed,lnewanl,       &
     & mtxd,neig,nexact,nitold,niter,nconv,                               &
     & psi,hpsi,bas,hbas,                                                 &
     & eg,lconv,                                                          &
     & ekpg,isort,vscr,kmscr,                                             &
     & ng,kgv,                                                            &
     & anlga,xnlkb,nanl,                                                  &
     & mxddim,mxdsml,mxdbnd,mxdgve,mxdscr,mxdanl)

!      Adapted from the 1990 code.
!      Written 12 August 2015. JLM
!      Modified kmscr, 28 October 2015. JLM
!      Modified, documentation, January 2020. JLM
!      Modified, qmod-->ekpg in hk_psi. 13 February 2021. JLM
!      copyright  Jose Luis Martins/INESC-MN

!      version 4.99

       implicit none
       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxddim                     !<  array dimension of plane-waves
       integer, intent(in)                ::  mxdsml                     !<  array dimension of small hamiltonian
       integer, intent(in)                ::  mxdbnd                     !<  array dimension for number of bands
       integer, intent(in)                ::  mxdgve                     !<  array dimension for g-space vectors
       integer, intent(in)                ::  mxdscr                     !<  array dimension for fft transform
       integer, intent(in)                ::  mxdanl                     !<  array dimension for number of KB projectors

       integer, intent(in)                ::  ipr                        !<  print control. 0 no printing, 2 lots of stuff
       real(REAL64), intent(in)           ::  epsa                       !<  eigenvalue energy scale

       integer, intent(in)                ::  mtxd                       !<  dimension of the hamiltonian
       integer, intent(in)                ::  nexact                     !<  number of basis functions in bas that are never changed
       integer, intent(in)                ::  nitold                     !<  old number of vectors in iterative space

       real(REAL64), intent(in)           ::  ekpg(mxddim)                    !<  kinetic energy (hartree) of k+g-vector of row/column i
       integer, intent(in)                ::  isort(mxddim)              !<  g-vector associated with row/column i of hamiltonian
       real(REAL64), intent(in)           ::  vscr(mxdscr)               !<  screened potential in the fft real space mesh
       integer, intent(in)                ::  kmscr(7)                   !<  max value of kgv(i,n) used for the potential fft mesh and fft mesh size

       integer, intent(in)                ::  ng                         !<  total number of g-vectors with length less than gmax
       integer, intent(in)                ::  kgv(3,ng)                  !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

       complex(REAL64), intent(in)        ::  anlga(mxddim,mxdanl)       !<  KB projectors
       real(REAL64), intent(in)           ::  xnlkb(mxdanl)              !<  KB projector normalization
       integer, intent(in)                ::  nanl                       !<  half of number of projectors without spin

!      output

       real(REAL64),intent(out)           ::  eg(mxddim)                 !<  eigenvalues CORRECT DIMENSIONS LATER
       logical,intent(out)                ::  lconv                      !<  indicates convergence (should never occur)
       logical, intent(in)                ::  lsafe                      !<  indicates that no problems were identified in the iterative basis
       integer, intent(out)               ::  ihkp                       !<  number of H |psi> products

!      input and output

       integer, intent(inout)             ::  iseed                      !<  random number generator seed
       integer, intent(inout)             ::  neig                       !<  number of eigenvectors
       integer, intent(inout)             ::  niter                      !<  number of vectors in iterative space
       integer, intent(inout)             ::  nconv                      !<  number of converged vectors

       logical, intent(inout)             ::  lnewanl                    !<  indicates that anlga has been recalculated (not used in default implementation)

       complex(REAL64), intent(inout)     ::  psi(mxddim,mxdbnd)         !<  |psi>
       complex(REAL64), intent(inout)     ::  hpsi(mxddim,mxdbnd)        !<  H |psi>
       complex(REAL64), intent(inout)     ::  bas(mxddim,mxdsml)         !<  |bas>
       complex(REAL64), intent(inout)     ::  hbas(mxddim,mxdsml)        !<  H |bas>

!      work arrays

       complex(REAL64),allocatable        ::  vsm(:,:)
       complex(REAL64),allocatable        ::  hamsm(:,:)
       integer,allocatable                ::  irow(:)

!      local variables

       integer       :: ired, ndeg, nmax, nsize, nlow, kmax
       real(REAL64)  :: diff, xn, xr, xi
       logical       :: ljump

       integer             ::  info

!      counters

       integer       :: i, k, n

!      parameters

       real(REAL64), parameter :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
       complex(REAL64), parameter  ::  C_UM = cmplx(UM,ZERO,REAL64)
       complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)

!      external functions

       real(REAL64),external      :: ran2
       complex(REAL64),external   :: zdotc


       lconv = .FALSE.

       nsize = nexact + nitold + niter - nconv
       nlow = nexact + nitold - nconv
       if(nsize > mxdsml) then
         write(6,'("   stopped in ditspc_ritz,  change mxdsml to",       &
     &        " at least  ",i5)') nsize

         stop

       endif

       if(nconv < niter) then
         do n=nconv+1,niter
           call zcopy(mtxd,psi(1,n),1,bas(1,nlow+n),1)
         enddo
       endif

!      adds random guess vector in special occasions

       if(.not. lsafe .or. niter == nconv) then
         if(ipr > 2) write(6,'("  random vector added to basis")')
         nsize = nsize + 1
         if(nsize > mxdsml) then
           write(6,'("   stopped in ditspc_ritz,  change mxdsml to",     &
     &          " at least  ",i5)') nsize

           stop

         endif
         do k=1,mtxd
           xr = ran2(iseed)
           xi = ran2(iseed)
           bas(k,nsize) = cmplx(xr,xi,REAL64)
         enddo
       endif

!      Gram-Schmidt orthogonalization

       allocate(irow(mxdsml))

       call grsch_loop_c16(.FALSE.,bas,hbas,mtxd,nexact+nconv,nsize,irow,    &
     & mxddim)

       call grsch_shift_c16(bas,mtxd,nexact+nconv,nsize,irow,ired,       &
     & mxddim)

       deallocate(irow)

       if(ipr > 2 .and. ired > 0) write(6,'("  orthogonalization",       &
     &     " found",i6," linearly dependent eigenvectors")') ired

       nsize = nsize - ired
       if(nsize > mtxd) then
         write(6,'("   stopped in ditspc_ritz,  problems with",          &
     &        " Gram-Schmidt.  nsize =  ",i5)') nsize

         stop

       endif

!      this step is very important to avoid propagation of
!      numerical instabilities even if it appears that
!      computing time is being wasted.

       call hk_psi_c16(mtxd,nsize-nexact-nconv,                          &
     & bas(1,nexact+nconv+1),hbas(1,nexact+nconv+1),lnewanl,             &
     & ng,kgv,                                                           &
     & ekpg,isort,vscr,kmscr,                                            &
     & anlga,xnlkb,nanl,                                                 &
     & mxddim,mxdsml-(nexact+nconv),mxdanl,mxdgve,mxdscr)

       ihkp = nsize - nexact - nconv

!      small matrix construction and diagonalization.
!      only triangular part is calculated. not!!!!

       allocate(hamsm(mxdsml,mxdsml))

       call zgemm('c','n',nsize,nsize,mtxd,C_UM,hbas,mxddim,             &
     &            bas,mxddim,C_ZERO,hamsm,mxdsml)

       allocate(vsm(mxdsml,mxdsml))

       call diag_c16(nsize,hamsm,eg,vsm,mxdsml,info)


       if(info /= 0) stop


       deallocate(hamsm)

!      recalculates ndeg taking into account degeneracies

       ndeg = neig
       if(nsize > neig) then
         do i=neig+1,nsize
           diff = abs(eg(i) - eg(i-1))

           if(diff > epsa) exit

           ndeg = i
         enddo
       endif

       if(ndeg > mxdbnd) ndeg = mxdbnd

!      recalculates niter

       niter = ndeg + 1
       if(nsize > niter) then
         do i=ndeg+2,nsize
           diff = abs(eg(i) - eg(i-1))

           if(diff > epsa) exit

           niter = i
         enddo
       endif


       nmax = (mxdsml - nexact + nconv) / 2
       if(niter > nmax) niter = nmax
       if(niter > nexact) niter = nexact
       if(niter > mxdbnd) niter = mxdbnd
       if(niter < ndeg) niter = ndeg

!      tests that converged eigenvectors did not change
!      and that there are no unexpected degeneracies

       kmax = nconv
       if(nconv > 0) then
         do k=1,kmax
           nmax = nconv
           nconv = 0
           ljump = .TRUE.
           do n=1,nmax
             xn = real(zdotc(nexact+nmax,vsm(1,n),1,              &
     &                                   vsm(1,n),1),REAL64)
             if(xn < 0.99) then
               ljump = .FALSE.

               exit

             else
               nconv = n
             endif
           enddo

           if(ljump) then
             diff = abs(eg(nconv+1)-eg(nconv))
             if(diff < epsa) then
               nconv = nconv-1
               ljump = .FALSE.
             endif
           endif

           if(ljump .or. nconv == 0) exit

         enddo
       endif

       if(ipr > 2) write(6,'("  new values niter=",i6,"  nconv=",        &
     &      i6,"  ndeg=",i6,"  eg=",/,20(2x,5e13.5,/))') niter,nconv,    &
     &      ndeg,(eg(n),n=1,ndeg)

!      constructs eigenvector

       call zgemm('n','n',mtxd,niter,nsize,C_UM,bas,mxddim,              &
     &            vsm,mxdsml,C_ZERO,psi,mxddim)
       call zgemm('n','n',mtxd,niter,nsize,C_UM,hbas,mxddim,             &
     &            vsm,mxdsml,C_ZERO,hpsi,mxddim)

       deallocate(vsm)

!      the copy to hbas was not on original code...

       do n=nconv+1,niter
         call zcopy(mtxd,psi(1,n),1,bas(1,nexact+n),1)
         call zcopy(mtxd,hpsi(1,n),1,hbas(1,nexact+n),1)
       enddo

       if(nconv >= ndeg) then

!        this should not occur here... just being paranoid

         neig = ndeg

         lconv = .TRUE.

       endif

       return

       end subroutine ditsp_c16_ritz
