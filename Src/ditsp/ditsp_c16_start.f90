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

!>     Start-up of the iterative diagonalization.

       subroutine ditsp_c16_start(ipr,iguess,epspsi,epsa,                &
     & lortho,ihkp,lnewanl,                                              &
     & neig,mtxd,mtxds,nexact,niter,nitold,ndeg,                         &
     & psi,hpsi,bas,hbas,                                                &
     & ekpg,isort,vscr,kmscr,                                            &
     & ng,kgv,                                                           &
     & anlga,xnlkb,nanl,                                                 &
     & hamsm,                                                            &
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
       integer, intent(in)                ::  iguess                     !<  if =1 uses neig guess eigenvectors in psir
       real(REAL64), intent(in)           ::  epspsi                     !<  criteria for convergence

       integer, intent(in)                ::  neig                       !<  number of eigenvectors
       integer, intent(in)                ::  mtxd                       !<  dimension of the hamiltonian
       integer, intent(in)                ::  mtxds                      !<  dimension of the small hamiltonian

       real(REAL64), intent(in)           ::  ekpg(mxddim)                    !<  kinetic energy (hartree) of k+g-vector of row/column i
       integer, intent(in)                ::  isort(mxddim)              !<  g-vector associated with row/column i of hamiltonian
       real(REAL64), intent(in)           ::  vscr(mxdscr)               !<  screened potential in the fft real space mesh
       integer, intent(in)                ::  kmscr(7)                   !<  max value of kgv(i,n) used for the potential fft mesh and fft mesh size

       integer, intent(in)                ::  ng                         !<  total number of g-vectors with length less than gmax
       integer, intent(in)                ::  kgv(3,ng)                  !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

       complex(REAL64), intent(in)        ::  anlga(mxddim,mxdanl)       !<  KB projectors
       real(REAL64), intent(in)           ::  xnlkb(mxdanl)              !<  KB projector normalization
       integer, intent(in)                ::  nanl                       !<  half of number of projectors without spin

       complex(REAL64), intent(in)        ::  hamsm(mxdsml,mxdsml)       !<  small hamiltonian

!      output

       real(REAL64), intent(out)          ::  epsa                       !<  eigenvalue energy scale
       integer, intent(out)               ::  nexact                     !<  number of basis functions in bas that are never changed
       integer, intent(out)               ::  nitold, niter              !<  number of vectors in iterative space
       integer, intent(out)               ::  ndeg                       !<  number of eigenvectors including degeneracy
       logical, intent(out)               ::  lortho                     !<  indicates that the |psi> are orthogonal
       integer, intent(out)               ::  ihkp                       !<  numbe rof H |psi> products

       complex(REAL64), intent(out)       ::  hpsi(mxddim,mxdbnd)        !<  H |psi>
       complex(REAL64), intent(out)       ::  bas(mxddim,mxdsml)         !<  |bas>
       complex(REAL64), intent(out)       ::  hbas(mxddim,mxdsml)        !<  H |bas>

!      input and output

       complex(REAL64), intent(inout)     ::  psi(mxddim,mxdbnd)         !<  |psi>

       logical, intent(inout)             ::  lnewanl                    !<  indicates that anlga has been recalculated (not used in default implementation)

!      local variables

       complex(REAL64), allocatable        ::  tvec(:,:)
       real(REAL64), allocatable           ::  et(:)

       integer       :: imax
       real(REAL64)  :: diff, ecut

       integer             ::  info

!      counters

       integer       :: i, j, n

!      parameters

       real(REAL64), parameter :: ZERO = 0.0_REAL64
       complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)

!      initial allocations

       allocate(tvec(mxdsml,mxdsml))
       allocate(et(mxdsml))

!      diagonalizes small matrix

       call diag_c16(mtxds,hamsm,et,tvec,mxdsml,info)


       if(info /= 0) stop


       if(ipr > 2) write(6,'("  results of small matrix ",               &
     &    "diagonalization  mtxds = ",i7,/,20(2x,5e13.5,/))')            &
     &          mtxds,(et(n),n=1,mtxds)

!      criteria for eigenvalue quasi-degeneracy

       epsa = max(epspsi,abs((et(neig)-et(1))/real(50*neig)))

!      finds initial value of ndeg. ndeg takes into account
!      the current pattern of degeneracies.

       ndeg = neig
       if(neig < mtxds) then
         do i = neig+1,mtxds
           diff = abs(et(i) - et(i-1))

           if(diff > epsa) exit

           ndeg = i
         enddo
       endif
       if(ndeg > mxdbnd) ndeg = mxdbnd

!      finds the initial value of niter. niter is the number of
!      vectors used in the iteration.

       niter = ndeg + 1
       if(niter < mtxds) then
         do i = ndeg+2,mtxds
           diff = abs(et(i) - et(i-1))

           if(diff > epsa) exit

           niter = i
         enddo
       else
         niter = mtxds
       endif
       if(niter > mxdbnd) niter = mxdbnd

!      finds nexact. nexact is the number of eigenvectors
!      of the small matrix that are included in the expansion basis
!      set. nexact is never changed, and should be greater than niter.

       nexact = neig
       ecut = 1.2*et(neig) - 0.2*et(1)
       if(neig < mtxds) then
         do i = neig+1,mtxds

           if(et(i) > ecut) exit

           nexact = i
         enddo
       endif
       nexact = min(nexact,mxdsml/3,neig+2+neig/5,mtxds/2) + 3

       imax = nexact-1
       do i = 1,imax
         nexact = nexact - 1
         diff = abs(et(imax-i+2) - et(imax-i+1))

         if(diff > epsa) exit

       enddo
       nexact = max(niter,nexact)

!      this does not make sense... should be checked elsewhere.

       imax = nexact + niter + niter/5 + 2
       if(imax > mxdsml) then
         write(6,'("  stopped in ditspc_start.  Increase the value ",    &
     &       "of mxdsml.  Suggested value: ", i7)') imax

         stop

       endif

       if(ipr > 2) write(6,'("  initial values:  nexact = ",i6,          &
     &      "  niter = ",i6,"  ndeg = ",i6)') nexact,niter,ndeg

!      first nexact basis vectors are from small matrix diagonalization.
!      calculate hbas=hamk*bas. ihkp is a counter for vector matrix products.

       do n = 1,nexact
         call zcopy(mtxds,tvec(1,n),1,bas(1,n),1)
         if(mtxds .lt. mtxd) then
           do j = mtxds+1,mtxd
             bas(j,n) = C_ZERO
           enddo
         endif
       enddo

       deallocate(tvec)
       deallocate(et)

       call hk_psi_c16(mtxd,nexact,bas,hbas,lnewanl,                     &
     & ng,kgv,                                                           &
     & ekpg,isort,vscr,kmscr,                                            &
     & anlga,xnlkb,nanl,                                                 &
     & mxddim,mxdsml,mxdanl,mxdgve,mxdscr)

       ihkp = nexact

!      guess eigenvectors. uses either the results of the small
!      matrix diagonalization or (neig) user supplied guesses.
!      nitold user supplied guesses are added to the iterative basis.

       if(iguess /= 1) then
         do n = 1,niter
           call zcopy(mtxd,bas(1,n),1,psi(1,n),1)
           call zcopy(mtxd,hbas(1,n),1,hpsi(1,n),1)
         enddo
         lortho = .TRUE.
         nitold = 0
       else

         call hk_psi_c16(mtxd,neig,psi,hpsi,lnewanl,                     &
     &   ng,kgv,                                                         &
     &   ekpg,isort,vscr,kmscr,                                          &
     &   anlga,xnlkb,nanl,                                               &
     &   mxddim,mxdbnd,mxdanl,mxdgve,mxdscr)

         ihkp = ihkp + neig
         do n = 1,neig
           call zcopy(mtxd,psi(1,n),1,bas(1,nexact+n),1)
           call zcopy(mtxd,hpsi(1,n),1,hbas(1,nexact+n),1)
         enddo
         lortho = .FALSE.
         nitold = neig
         if(niter > neig) then
           do n = neig+1,niter
             call zcopy(mtxd,bas(1,n),1,psi(1,n),1)
             call zcopy(mtxd,hbas(1,n),1,hpsi(1,n),1)
           enddo
         endif
       endif

       return
       end subroutine ditsp_c16_start
