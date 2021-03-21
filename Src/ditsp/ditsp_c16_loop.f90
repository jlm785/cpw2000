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

!>     Loop step of the iterative diagonalization

       subroutine ditsp_c16_loop(ipr,epspsi,epsa,                        &
     & lortho,lsafe,ihkp,lnewanl,                                        &
     & neig,mtxd,nexact,niter,ndeg,nconv,                                &
     & mcref,mc,xmax,                                                    &
     & psi,hpsi,bas,hbas,                                                &
     & eg,lconv,lexit,                                                   &
     & ekpg,isort,vscr,kmscr,                                            &
     & ng,kgv,                                                           &
     & anlga,xnlkb,nanl,                                                 &
     & hdiag,                                                            &
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
       real(REAL64), intent(in)           ::  epspsi                     !<  criteria for convergence
       real(REAL64), intent(in)           ::  epsa                       !<  eigenvalue energy scale
       integer, intent(in)                ::  mc                         !<  index of current inner iteration

       integer, intent(in)                ::  mtxd                       !<  dimension of the hamiltonian
       integer, intent(in)                ::  nexact                     !<  number of basis functions in bas that are never changed
       integer, intent(in)                ::  ndeg                       !<  number of eigenvectors including degeneracy

       real(REAL64), intent(in)           ::  ekpg(mxddim)                    !<  kinetic energy (hartree) of k+g-vector of row/column i
       integer, intent(in)                ::  isort(mxddim)              !<  g-vector associated with row/column i of hamiltonian
       real(REAL64), intent(in)           ::  vscr(mxdscr)               !<  screened potential in the fft real space mesh
       integer, intent(in)                ::  kmscr(7)                   !<  max value of kgv(i,n) used for the potential fft mesh and fft mesh size

       real(REAL64), intent(in)           ::  hdiag(mxddim)              !<  hamiltonian diagonal

       integer, intent(in)                ::  ng                         !<  total number of g-vectors with length less than gmax
       integer, intent(in)                ::  kgv(3,ng)                  !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

       complex(REAL64), intent(in)        ::  anlga(mxddim,mxdanl)       !<  KB projectors
       real(REAL64), intent(in)           ::  xnlkb(mxdanl)              !<  KB projector normalization
       integer, intent(in)                ::  nanl                       !<  half of number of projectors without spin

!      output

       integer, intent(out)               ::  neig                       !<  number of eigenvectors
       integer, intent(out)               ::  mcref                      !<  reference number of inner iterations
       real(REAL64),intent(out)           ::  xmax                       !<  maximum error: | H |psi> - |psi><psi|H|psi> |^2

       real(REAL64),intent(out)           ::  eg(mxddim)                 !<  eigenvalues CORRECT DIMENSIONS LATER
       logical,intent(out)                ::  lconv                      !<  indicates convergence
       logical,intent(out)                ::  lexit                      !<  indicates break from loop
       integer, intent(out)               ::  ihkp                       !<  number of H |psi> products

!      input and output

       logical, intent(inout)             ::  lortho                     !<  indicates that the |psi> are orthogonal
       logical, intent(inout)             ::  lsafe                      !<  indicates that no problems were identified in the iterative basis

       integer, intent(inout)             ::  niter                      !<  number of vectors in iterative space
       integer, intent(inout)             ::  nconv                      !<  number of converged vectors

       complex(REAL64), intent(inout)     ::  psi(mxddim,mxdbnd)         !<  |psi>
       complex(REAL64), intent(inout)     ::  hpsi(mxddim,mxdbnd)        !<  H |psi>
       complex(REAL64), intent(inout)     ::  bas(mxddim,mxdsml)         !<  |bas>
       complex(REAL64), intent(inout)     ::  hbas(mxddim,mxdsml)        !<  H |bas>

       logical, intent(inout)             ::  lnewanl                    !<  indicates that anlga has been recalculated (not used in default implementation)

!      local variables

       real(REAL64)  ::  xmold                                           !  old value of maximum error: | H |psi> - |psi><psi|H|psi> |^2
       integer       ::  nmax, ntst
       integer       ::  ired, nac
       real(REAL64)  ::  emax, emin, xn
       logical       ::  lxn

!      temps/work

       complex(REAL64),allocatable        ::  vsm(:,:)
       complex(REAL64),allocatable        ::  xerror(:)

       integer, allocatable               ::  irow(:)
       real(REAL64), allocatable          ::  work(:)


!      counters

       integer       :: i, n, nn

!      parameters

       real(REAL64), parameter :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
       complex(REAL64), parameter  ::  C_UM = cmplx(UM,ZERO,REAL64)
       complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)

!      external functions

       complex(REAL64),external   :: zdotc

!      allocations

       allocate(irow(mxdsml))
       allocate(work(mxdsml))

!      state at return from subroutine

       lconv = .FALSE.
       lexit = .FALSE.

       ihkp = 0

       allocate(xerror(mxddim))

       do n = nconv+1,niter
         eg(n) = real(zdotc(mtxd,psi(1,n),1,hpsi(1,n),1),REAL64)
         call zcopy(mtxd,hpsi(1,n),1,xerror,1)
         call zaxpy(mtxd,cmplx(-eg(n),ZERO,REAL64),                      &
     &                            psi(1,n),1,xerror,1)
         work(n) = real(zdotc(mtxd,xerror,1,xerror,1),REAL64)
       enddo

       deallocate(xerror)

       xmold = xmax
       xmax = zero
       nmax = min(ndeg,niter)
       do n = nconv+1,nmax
         if(work(n) > xmax) xmax = work(n)
       enddo
       if(ipr > 3) write(6,'("  In ditspc_loop: current error values",   &
     &                   " of error work(n) = ",/,20(2X,5E13.5,/))')     &
     &                  (work(n),n=nconv+1,nmax)

!      if the vectors are not improving exits relaxation loop

       if(mc > 1 .and. xmax > 0.2*xmold) then
         mcref = mc
         if(ipr > 3) write(6,*) '  exiting ditspc_loop 1'

         lexit = .TRUE.

       endif

!      test for global convergence (|(h-e)*x|<epspsi for ndeg vectors)
!      lortho indicates that eigenvectors are orthogonal,
!      lsafe indicates that no linear dependencies were found
!      since the previous diagonalization step.

       if(.not. lexit) then
         if(xmax < epspsi*epspsi .and. lsafe .and. niter >= ndeg) then

           if(lortho) then

             neig = ndeg
             if(ipr > 2) write(6,*) '  returning ditspc_loop 2'

             lconv = .TRUE.

           else

!            orthogonalizes before retesting and
!            removes linear dependent vectors from iteration


             call grsch_loop_c16(.TRUE.,psi,hpsi,mtxd,nconv,niter,irow,     &
     &       mxddim)

             call grsch_shift_c16(psi,mtxd,nconv,niter,irow,ired,        &
     &       mxddim)

             call grsch_shift_c16(hpsi,mtxd,nconv,niter,irow,ired,       &
     &       mxddim)

             if(ired > 0) lsafe = .FALSE.

             n = nconv+1
             do nac = nconv+1,niter
               if(irow(nac) /= 0) then
                 if(n < nac) then
                   eg(n) = eg(nac)
                   work(n) = work(nac)
                 endif
                 n = n + 1
               endif
             enddo

             if(ipr > 3 .and. ired > 0) write(6,'("  In ditspc_loop ",   &
     &      "the niter orthogonalization found",i6," linearly ",         &
     &      "dependent eigenvectors")') ired

             niter = niter-ired
             lortho = .TRUE.

             if(nconv < niter) then

               xmax = zero
               nmax = min(niter,ndeg)
               do n = nconv+1,nmax
                 if(work(n) > xmax) xmax = work(n)
               enddo
               if(xmax < epspsi*epspsi .and. lsafe .and.                 &
     &                                     niter >= ndeg) then
                 neig = ndeg
                 if(ipr > 2) write(6,*) '  returning  ditspc_loop 3'

                 lconv = .TRUE.

               endif
             endif

           endif                                                         ! if(lortho)

         endif
       endif

!      end of global convergence test
!      break out of inner loop

       if(.not. lconv .and. (.not. lexit) ) then
         if(niter == nconv) then
           if(ipr > 3) write(6,*) '  exiting ditspc_loop 4'

           lexit = .TRUE.

         endif
       endif

!      test for individual convergence.
!      degenerate eigenvectors must be accepted together,
!      and accepted eigenvalues must be the lowest.
!      ntst keeps track of candidate eigenvectors (ntst-nconv)


       if(.not. lconv .and. (.not. lexit) ) then

         ntst = nconv
         do n = nconv+1,niter

           if(work(n) > 0.01*epspsi*epspsi) exit

           ntst = n
         enddo

         if(ntst > nconv) then

!          adjusts ntst for degeneracies loop (i) is a while...

           do i = nconv,ntst
             emax = eg(nconv+1)
             do n = nconv+1,ntst
               if(eg(n) > emax) emax = eg(n)
             enddo
             emin = emax + 10.0*epsa
             if(ntst < niter) then
               do n = ntst+1,niter
                 if(eg(n) < emin) emin = eg(n)
               enddo
             endif
             if(emax + epsa > emin) then
               ntst = ntst - 1

               if(ntst <= nconv) exit

             endif
           enddo
         endif

       endif

!      if there are candidate vectors and they are not orthogonal
!      they are orthogonalized

       if(.not. lconv .and. (.not. lexit) ) then
         if(ntst > nconv .and. .not. lortho) then

           call grsch_loop_c16(.TRUE.,psi,hpsi,mtxd,nconv,ntst,irow,        &
     &     mxddim)

           do n = ntst+1,niter
             irow(n) = 1
           enddo

           call grsch_shift_c16(psi,mtxd,nconv,niter,irow,ired,          &
     &     mxddim)

           call grsch_shift_c16(hpsi,mtxd,nconv,niter,irow,ired,         &
     &     mxddim)

           if(ired > 0) lsafe = .FALSE.

           n = nconv+1
           do nac = nconv+1,niter
             if(irow(nac) /= 0) then
               if(n < nac) then
                 eg(n) = eg(nac)
                 work(n) = work(nac)
               endif
               n = n + 1
             endif
           enddo

           if(ipr > 3 .and. ired > 0) write(6,'("  In ditspc_loop the ", &
     &       "ntst orthogonalization found",i6," linearly dependent ",   &
     &          " eigenvectors")') ired
           lortho = .TRUE.
           ntst = ntst - ired
           niter = niter - ired

!          retests the eigensolutions

           if(ntst > nconv) then
             nmax = ntst
             ntst = nconv
             do n = nconv+1,nmax

               if(work(n) > 0.01*epspsi*epspsi) exit

               ntst = n
             enddo

             if(ntst > nconv) then
               do i = nconv,ntst
                 emax = eg(nconv+1)
                 do n = nconv+1,ntst
                   if(eg(n) > emax) emax = eg(n)
                 enddo
                 emin = emax + 10.0*epsa
                 if(ntst < niter) then
                   do n = ntst+1,niter
                     if(eg(n) < emin) emin = eg(n)
                   enddo
                 endif
                 if(emax + epsa > emin) then
                   ntst = ntst - 1

                   if(ntst <= nconv) exit

                 endif
               enddo
             endif
           endif

         endif

       endif

!      break out of inner loop


       if(.not. lconv .and. (.not. lexit) ) then
         if(niter == nconv) then
           if(ipr > 3) write(6,*) '  exiting ditspc_loop 3'

           lexit = .TRUE.

         endif
       endif

!      at this point there are ntst-nconv candidates.
!      spanned subspace should be similar for replacement
!      of vectors in the iterative subspace. (avoids
!      missing eigenvectors.)

       if(.not. lconv .and. (.not. lexit) ) then

         if(ntst > nconv) then

           allocate(vsm(mxdsml,mxdsml))

           call zgemm('c','n',ntst-nconv,ntst-nconv,mtxd,C_UM,           &
     &                bas(1,nexact+nconv+1),mxddim,                      &
     &                psi(1,nconv+1),mxddim,C_ZERO,vsm,mxdsml)

!          span test. if successful ntst does not change

           nn = 0
           do i = 1,10*mxdbnd
             nmax = ntst
             ntst = nconv
             lxn = .FALSE.
             do n = 1,nmax-nconv
               xn = real(zdotc(nmax-nconv,vsm(1,n),1,vsm(1,n),1),REAL64)

               if(xn < 0.95) then
                 lxn = .TRUE.

                 exit

               else
                 ntst = nconv + n
               endif
             enddo

             if(lxn) then

               if(ntst <= nconv) exit

             else

!              degeneracy check

               emax = eg(nconv+1)
               do n = nconv+1,ntst
                 if(eg(n) > emax) emax = eg(n)
               enddo
               emin = emax + 10.0*epsa
               if(ntst < niter) then
                 do n = ntst+1,niter
                   if(eg(n) < emin) emin = eg(n)
                 enddo
               endif

               if(emax + epsa > emin) then
                 ntst = ntst - 1

                 if(ntst <= nconv) exit

               else

                 exit

               endif
             endif

             nn = i
           enddo

           if(nn == 10*mxdbnd) then
             write(6,*) '  stopped in old 72 loop of ditspc_loop'

             stop

           endif

           deallocate(vsm)

         endif

       endif

!      replaces basis vectors

       if(.not. lconv .and. (.not. lexit) ) then
         if(ntst > nconv) then

           do n = nconv+1,ntst
             call zcopy(mtxd,psi(1,n),1,bas(1,nexact+n),1)
           enddo

           call grsch_loop_c16(.FALSE.,bas,hbas,mtxd,nexact+nconv,           &
     &     nexact+ntst,irow,                                             &
     &     mxddim)

           call grsch_shift_c16(bas,mtxd,nexact+nconv,nexact+ntst,       &
     &     irow,ired,                                                    &
     &     mxddim)

!          removes the corresponding test vector

           if(ired > 0) then

             lsafe = .FALSE.

             do n = 1,nconv
               irow(n) = 1
             enddo
             do n = nconv+1,ntst
               irow(n) = irow(nexact+n)
             enddo
             do n = ntst+1,niter
               irow(n) = 1
             enddo

             call grsch_shift_c16(psi,mtxd,nconv,niter,irow,ired,        &
     &       mxddim)

             call grsch_shift_c16(hpsi,mtxd,nconv,niter,irow,ired,       &
     &       mxddim)

             n = nconv+1
             do nac = nconv+1,niter
               if(irow(nac) /= 0) then
                 if(n < nac) then
                   eg(n) = eg(nac)
                   work(n) = work(nac)
                 endif
                 n = n + 1
               endif
             enddo

           endif


           call hk_psi_c16(mtxd,ntst-ired-nconv,                         &
     &     bas(1,nexact+nconv+1),hbas(1,nexact+nconv+1),lnewanl,         &
     &     ng,kgv,                                                       &
     &     ekpg,isort,vscr,kmscr,                                        &
     &     anlga,xnlkb,nanl,                                             &
     &     mxddim,mxdbnd,mxdanl,mxdgve,mxdscr)

           ihkp = ihkp + ntst - ired - nconv
           nconv = ntst - ired
           niter = niter - ired
         endif

!        relaxation step for not yet converged eigenvectors

         do n = nconv+1,niter
           call zaxpy(mtxd,cmplx(-eg(n),ZERO,REAL64),                    &
     &                                  psi(1,n),1,hpsi(1,n),1)
         enddo

         call rq_jac_c16(hpsi(1,nconv+1),eg(nconv+1),hdiag,mtxd,         &
     &   niter-nconv,                                                    &
     &   mxddim,mxdbnd)

         do n = nconv+1,niter
           call zaxpy(mtxd,C_UM,hpsi(1,n),1,psi(1,n),1)
           xn = real(zdotc(mtxd,psi(1,n),1,psi(1,n),1),REAL64)
           xn = UM/sqrt(xn)
           call zscal(mtxd,cmplx(xn,ZERO,REAL64),psi(1,n),1)
         enddo

         lortho = .FALSE.

         call hk_psi_c16(mtxd,niter-nconv,                               &
     &   psi(1,nconv+1),hpsi(1,nconv+1),lnewanl,                         &
     &   ng,kgv,                                                         &
     &   ekpg,isort,vscr,kmscr,                                          &
     &   anlga,xnlkb,nanl,                                               &
     &   mxddim,mxdbnd,mxdanl,mxdgve,mxdscr)

         ihkp = ihkp + niter - nconv

       endif

       deallocate(irow)
       deallocate(work)

       return
       end subroutine ditsp_c16_loop
