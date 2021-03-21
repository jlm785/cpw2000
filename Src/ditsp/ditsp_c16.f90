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

!>     Iterative diagonalization of a matrix with a leading
!>     submatrix. It is inspired on the ritzit procedure
!>     (H. Ruthishauser Numer. Math. 16, 205 (1970))
!>     but uses Jacobi relaxation. Complex version.
!>     The non-local pseudopotential is separable.

       subroutine ditsp_c16(ipr,ifail,icmax,iguess,epspsi,lnewanl,       &
     & neig,mtxd,mtxds,                                                  &
     & psi,hpsi,ei,                                                      &
     & ekpg,isort,vscr,kmscr,                                            &
     & ng,kgv,                                                           &
     & anlga,xnlkb,nanl,                                                 &
     & hamsm,hdiag,                                                      &
     & mxddim,mxdsml,mxdbnd,mxdgve,mxdscr,mxdanl)

!      Written february 1990. JLM
!      Modified 17 january 2007 for new orthogonalization
!      and dgemm hamiltonian calculation. JLM
!      Split in three parts and full complex version,
!      12 August 2015. JLM
!      Modified kmscr, 28 October 2015. JLM
!      Modified, documentation, January 2020. JLM
!      Modified, ei dimensions finaly corrected, ifail,icmax. 18 February 2020. JLM
!      Modified, debugging print statements, 12 June 2020. JLM
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

       integer, intent(in)                ::  ipr                        !<  print control. 0 no printing, 2, 3, 4 lots of stuff

       integer, intent(in)                ::  icmax                      !<  maximum value of outer iteration
       integer, intent(in)                ::  iguess                     !<  if =1 uses neig guess eigenvectors in xvecr
       real(REAL64), intent(in)           ::  epspsi                     !<  criteria for convergence

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
       real(REAL64), intent(in)           ::  hdiag(mxddim)              !<  hamiltonian diagonal

!      output

       integer, intent(out)               ::  ifail                      !<  if ifail=0 the subroutine was successfull. Otherwise ifail indicates the number of correct digits.
       real(REAL64),intent(out)           ::  ei(mxdbnd)                 !<  eigenvalues

!      input and output

       complex(REAL64), intent(inout)     ::  psi(mxddim,mxdbnd)         !<  |psi>
       complex(REAL64), intent(inout)     ::  hpsi(mxddim,mxdbnd)        !<  H |psi>

       integer, intent(inout)             ::  neig                       !<  number of eigenvectors

       logical, intent(inout)             ::  lnewanl                    !<  indicates that anlga has been recalculated (not used in default implementation)

!      local variables

       integer                ::  mcref                                  !  reference number of inner iterations
       integer                ::  mcmax                                  !  maximum value of inner iteration
       real(REAL64)           ::  xmax                                   !  maximum error: | H |psi> - |psi><psi|H|psi> |^2

       logical                ::  lnoiter                                !  indicates that the eigensolution was obtained by conventional algorithms
       logical                ::  lconv                                  !  indicates convergence
       logical                ::  lexit                                  !  indicates break from loop
       logical                ::  lortho                                 !  indicates that the |psi> are orthogonal
       logical                ::  lsafe                                  !  indicates that no problems were identified in the iterative basis

       integer                ::  ihkp                                   !  number of H |psi> products
       integer                ::  ihkpsum                                !  total number of H |psi> products
       integer                ::  niter                                  !  number of vectors in iterative space
       integer                ::  nconv                                  !  number of converged vectors
       integer                ::  nexact                                 !  number of basis functions in bas that are never changed
       integer                ::  ndeg                                   !  number of eigenvectors including degeneracy
       integer                ::  iseed                                  !  random number generator seed
       integer                ::  nitold                                 !  previous number of vectors in iterative space
       real(REAL64)           ::  epsa                                   !  eigenvalue quasi-degeneracy criterium

       integer                ::  info

       real(REAL64)           ::  xmax2

!      allocatable arrays

       real(REAL64),allocatable           ::  eg(:)

       complex(REAL64),allocatable        ::  bas(:,:)
       complex(REAL64),allocatable        ::  hbas(:,:)
       complex(REAL64),allocatable        ::  vsm(:,:)

!      counters

       integer       :: i, n, ic, mc


!      initializations

       allocate(eg(mxddim))

       ifail = 0
       xmax = 1.0
       iseed = -1
       mcref = 2
       ihkpsum = 0

       lsafe = .FALSE.
       lortho = .FALSE.
       lconv = .FALSE.
       lexit = .FALSE.
       lnoiter = .FALSE.

!      tests for array sizes, etc...

       if(neig < 1) then
         write(6,'(/,"   WARNING:   ditsp_c16 called with number",       &
     &        " of eigenvalues = ",i8)') neig

         return

       endif

       if(neig > mtxds) then
         write(6,'("  Stopped in ditsp_c16. neig = ",i8," is greater",   &
     &        " than mtxds = ",i8)')

         stop

       endif

       if(neig > mxdbnd) then
         write(6,'("  Stopped in ditsp_c16. neig = ",i8," is greater",   &
     &        " than mxdbnd = ",i8)')

         stop

       endif

       if(3*neig > mxdsml) then
         write(6,'("  Stopped in ditsp_c16. 3*neig = ",i8," is greater", &
     &        " than mxdsml = ",i8)')

         stop

       endif

       if(mtxds > mxdsml) then
         write(6,'("  Stopped in ditsp_c16. mtxds = ",i8," is greater",  &
     &        " than mxdsml = ",i8)')

         stop

       endif

       if(mtxd > mxddim) then
         write(6,'("  Stopped in ditsp_c16. mtxd = ",i8," is greater",   &
     &        " than mxddim = ",i8)')

         stop

       endif

       if(nanl > mxdanl) then
         write(6,'("  Stopped in ditsp_c16. nanl = ",i8," is greater",   &
     &        " than mxdanl = ",i8)')

         stop

       endif

       if(ng > mxdgve) then
         write(6,'("  Stopped in ditsp_c16. ng = ",i8," is greater",     &
     &        " than mxdgve = ",i8)')

         stop

       endif

!      tests for non-iterative diagonalization

       if(mtxds >= mtxd) then

         allocate(vsm(mxdsml,mxdsml))

         call diag_c16(mtxds,hamsm,eg,vsm,mxdsml,info)


         if(info /= 0) then
           write(6,*)
           write(6,*) '  ditsp_c16 failure in non-iterative diag.'
           write(6,*)

           stop

         endif

         do n = 1,neig
         do i = 1,mtxd
           psi(i,n) = vsm(i,n)
         enddo
         enddo

         lnoiter = .TRUE.

         deallocate(vsm)

       else

         if(iguess /= 0.and. ipr > 2) then
           call ditsp_error(xmax2, neig, mtxd, psi, hpsi,                   &
     &       mxddim, mxdbnd)
           write(6,*)
           write(6,*) '  From the guess eigenvectors '
           write(6,*) '  The estimated error in energy has an accuracy'
           write(6,'("  of ",f8.1," digits")') -log10(xmax2)
           write(6,*)
         endif

         allocate(bas(mxddim,mxdsml))
         allocate(hbas(mxddim,mxdsml))

         call ditsp_c16_start(ipr,iguess,epspsi,epsa,lortho,ihkp,lnewanl,   &
     &   neig,mtxd,mtxds,nexact,niter,nitold,ndeg,                          &
     &   psi,hpsi,bas,hbas,                                                 &
     &   ekpg,isort,vscr,kmscr,                                             &
     &   ng,kgv,                                                            &
     &   anlga,xnlkb,nanl,                                                  &
     &   hamsm,                                                             &
     &   mxddim,mxdsml,mxdbnd,mxdgve,mxdscr,mxdanl)

         ihkpsum = ihkpsum + ihkp

         if(ipr > 2) then
           call ditsp_error(xmax2, neig, mtxd, psi, hpsi,                   &
     &      mxddim, mxdbnd)
           write(6,*)
           write(6,*) '  After ditsp_c16_start with iguess = ',iguess
           write(6,*) '  The estimated error in energy has an accuracy'
           write(6,'("  of ",f8.1," digits")') -log10(xmax2)
           write(6,*)
         endif

!        initial value of nconv is 0. nconv is the number of
!        eigenvectors that are considered converged (|(h-e)*x|<0.1*epspsi,
!        and all eigenvectors with lower energy are also converged).
!        it may be increased in the relaxation loop and
!        decreased in the main loop.

         nconv = 0

!        niter is recalculated at
!        the end of the main (outer) loop, and reduced when
!        linear dependencies appear in the iteration scheme.
!        niter is preferably greater than ndeg.

!        ndeg is modified
!        at the end of the main (outer) loop. before return neig
!        is changed to ndeg. ndeg is never smaller than input neig.

!        Begin main loop

!        contains the relaxation loop and the diagonalization
!        of the hamiltonian in the iterative basis.

         do ic = 1,icmax
           mcmax = mcref

!          Begin relaxation loop

!          evaluates the errors, decides if relaxation is improving,
!          checks for convergence, and finnally calculates new guesses.

           do mc=1,mcmax
             if(ipr > 3) write(6,'("  ic = ",i6,"  mc = ",i6,            &
     &         "  nconv = ",i6,"  niter = ",i6)') ic,mc,nconv,niter

             call ditsp_c16_loop(ipr,epspsi,epsa,                        &
     &       lortho,lsafe,ihkp,lnewanl,                                  &
     &       neig,mtxd,nexact,niter,ndeg,nconv,                          &
     &       mcref,mc,xmax,                                              &
     &       psi,hpsi,bas,hbas,                                          &
     &       eg,lconv,lexit,                                             &
     &       ekpg,isort,vscr,kmscr,                                      &
     &       ng,kgv,                                                     &
     &       anlga,xnlkb,nanl,                                           &
     &       hdiag,                                                      &
     &       mxddim,mxdsml,mxdbnd,mxdgve,mxdscr,mxdanl)

             ihkpsum = ihkpsum + ihkp

             if(ipr > 3) then
               call ditsp_error(xmax2, neig, mtxd, psi, hpsi,            &
    &          mxddim, mxdbnd)
               write(6,*)
               write(6,*) '  After ditsp_c16_loop with mc = ',mc
               write(6,*) '  The estimated error in energy has an accuracy'
               write(6,'("  of ",f8.1," digits")') -log10(xmax2)
               write(6,*)
             endif

             if(lexit .or. lconv) exit

           enddo

           if(lconv) then
             if(ipr > 2) then
               call ditsp_error(xmax2, neig, mtxd, psi, hpsi,            &
    &          mxddim, mxdbnd)
               write(6,*)
               write(6,*) '  After ditsp_c16_loop convergence '
               write(6,*) '  The estimated error in energy has an accuracy'
               write(6,'("  of ",f8.1," digits")') -log10(xmax2)
               write(6,*)
             endif

           endif

           if(.not. lsafe) mcref = 3

           call ditsp_c16_ritz(lsafe,ipr,ihkp,epsa,iseed,lnewanl,        &
     &     mtxd,neig,nexact,nitold,niter,nconv,                          &
     &     psi,hpsi,bas,hbas,eg,lconv,                                   &
     &     ekpg,isort,vscr,kmscr,                                        &
     &     ng,kgv,                                                       &
     &     anlga,xnlkb,nanl,                                             &
     &     mxddim,mxdsml,mxdbnd,mxdgve,mxdscr,mxdanl)

           call ditsp_error(xmax, neig, mtxd, psi, hpsi,              &
    &      mxddim, mxdbnd)

           if(ipr > 2) then
             write(6,*)
             write(6,*) '  After ditsp_c16_ritz with ic = ',ic
             write(6,*) '  The estimated error in energy has an accuracy'
             write(6,'("  of ",f8.1," digits")') -log10(xmax)
             write(6,*)
           endif

           ihkpsum = ihkpsum + ihkp

           lconv = .FALSE.
           if(xmax < epspsi*epspsi) lconv = .TRUE.

           if(lconv) exit

           nitold = niter
           lortho = .TRUE.
           lsafe = .TRUE.

           if(xmax < 0.01) mcref = min(mcref + 2, 8)


         enddo

         deallocate(bas)
         deallocate(hbas)

       endif

       if(lnoiter) then
         if(ipr > 1) then
           write(6,*)
           write(6,'("  Using standard diagonalization  mtxd = ",i8)')   &
     &          mtxd
           write(6,*)
         endif
       elseif(lconv) then
         if(ipr > 1) then
           write(6,*)
           write(6,'("  Successful diagonalization after",i8,            &
     &           " matrix-vector multiplications")') ihkpsum
           write(6,*)
         endif
       else
         ifail = 1+int(-log10(xmax))
         if(ifail == 0) ifail = 1
         if(ipr > 0) then
           write(6,*)
           write(6,'("  Unsuccessful diagonalization after",i8,          &
     &          " matrix-vector multiplications,  xmax = ",e12.4)')      &
     &          ihkpsum,xmax
           write(6,*)
         endif
       endif

       do i = 1,neig
         ei(i) = eg(i)
       enddo

       deallocate(eg)

       return
       end subroutine ditsp_c16
