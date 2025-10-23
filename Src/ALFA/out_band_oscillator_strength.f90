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

!>  calculates and prints the oscillator strengths
!>
!>  \author       Jose Luis Martins
!>  \version      5.12
!>  \date         2 July 2014.  23 October2025.
!>  \copyright    GNU Public License v2

subroutine out_band_oscillator_strength(neig, ei, dh0drk, adot,          &
               lpair, lexcit, lxyz, rdircar,                             &
               ninitbeg, ninitend, nfinalbeg, nfinalend,                 &
               mxdbnd)

! Written July 2, 2014. JLM
! Modified, 4 March 2020, documentation. JLM
! Modified, name, indentation, new grouping of excitation levels. 16 May 2024. JLM
! Modified to be more flexible. 14 May 2025. JLM
! Modified to allow values in a giben direction. 23 October2025. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands

  integer, intent(in)                ::  neig                            !<  number of wavefunctions
  real(REAL64), intent(in)           ::  ei(mxdbnd)                      !<  eigenvalues (Hartree) for rk0
  complex(REAL64), intent(in)        ::  dh0drk(mxdbnd,mxdbnd,3)         !<  d <Psi|H|Psi> d k
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in real space

  logical, intent(in)                ::  lpair                           !<  prints the oscillator strengths for pairs of bands
  logical, intent(in)                ::  lexcit                          !<  prints the oscillator strengths by excitation energies
  logical, intent(in)                ::  lxyz                            !<  prints osc. str. in x y z directions.  Otherwise in rdircar direction
  real(REAL64), intent(in)           ::  rdircar(3)                      !<  choice of direction (cartesian coordinates)

  integer, intent(in)                ::  ninitbeg, ninitend              !<  begin and end of initial state index
  integer, intent(in)                ::  nfinalbeg, nfinalend            !<  begin and end of final state index

! local allocatable arrays

  real(REAL64), allocatable          ::  efmei(:)                        !  E_f - E_i
  complex(REAL64), allocatable       ::  ematif(:)                       !  |<i|d H / d k|f>|**2
  complex(REAL64), allocatable       ::  ematcar(:,:)                    !  |<i|d H / d k_x|f>|**2 cartesian

  integer, allocatable               ::  levdeg(:)                       !  degeneracy of energy level
  integer, allocatable               ::  leveigs(:,:)                    !  states belonging to level

  integer, allocatable               ::  indx(:)                         !  index for sorting
  integer, allocatable               ::  invindx(:)                      !  inverse of index for sorting
  real(REAL64), allocatable          ::  efmeisort(:)                    !  sorted E_f - E_i in increasing order

  logical, allocatable               ::  lnodup(:)                       !  allows skipping duplicate pairs

! local variables

  integer           ::  nfinal, ninit                                    !  number of initial and final states
  real(REAL64)      ::  avec(3,3)                  !  primitive lattice vectors that generate adot in canonical orientation
  real(REAL64)      ::  bvec(3,3)                  !  reciprocal lattice vectors

  integer           ::  mxdlev                                           !  array dimension for number of levels
  integer           ::  mxddeg                                           !  array dimension for number of levels

  integer           ::  nlevel                     !  number of energy levels (not counting degeneracies)
  integer           ::  maxdeg                     !  maximum number of degeneracies
  integer           ::  neigin                     !  number of pairs of levels
  integer           ::  nlevpr                     !  number of levels that are printed

  real(REAL64)      ::  sm(3)                      !  sum of M^2
  real(REAL64)      ::  so                         !  sum of F_if

  real(REAL64)      ::  xsum
  real(REAL64)      ::  rdirnorm(3)                !  normalized rdircar
  real(REAL64)      ::  ematdirsq                  !  emat in a given direction squared

! constants

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter     ::  PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter     ::  HARTREE = 27.21138386_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
  real(REAL64), parameter     ::  TOL = 1.0E-5_REAL64

! counters

  integer    ::  i, j, k, n, m
  integer    ::  ij, ji, ii, jj


! checks input variables

  if(neig > mxdbnd) then
    write(6,*)
    write(6,*) '   Stopped in out_band_oscillator_strength'
    write(6,*) '   Number of bands', neig,' greater than dimensions', mxdbnd
    write(6,*)

    stop

  endif

  if(ninitbeg < 1 .or. ninitend > neig .or. ninitbeg > ninitend) then
    write(6,*)
    write(6,*) '   Stopped in out_band_oscillator_strength'
    write(6,*) '   Inconsistent ninitbeg, ninitend, neig', ninitbeg, ninitend, neig
    write(6,*)

    stop

  endif

  if(nfinalbeg < 1 .or. nfinalend > neig .or. nfinalbeg > nfinalend) then
    write(6,*)
    write(6,*) '   Stopped in out_band_oscillator_strength'
    write(6,*) '   Inconsistent nfinalbeg, nfinalend, neig', nfinalbeg, nfinalend, neig
    write(6,*)

    stop

  endif

  if(.not. lxyz) then
    xsum = ZERO
    do i = 1,3
      xsum = xsum + rdircar(i)*rdircar(i)
    enddo
    if(xsum > TOL*TOL*TOL) then
      do i = 1,3
        rdirnorm(i) = rdircar(i) / sqrt(xsum)
      enddo
    else
      write(6,*)
      write(6,*) '   Small or zero direction vector'
      write(6,*) '   using x-direction'
      write(6,*)
      rdirnorm(1) = UM
      rdirnorm(2) = ZERO
      rdirnorm(3) = ZERO
    endif
  endif


  call adot_to_avec_sym(adot,avec,bvec)

  write(6,*)
  write(6,*) '  Dipole matrix elements and Oscillator strengths '
  write(6,*)
  write(6,*) '  The cartesian matrix elements are defined'
  write(6,*) '  with respect to the following orientation of'
  write(6,*) '  the crystal primitive lattice vectors'
  write(6,*)
  write(6,'("  a1 = (",f9.4,")    a2 = (",f9.4,")    a3 = (",f9.4,")")') &
                 (avec(1,k),k=1,3)
  write(6,'("       (",f9.4,")         (",f9.4,")         (",f9.4,")")') &
                 (avec(2,k),k=1,3)
  write(6,'("       (",f9.4,")         (",f9.4,")         (",f9.4,")")') &
                 (avec(3,k),k=1,3)

  ninit = ninitend - ninitbeg + 1
  nfinal = nfinalend - nfinalbeg + 1

  allocate(efmei(ninit*nfinal))
  allocate(ematif(ninit*nfinal))
  allocate(ematcar(ninit*nfinal,3))
  allocate(lnodup(ninit*nfinal))

  do i = 1,ninit
  do j = 1,nfinal
    ij = (i-1)*nfinal + j
    ii = ninitbeg + i - 1
    jj = nfinalbeg + j - 1
    efmei(ij) = ei(jj) - ei(ii)
    ematif(ij) = C_ZERO
    lnodup(ij) = .TRUE.
    do n = 1,3
    do m = 1,3
      ematif(ij) = ematif(ij) + dh0drk(jj,ii,n)*adot(n,m)*dh0drk(ii,jj,m)
    enddo
    enddo
    ematif(ij) = ematif(ij) / (4*PI*PI)
    do k = 1,3
      ematcar(ij,k) = C_ZERO
    enddo
    do n = 1,3
    do k = 1,3
      ematcar(ij,k) = ematcar(ij,k) + avec(k,n)*dh0drk(jj,ii,n)
    enddo
    enddo
    do k = 1,3
      ematcar(ij,k) = ematcar(ij,k) / (2*PI)
    enddo
  enddo
  enddo

  if(lpair) then

    write(6,*)
    write(6,*) '  Oscillator strengths by band pairs'
    write(6,*)

    write(6,*)
    if(lxyz) then
      write(6,'(4x,"i",3x,"E_i(eV)",4x,"f",3x,"E_f(eV)",5x,              &
            &   "|<i|dH/dk|f>|**2",8x,"F_if",9x,"|M_x|**2",4x,           &
            &   "|M_y|**2",4x,"|M_z|**2")')
    else
      write(6,'(4x,"i",3x,"E_i(eV)",4x,"f",3x,"E_f(eV)",5x,              &
            &   "|<i|dH/dk|f>|**2",8x,"F_if",9x,"|M_dir|**2")')
    endif
    write(6,*)

    do i = 1,ninit
    do j = 1,nfinal
      ij = (i-1)*nfinal + j
      ji = (j-1)*nfinal + i
      if(lnodup(ij)) then
        ii = ninitbeg + i - 1
        jj = nfinalbeg + j - 1
        lnodup(ji) = .FALSE.
        if(abs(efmei(ij)) > TOL) then
          if(lxyz) then
            write(6,'(i5,f10.3,i5,f10.3,6x,f12.5,4x,f12.5,4x,3f12.5)')   &
                ii, ei(ii)*HARTREE, jj, ei(jj)*HARTREE,                  &
                real(ematif(ij)),                                        &
                (2.0/3.0)*real(ematif(ij))/efmei(ij),                    &
                (abs(ematcar(ij,k))**2,k=1,3)
          else
            ematdirsq = abs(ematcar(ij,1)*rdirnorm(1) +                  &
                            ematcar(ij,2)*rdirnorm(2) +                  &
                            ematcar(ij,3)*rdirnorm(3))
            ematdirsq = ematdirsq*ematdirsq
            write(6,'(i5,f10.3,i5,f10.3,6x,f12.5,4x,f12.5,4x,3f12.5)')   &
                ii, ei(ii)*HARTREE, jj, ei(jj)*HARTREE,                  &
                real(ematif(ij)),                                        &
                (2.0/3.0)*real(ematif(ij))/efmei(ij),                    &
                ematdirsq
          endif
        endif
      endif
    enddo
    write(6,*)
    enddo
    write(6,*)

  endif

! now finds degeneracies but first sorts the excitation levels

  if(lexcit) then
    allocate(indx(ninit*nfinal))
    allocate(invindx(ninit*nfinal))
    allocate(efmeisort(ninit*nfinal))

    call sort(ninit*nfinal, efmei, indx)

    do ij = 1,ninit*nfinal
      efmeisort(ij) = efmei(indx(ij))
    enddo
    do ij = 1,ninit*nfinal
      invindx(indx(ij)) = ij
    enddo

    neigin = ninit*nfinal

!   number of levels that we know there no other levels with smaller excitations

    nlevpr = 1

    do i = 1,neigin
      if( efmeisort(i) > ei(neig) - ei(ninit) ) exit
      nlevpr = i
    enddo

    allocate(levdeg(1))
    allocate(leveigs(1,1))

    call berry_degeneracy(.TRUE., neigin, nlevpr, efmeisort, TOL,          &
           nlevel, maxdeg, levdeg, leveigs,                                &
           ninit*nfinal, 1, 1)

    mxdlev = nlevel
    mxddeg = maxdeg

    deallocate(levdeg)
    deallocate(leveigs)

    allocate(levdeg(mxdlev))
    allocate(leveigs(mxdlev,mxddeg))

!   fills the information

    call berry_degeneracy(.FALSE., neigin, nlevpr, efmeisort, TOL,         &
           nlevel, maxdeg, levdeg, leveigs,                                &
           ninit*nfinal, mxdlev, mxddeg)

    write(6,*)
    write(6,*) '  Oscillator strengths by excitation energies'
    write(6,*)

    do n = 1,nlevel

     if(abs(efmei(indx(leveigs(n,1)))) > TOL) then

        write(6,*)
        write(6,'("  excitation:",i5,6x,"degeneracy:",i5,6x,"energy(eV):"f12.5)')   &
              n, levdeg(n), efmei(indx(leveigs(n,1)))*HARTREE
        write(6,*)
        if(lxyz) then
          write(6,'(1x,"E_i-E_f(eV)",2x,"i",3x,"E_i(eV)",4x,"f",3x,"E_f(eV)",3x,    &
              &   "|<i|dH/dk|f>|**2",6x,"F_if",9x,                                  &
              &   "|M_x|**2",4x,"|M_y|**2",4x,"|M_z|**2",9x,                        &
              &   "|Fij_x|**2",2x,"|Fij_y|**2",2x,"|Fij_z|**2")')
        else
          write(6,'(1x,"E_i-E_f(eV)",2x,"i",3x,"E_i(eV)",4x,"f",3x,"E_f(eV)",3x,    &
              &   "|<i|dH/dk|f>|**2",6x,"F_if",9x,                                  &
              &   "|M_dir|**2",6x,"|Fij_dir|**2")')
        endif
        write(6,*)
        do k = 1,3
          sm(k) = ZERO
        enddo
        so = ZERO
        xsum = ZERO

        do m = 1,levdeg(n)

          ij = indx(leveigs(n,m))
          i = (ij-1) / nfinal + 1
          j = ij - (i-1)*nfinal
          ii = ninitbeg + i - 1
          jj = nfinalbeg + j - 1
          if(lxyz) then
            do k = 1,3
              sm(k) = sm(k) + abs(ematcar(ij,k))**2
            enddo
          else
            ematdirsq = abs(ematcar(ij,1)*rdirnorm(1) +                    &
                            ematcar(ij,2)*rdirnorm(2) +                    &
                            ematcar(ij,3)*rdirnorm(3))
            ematdirsq = ematdirsq*ematdirsq
            xsum = xsum + ematdirsq
          endif
          so = so + (2.0/3.0)*real(ematif(ij))/efmei(ij)

          if(lxyz) then
            write(6,'(f10.3,i5,f10.3,i5,f10.3,2x,f12.5,4x,f12.5,4x,        &
                &   3f12.5,4x,3f12.5)')                                    &
              efmei(ij)*HARTREE, ii, ei(ii)*HARTREE, jj, ei(jj)*HARTREE,   &
              real(ematif(ij)), (2.0/3.0)*real(ematif(ij))/efmei(ij),      &
              (abs(ematcar(ij,k))**2,k=1,3),                               &
              (abs(ematcar(ij,k))**2/efmei(ij),k=1,3)
          else
            write(6,'(f10.3,i5,f10.3,i5,f10.3,2x,f12.5,4x,f12.5,4x,        &
                &   f12.5,4x,f12.5)')                                      &
              efmei(ij)*HARTREE, ii, ei(ii)*HARTREE, jj, ei(jj)*HARTREE,   &
              real(ematif(ij)), (2.0/3.0)*real(ematif(ij))/efmei(ij),      &
              ematdirsq, ematdirsq/efmei(ij)
          endif

        enddo
        write(6,*)
        if(lxyz) then
          write(6,'(52x,"sum = ",f12.5,4x,3f12.5,4x,3f12.5)')              &
                so, (sm(k),k=1,3), (sm(k)/efmei(ij),k=1,3)
        else
          write(6,'(52x,"sum = ",f12.5,4x,f12.5,4x,f12.5)')                &
                so, xsum, xsum/efmei(ij)
        endif
        write(6,*)
        write(6,*)

      endif

    enddo


    deallocate(levdeg)
    deallocate(leveigs)

    deallocate(indx)
    deallocate(invindx)
    deallocate(efmeisort)

  endif

  deallocate(efmei)
  deallocate(ematif)
  deallocate(ematcar)
  deallocate(lnodup)

  return

end subroutine out_band_oscillator_strength
