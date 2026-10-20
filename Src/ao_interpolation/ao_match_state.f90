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

!>  Matches the eigenstates of an atomic orbital calculation
!>  with those of a plane-wave calculation at the same k-point,
!>  and gives each atomic orbital state the energy of its
!>  plane-wave partner (modified tight-binding).
!>
!   The mapping is not symmetric.  The atomic orbital states are the
!   targets, and an atomic orbital level (states with the same energy)
!   only gets plane-wave energies if all its states are well
!   represented by plane-wave states.  Otherwise it keeps its atomic
!   orbital energy.  Plane-wave states without a partner are ignored.
!
!   The weight of an atomic orbital level in a plane-wave level is the
!   sum of the |<pw|ao>|**2 of their states.  It does not depend on the
!   choice of basis inside degenerate levels.  The assignment of states
!   maximizes the sum of the weights.
!
!   Levels that mix (weight per atomic orbital state above XLINK) form
!   clusters, for example near-degenerate levels mixed by spin-orbit,
!   or a degenerate plane-wave level split by numerical noise.
!   An atomic orbital level is accepted if its weight per state in the
!   plane-wave levels of its cluster is at least xmin, and its partners
!   are in the same cluster.  All the states of an atomic orbital
!   level get the same energy, the average of their partners.
!
!   The last plane-wave level is not used, as the number of
!   plane-wave states may have cut a degenerate level.
!
!   It can be used with spinor wave-functions (spin-orbit),
!   with mtxd and mxddim replaced by 2*mtxd and 2*mxddim.
!>
!>  \author       Jose Luis Martins
!>  \version      5.13
!>  \date         29 October 2024, 19 October 2026.
!>  \copyright    GNU Public License v2

subroutine ao_match_state(ipr, xmin, mtxd,                               &
      neig_pw, psi_pw, ei_pw, neig_ao, psi_ao, ei_ao,                    &
      ei_match, imatch, xover,                                           &
      mxddim, mxdbnd_pw, mxdbnd_ao)


! Adapted from out_band_match_state, 29 October 2024. JLM
! Constants updated to CODATA 2022 (HARTREE, EV replaced by HARTREE). 3 October 2026. JLM+claude
! Matching of degenerate levels, assignment that maximizes the overlap (num_assign_max),
! ei_match, ipr, xmin. 19 October 2026. JLM+claude WARNING NEW API


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd_pw                       !<  array dimension for number of plane-wave bands
  integer, intent(in)                ::  mxdbnd_ao                       !<  array dimension for number of atomic orbital bands

  integer, intent(in)                ::  ipr                             !<  print control: 0 nothing, 1 warnings, 2 table of overlaps
  real(REAL64), intent(in)           ::  xmin                            !<  minimum weight per state of an atomic orbital level in its cluster (1/2 recommended)

  integer, intent(in)                ::  mtxd                            !<  dimension of the hamiltonian (basis set size)

  integer, intent(in)                ::  neig_pw                         !<  number of plane-wave eigenvectors
  complex(REAL64), intent(in)        ::  psi_pw(mxddim, mxdbnd_pw)       !<  component j of plane-wave eigenvectors
  real(REAL64), intent(in)           ::  ei_pw(mxdbnd_pw)                !<  eigenvalue of plane-wave eigenvectors (hartree)

  integer, intent(in)                ::  neig_ao                         !<  number of atomic orbital eigenvectors
  complex(REAL64), intent(in)        ::  psi_ao(mxddim, mxdbnd_ao)       !<  component j of atomic orbital eigenvectors
  real(REAL64), intent(in)           ::  ei_ao(mxdbnd_ao)                !<  eigenvalue of atomic orbital eigenvectors (hartree)

! output

  real(REAL64), intent(out)          ::  ei_match(mxdbnd_ao)             !<  energy for atomic orbital state i, from plane-waves if matched (hartree)
  integer, intent(out)               ::  imatch(mxdbnd_ao)               !<  plane-wave partner of atomic orbital state i, 0 if not matched
  real(REAL64), intent(out)          ::  xover(mxdbnd_ao)                !<  weight per state of the level of state i in the plane-wave levels of its cluster

! local allocatable arrays

  complex(REAL64), allocatable       ::  prod(:,:)                       !  <pw|ao>
  real(REAL64), allocatable          ::  prodsq(:,:)                     !  |<pw|ao>|**2
  real(REAL64), allocatable          ::  aonorm(:)                       !  norm of plane-wave state in the space of atomic orbital states

  integer, allocatable               ::  levdeg_pw(:), levdeg_ao(:)      !  degeneracy of energy level
  integer, allocatable               ::  leveigs_pw(:,:)                 !  states belonging to level
  integer, allocatable               ::  leveigs_ao(:,:)                 !  states belonging to level
  integer, allocatable               ::  ilev_pw(:), ilev_ao(:)          !  level of state

  real(REAL64), allocatable          ::  wlev(:,:)                       !  weight of atomic orbital level in plane-wave level
  real(REAL64), allocatable          ::  score(:,:)                      !  weight of atomic orbital state in level of plane-wave state
  integer, allocatable               ::  iassign(:)                      !  plane-wave state assigned to atomic orbital state
  integer, allocatable               ::  iprt(:)                         !  plane-wave states printed in the table
  integer, allocatable               ::  iclus(:)                        !  cluster (lowest node) of atomic orbital or plane-wave level

! local variables

  integer           ::  nlev_pw, nlev_ao                                 !  number of energy levels
  integer           ::  maxdeg_pw, maxdeg_ao                             !  maximum degeneracy
  integer           ::  neig                                             !  dummy for berry_degeneracy
  integer           ::  nused                                            !  number of plane-wave states used in the matching
  integer           ::  nlevu                                            !  number of plane-wave levels used in the matching
  integer           ::  ia, ip, la, lp
  logical           ::  lok                                              !  atomic orbital level is well matched
  real(REAL64)      ::  xsum, eavg

! constants

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter  ::  C_UM = cmplx(UM,ZERO,REAL64)
  real(REAL64), parameter     ::  EPS = 1.0E-6_REAL64                    !  tolerance for the normalization check
  real(REAL64), parameter     ::  TOL = 1.0E-5_REAL64                    !  tolerance for degenerate energies (hartree)
  real(REAL64), parameter     ::  XNORM = 0.9_REAL64                     !  warns if plane-wave partner is less than that in atomic orbital space
  real(REAL64), parameter     ::  XLINK = 0.1_REAL64                     !  minimum weight of atomic orbital state for levels in the same cluster
  real(REAL64), parameter     ::  XPRINT = 0.05_REAL64                   !  smallest overlap in the table
  real(REAL64), parameter     ::  HARTREE = 27.211386246_REAL64

! counters

  integer   :: i, j, n


! overlaps

  allocate(prod(neig_pw,neig_ao))
  allocate(prodsq(neig_pw,neig_ao))
  allocate(aonorm(neig_pw))

  call zgemm('C', 'N', neig_pw, neig_ao, mtxd, C_UM, psi_pw, mxddim,     &
          psi_ao, mxddim, C_ZERO, prod, neig_pw)

  do j = 1,neig_ao
  do i = 1,neig_pw
    prodsq(i,j) = real(prod(i,j)*conjg(prod(i,j)),REAL64)
  enddo
  enddo

  do i = 1,neig_pw
    aonorm(i) = ZERO
    do j = 1,neig_ao
      aonorm(i) = aonorm(i) + prodsq(i,j)
    enddo
  enddo

! paranoid check that things are consistent

  do j = 1,neig_ao
    xsum = ZERO
    do i = 1,neig_pw
      xsum = xsum + prodsq(i,j)
    enddo

    if(xsum > UM+EPS) then
      write(6,*)
      write(6,'("   STOPPED in ao_match_state: ")')
      write(6,'("   projection of atomic orbital state ",i5,             &
          &     " is greater than one: ",f14.8)') j, xsum
      write(6,*)

      stop

    endif

  enddo

! degenerate levels of the plane-wave states, first the dimensions

  allocate(levdeg_pw(1))
  allocate(leveigs_pw(1,1))

  neig = neig_pw
  call berry_degeneracy(.TRUE., neig_pw, neig, ei_pw, TOL,               &
         nlev_pw, maxdeg_pw, levdeg_pw, leveigs_pw,                      &
         mxdbnd_pw, 1, 1)

  deallocate(levdeg_pw)
  deallocate(leveigs_pw)

  allocate(levdeg_pw(nlev_pw))
  allocate(leveigs_pw(nlev_pw,maxdeg_pw))

  call berry_degeneracy(.FALSE., neig_pw, neig, ei_pw, TOL,              &
         nlev_pw, maxdeg_pw, levdeg_pw, leveigs_pw,                      &
         mxdbnd_pw, nlev_pw, maxdeg_pw)

! same for the atomic orbital states

  allocate(levdeg_ao(1))
  allocate(leveigs_ao(1,1))

  neig = neig_ao
  call berry_degeneracy(.TRUE., neig_ao, neig, ei_ao, TOL,               &
         nlev_ao, maxdeg_ao, levdeg_ao, leveigs_ao,                      &
         mxdbnd_ao, 1, 1)

  deallocate(levdeg_ao)
  deallocate(leveigs_ao)

  allocate(levdeg_ao(nlev_ao))
  allocate(leveigs_ao(nlev_ao,maxdeg_ao))

  call berry_degeneracy(.FALSE., neig_ao, neig, ei_ao, TOL,              &
         nlev_ao, maxdeg_ao, levdeg_ao, leveigs_ao,                      &
         mxdbnd_ao, nlev_ao, maxdeg_ao)

! level of each state

  allocate(ilev_pw(neig_pw))
  allocate(ilev_ao(neig_ao))

  do lp = 1,nlev_pw
    do i = 1,levdeg_pw(lp)
      ilev_pw(leveigs_pw(lp,i)) = lp
    enddo
  enddo

  do la = 1,nlev_ao
    do j = 1,levdeg_ao(la)
      ilev_ao(leveigs_ao(la,j)) = la
    enddo
  enddo

! the last plane-wave level may be incomplete

  nused = neig_pw
  if(nlev_pw > 1) nused = neig_pw - levdeg_pw(nlev_pw)

! weights of atomic orbital levels in plane-wave levels

  allocate(wlev(nlev_ao,nlev_pw))

  do lp = 1,nlev_pw
  do la = 1,nlev_ao
    wlev(la,lp) = ZERO
    do j = 1,levdeg_ao(la)
    do i = 1,levdeg_pw(lp)
      wlev(la,lp) = wlev(la,lp) + prodsq(leveigs_pw(lp,i),leveigs_ao(la,j))
    enddo
    enddo
  enddo
  enddo

! average weight of an atomic orbital state in the level of a plane-wave state

  allocate(score(neig_ao,nused))
  allocate(iassign(neig_ao))

  do ip = 1,nused
  do ia = 1,neig_ao
    la = ilev_ao(ia)
    score(ia,ip) = wlev(la,ilev_pw(ip)) / levdeg_ao(la)
  enddo
  enddo

  call num_assign_max(neig_ao, nused, score, iassign, neig_ao)

! clusters of atomic orbital and plane-wave levels that mix
! (graph nodes are the atomic orbital levels followed by the plane-wave levels)

  nlevu = nlev_pw
  if(nlev_pw > 1) nlevu = nlev_pw - 1

  allocate(iclus(nlev_ao+nlevu))

  do n = 1,nlev_ao+nlevu
    iclus(n) = n
  enddo

  do lp = 1,nlevu
  do la = 1,nlev_ao
    if(wlev(la,lp) >= XLINK*levdeg_ao(la)) then
      i = la
      do while(iclus(i) /= i)
        i = iclus(i)
      enddo
      j = nlev_ao + lp
      do while(iclus(j) /= j)
        j = iclus(j)
      enddo
      if(i /= j) iclus(max(i,j)) = min(i,j)
    endif
  enddo
  enddo

  do n = 1,nlev_ao+nlevu
    i = n
    do while(iclus(i) /= i)
      i = iclus(i)
    enddo
    iclus(n) = i
  enddo

! accepts atomic orbital levels where all states are well matched

  do la = 1,nlev_ao

    lok = .TRUE.
    eavg = ZERO
    do j = 1,levdeg_ao(la)
      ia = leveigs_ao(la,j)
      if(iassign(ia) == 0) then
        lok = .FALSE.
      elseif(iclus(nlev_ao+ilev_pw(iassign(ia))) /= iclus(la)) then
        lok = .FALSE.
      else
        eavg = eavg + ei_pw(iassign(ia))
      endif
    enddo

!   weight of the level in the plane-wave levels of its cluster.
!   It does not depend on the choice of basis inside degenerate plane-wave levels,
!   adds the parts of a plane-wave level split by numerical noise,
!   and the plane-wave levels that mix with it (spin-orbit).

    xsum = ZERO
    do lp = 1,nlevu
      if(iclus(nlev_ao+lp) == iclus(la)) xsum = xsum + wlev(la,lp)
    enddo
    xsum = xsum / levdeg_ao(la)
    if(xsum < xmin) lok = .FALSE.

    if(lok) then
      eavg = eavg / levdeg_ao(la)
      do j = 1,levdeg_ao(la)
        ia = leveigs_ao(la,j)
        ei_match(ia) = eavg
        imatch(ia) = iassign(ia)
        xover(ia) = xsum
      enddo
    else
      do j = 1,levdeg_ao(la)
        ia = leveigs_ao(la,j)
        ei_match(ia) = ei_ao(ia)
        imatch(ia) = 0
        xover(ia) = xsum
      enddo
    endif

  enddo

! warnings

  if(ipr > 0) then
    do ia = 1,neig_ao
      if(imatch(ia) > 0) then
        ip = imatch(ia)
        if(abs(ei_match(ia) - ei_ao(ia)) > UM/HARTREE .or.               &
           aonorm(ip) < XNORM) then
          write(6,'("   WARNING ao_match_state: ao state",i5," (",f10.4, &
              &     ") -> pw state",i5," (",f10.4,") eV, weight",f6.3,   &
              &     ", pw state in ao space",f6.3)')                     &
              ia, ei_ao(ia)*HARTREE, ip, ei_pw(ip)*HARTREE,              &
              xover(ia), aonorm(ip)
        endif
      endif
    enddo
  endif

! table of overlaps

  if(ipr > 1) then
    allocate(iprt(neig_pw))
    write(6,*)
    write(6,'("   ao state   e_ao      e_match   pw  weight  |",         &
        &     "  pw states with |<pw|ao>|**2 >",f5.2)') XPRINT
    do ia = 1,neig_ao
      n = 0
      do ip = 1,neig_pw
        if(prodsq(ip,ia) > XPRINT) then
          n = n + 1
          iprt(n) = ip
        endif
      enddo
      if(n > 0 .or. ia <= nused) then
        write(6,'(3x,i5,2f10.4,i5,f8.3,"  |",20(i5,":",f5.2))')          &
            ia, ei_ao(ia)*HARTREE, ei_match(ia)*HARTREE, imatch(ia),     &
            xover(ia), (iprt(i), prodsq(iprt(i),ia), i = 1,n)
      endif
    enddo
    write(6,*)
    write(6,'("   pw state   e_pw      norm in ao space")')
    do ip = 1,neig_pw
      write(6,'(3x,i5,f10.4,f10.3)') ip, ei_pw(ip)*HARTREE, aonorm(ip)
    enddo
    write(6,*)
    deallocate(iprt)
  endif

  deallocate(prod)
  deallocate(prodsq)
  deallocate(aonorm)

  deallocate(levdeg_pw, levdeg_ao)
  deallocate(leveigs_pw, leveigs_ao)
  deallocate(ilev_pw, ilev_ao)

  deallocate(wlev)
  deallocate(score)
  deallocate(iassign)
  deallocate(iclus)

  return

end subroutine ao_match_state
