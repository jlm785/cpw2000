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

!>   Step subroutine to calculate the equation of state
!>
!>  \author       Jose Luis Martins
!>  \version      1.6.2 of md
!>  \date         18 March 2026.
!>  \copyright    GNU Public License v2

subroutine  move_eos(energy, adot, lepi, lfinisheos)

! Written 18 March 2026. JLM

  implicit none

  integer, parameter :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)           ::  energy                          !<  energy in hartree
  logical, intent(in)                ::  lepi                            !<  if true, only adot(3,3) is scaled

! input and output

  real(REAL64), intent(inout)        ::  adot(3,3)                       !<  metric in real space

! output

  logical, intent(out)               ::  lfinisheos                      !<  if true, the EOS has been determined

! local saved variables

  integer, save                      ::  nstatus = 1                     !  1 first call, 2 second call
  real(REAL64), save                 ::  adot_init(3,3)                  !  initial adot

  real(REAL64), save                 ::  e1, e2, e3                      !  energies for different steps
  real(REAL64), save                 ::  a1, a2, a3                      !  lattice constant "c" for different steps
  real(REAL64), save                 ::  apred, epred                    !  predicted lattice constant and energy

  real(REAL64), save                 ::  as(9), es(9)                    !  scan of lattice constant and energy
  real(REAL64), save                 ::  vs(9)                           !  scan of volume
  integer, save                      ::  ns = 0

! local variables

  real(REAL64)          ::  atmp, etmp
  real(REAL64)          ::  fac
  real(REAL64)          ::  vcell, bdot(3,3)
  real(REAL64)          ::  volfac

! parameters

  real(REAL64), parameter :: EPS = 0.00001_REAL64
  real(REAL64), parameter :: DELTA = 0.1_REAL64
  real(REAL64), parameter :: UM = 1.0_REAL64

! counters

  integer     ::  i, j


  lfinisheos = .FALSE.

  if(nstatus == 1) then

!   initial point

    e1 = energy
    a1 = sqrt(adot(3,3))
    do j = 1,3
    do i = 1,3
      adot_init(i,j) = adot(i,j)
    enddo
    enddo

    if(lepi) then
      adot(3,3) = 1.1*adot(3,3)
    else
      do j = 1,3
      do i = 1,3
        adot(i,j) = 1.1*adot(i,j)
      enddo
      enddo
    endif

    nstatus = 2

  elseif(nstatus == 2) then

!   second point

    e2 = energy
    a2 = sqrt(adot(3,3))

    if(e1 > e2) then
      if(lepi) then
        adot(3,3) = 1.1*adot(3,3)
      else
        do j = 1,3
        do i = 1,3
          adot(i,j) = 1.1*adot(i,j)
        enddo
        enddo
      endif
    else
      if(lepi) then
        adot(3,3) = 0.9*adot_init(3,3)
      else
        do j = 1,3
        do i = 1,3
          adot(i,j) = 0.9*adot_init(i,j)
        enddo
        enddo
      endif
    endif

    nstatus = 3

  elseif(nstatus == 3) then

!    brackets minimum

    e3 = energy
    a3 = sqrt(adot(3,3))

    if(e1 > e2 .and. e2 > e3) then

      e1 = e2
      a1 = a2
      e2 = e3
      a2 = a3

      if(lepi) then
        adot(3,3) = 1.1*adot(3,3)
      else
        do j = 1,3
        do i = 1,3
          adot(i,j) = 1.1*adot(i,j)
        enddo
        enddo
      endif

    elseif(e3 < e1 .and. e1 < e2) then

      e2 = e1
      a2 = a1
      e1 = e3
      a1 = a3
      if(lepi) then
        adot(3,3) = 0.9*adot(3,3)
      else
        do j = 1,3
        do i = 1,3
          adot(i,j) = 0.9*adot(i,j)
        enddo
        enddo
      endif

    else

      nstatus = 4

      apred = -(a3*a3*(e1 - e2) + a1*a1*(e2 - e3) + a2*a2*(e3 - e1)) /    &
                   (2*(a3*(e2 - e1) + a1*(e3 - e2) + a2*(e1 - e3)))
      epred =  e1*(apred - a2)*(apred - a3) / ((a1 - a2)*(a1 - a3)) +     &
               e2*(apred - a3)*(apred - a1) / ((a2 - a3)*(a2 - a1)) +     &
               e3*(apred - a1)*(apred - a2) / ((a3 - a1)*(a3 - a2))

!     reorders

      if(e1 < e2) then
        atmp = a1
        etmp = e1
        a1 = a3
        e1 = e3
        a3 = atmp
        e3 = etmp
      else
        atmp = a2
        etmp = e2
        a2 = a3
        e2 = e3
        a3 = atmp
        e3 = etmp
      endif

!     now   a1 < a3 < a2 and e3 < e1 and e3 < e2

      if(apred < a1-EPS .or. apred > a2+EPS) then

        write(6,*)
        write(6,*) '   STOPPED in move_eos:  prediction outside interval'
        write(6,'("  a1, a2 , a3, apred = ",4i5 )') a1, a2 , a3, apred
        write(6,*)

        stop

      endif

      if(apred < a3) then

        a2 = a3
        e2 = e3
        a3 = apred

      else

        a1 = a3
        e1 = e3
        a3 = apred

      endif

      fac = a3*a3 / adot(3,3)
      if(lepi) then
        adot(3,3) = a3*a3
      else
        do j = 1,3
        do i = 1,3
          adot(i,j) = fac*adot(i,j)
        enddo
        enddo
      endif

    endif

  elseif(nstatus == 4) then

!   improves minimum

    e3 = energy
    a3 = sqrt(adot(3,3))

    apred = -(a3*a3*(e1 - e2) + a1*a1*(e2 - e3) + a2*a2*(e3 - e1)) /    &
                 (2*(a3*(e2 - e1) + a1*(e3 - e2) + a2*(e1 - e3)))
    epred =  e1*(apred - a2)*(apred - a3) / ((a1 - a2)*(a1 - a3)) +     &
             e2*(apred - a3)*(apred - a1) / ((a2 - a3)*(a2 - a1)) +     &
             e3*(apred - a1)*(apred - a2) / ((a3 - a1)*(a3 - a2))

    if(apred < a3) then

      a2 = a3
      e2 = e3
      a3 = apred

    else

      a1 = a3
      e1 = e3
      a3 = apred

    endif

    fac = a3*a3 / adot(3,3)
    if(lepi) then
      adot(3,3) = a3*a3
    else
      do j = 1,3
      do i = 1,3
        adot(i,j) = fac*adot(i,j)
      enddo
      enddo
    endif

    if((a2 - a1) < DELTA/10) then

!     good enough estimate of the lattice constant

      write(6,*)
      write(6,'(5x,f14.3,"   move_eos lattice constant estimate")') a3
      write(6,*)
      nstatus = 5
      do i = 1,9
        as(i) = a3 + (i-5)*DELTA
      enddo
      a3 = as(1)
      fac = a3*a3 / adot(3,3)
      if(lepi) then
        adot(3,3) = a3*a3
      else
        do j = 1,3
        do i = 1,3
          adot(i,j) = fac*adot(i,j)
        enddo
        enddo
      endif
      call adot_to_bdot(adot,vcell,bdot)
      vs(1) = vcell

    endif

  else

!   writes eqst.dat

    ns = ns+1
    es(ns) = energy

    if(ns < 9) then

      fac = as(ns+1)*as(ns+1) / adot(3,3)
      if(lepi) then
        adot(3,3) = as(ns+1)*as(ns+1)
      else
        do j = 1,3
        do i = 1,3
          adot(i,j) = fac*adot(i,j)
        enddo
        enddo
      endif
      call adot_to_bdot(adot,vcell,bdot)
      vs(ns+1) = vcell

    else

      open(unit = 12, file = 'eqst.dat', form = 'formatted')

      write(12,'("   BIRCH")')
      write(12,'(i5)') 1
      write(12,'(i5)') 9
      if(lepi) then
        write(12,'(f14.8,"     move\_eos")') -1.0
        do i = 1,9
          write(12,'(f18.8,f14.8)') vs(i), es(i)
        enddo
      else
        a3 = sqrt(adot(3,3))
        fac = vs(9) / (a3*a3*a3)
        write(12,'(f14.8,"     move\_eos")') fac
        do i = 1,9
          write(12,'(2f14.8)') as(i),es(i)
        enddo
      endif
      write(12,'(f8.3,i5)') 4.0, 101

      close(unit = 12)

      fac = as(5)*as(5) / adot(3,3)
      if(lepi) then
         adot(3,3) = as(5)*as(5)
      else
        do j = 1,3
        do i = 1,3
          adot(i,j) = fac*adot(i,j)
        enddo
        enddo
      endif

      lfinisheos = .TRUE.

    endif

  endif

  return

end subroutine  move_eos
