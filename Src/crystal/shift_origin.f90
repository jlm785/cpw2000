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

!>  This program finds the closest atom to the origin and gives the
!>  shortest shift that brings the nearest image of that atom to the origin.
!>  ratio is the ratio of the distances to the origin of the nearest and
!>  second nearest atoms (1 if there is only one atom).
!>
!>  \author       José Luís Martins
!>  \version      5.13
!>  \date         10 January 2017, 10 October 2026.
!>  \copyright    GNU Public License v2

subroutine shift_origin(adot, ntype, natom, rat, shift, ratio,           &
    mxdtyp, mxdatm)

! writen January 10, 2017. JLM
! Modified, documentation, August 2019. JLM
! Indentation. 28 September 2026. JLM+claude
! near_dist_image: exact nearest image, shortest shift, ratio of distances,
! no assumption on natom(1). 10 October 2026. JLM+claude


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of types of atoms

  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i

! output

  real(REAL64), intent(out)          ::  shift(3)                        !<  shortest shift that brings the nearest atom to the origin (lattice coordinates)
  real(REAL64), intent(out)          ::  ratio                           !<  ratio of the distances of the nearest and second nearest atoms (1 if only one atom)

! local variables

  real(REAL64)       ::  dist, distmin, distmin2
  real(REAL64)       ::  r0(3), rimg(3), rmin(3)
  integer            ::  natot
  integer            ::  ntmin, jmin
  logical            ::  lfound                                          !  a second atom was found

! parameters

  real(REAL64), parameter :: ZERO = 0.0_REAL64, UM = 1.0_REAL64

! counters

  integer i, j, nt


  natot = 0
  do nt = 1,ntype
    natot = natot + natom(nt)
  enddo

  do i = 1,3
    r0(i) = ZERO
    rmin(i) = ZERO
  enddo

! nearest atom (nearest image) to the origin

  ntmin = 0
  jmin = 0
  distmin = ZERO

  do nt = 1,ntype
  do j = 1,natom(nt)

    call near_dist_image(dist, rimg, adot, r0, rat(:,j,nt))

    if(ntmin == 0 .or. dist < distmin) then
      ntmin = nt
      jmin = j
      distmin = dist
      rmin(:) = rimg(:)
    endif

  enddo
  enddo

  do i = 1,3
    shift(i) = -rmin(i)
  enddo

! second nearest atom

  if(natot < 2) then

    ratio = UM

  else

    distmin2 = ZERO
    lfound = .FALSE.

    do nt = 1,ntype
    do j = 1,natom(nt)

      if(nt /= ntmin .or. j /= jmin) then

        call near_dist(dist, adot, r0, rat(:,j,nt))

        if(.not. lfound .or. dist < distmin2) then
          lfound = .TRUE.
          distmin2 = dist
        endif

      endif

    enddo
    enddo

    if(distmin2 > ZERO) then
      ratio = distmin/distmin2
    else
      ratio = UM
    endif

  endif

  return

end subroutine shift_origin
