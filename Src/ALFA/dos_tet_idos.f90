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

!>  calculates the contribution to the density of states
!>  and the integrated density of states from one tetrahedron
!>  implements tetra_QE.pdf

subroutine dos_tet_idos(vol,e,nhist,ehist,dhist,chist,chlow,lunif)

! written November 9 1987. jlm
! modified November 17, 2013. jlm
! Modified, documentation, 19 September 2020. JLM
! Modified, sorting, chlow, 19 October 2020. JLM
! copyright  J.L.Martins, INESC-MN.

! version 4.98 of cpw

  implicit none

  integer, parameter  :: REAL64 = selected_real_kind(12)

! input:

  real(REAL64), intent(in)        ::  vol                                !<  weighted volume of tetrahedron in k-space
  real(REAL64), intent(in)        ::  e(4)                               !<  energy levels at corners of tetrahedron
  integer, intent(in)             ::  nhist                              !<  number of points in histograms
  real(REAL64), intent(in)        ::  ehist(nhist)                       !<  energy levels in histograms
  logical, intent(in)             ::  lunif                              !<  true if grid is uniform

! input and output

  real(REAL64), intent(inout)     ::  dhist(nhist)                       !<  accumulated density of states
  real(REAL64), intent(inout)     ::  chist(nhist)                       !<  accumulated integrated density of states
  real(REAL64), intent(inout)     ::  chlow(nhist)                       !<  add to chist for n >= j at the end 

! local variables

  real(REAL64)    ::  e1, e2, e3, e4
  real(REAL64)    ::  d21, d31, d41, d32, d42, d43
  real(REAL64)    ::  dn1, dn2, dn3, dn4
  real(REAL64)    ::  ei, de, ce
  integer         ::  nlow, nhigh
  real(REAL64)    ::  dener
  real(REAL64)    ::  t

! counters

  integer    ::      i

! constants

  real(REAL64), parameter :: DELTA = 1.0E-9
  real(REAL64), parameter :: ZERO = 0.0_REAL64    

! sort energylevels at corners in increasing order

  e1 = e(1)
  e2 = e(2)
  e3 = e(3)
  e4 = e(4)

  if(e1 > e2) then
    t = e1
    e1 = e2
    e2 = t
  endif

  if(e3 > e4) then
    t = e4
    e4 = e3
    e3 = t
  endif

  if(e1 > e3) then
    t = e3
    e3 = e1
    e1 = t
  endif

  if(e2 > e4) then
    t = e4
    e4 = e2
    e2 = t
  endif

  if(e2 > e3) then
    t = e3
    e3 = e2
    e2 = t
  endif

! if energies slightly different set them equal

  if (abs(e2-e1) <= DELTA) e2 = e1
  if (abs(e3-e2) <= DELTA) e3 = e2
  if (abs(e4-e3) <= DELTA) e4 = e3

! finds range of relevant histogram points

  if(lunif) then
    if(nhist == 1) then
      nlow = 1
      nhigh = nhist
    else
      dener = (ehist(nhist) - ehist(1) ) / (nhist-1)
      nlow = floor( (e1 - ehist(1)) / dener )
      if(nlow < 1) nlow = 1
      if(nlow > nhist) nlow = nhist
      nhigh  = floor( (e4 - ehist(1) ) / dener ) + 2
      if(nhigh > nhist) nhigh = nhist
      if(nhigh < 1) nhigh = 1
    endif
  else
    nlow = 1
    nhigh = nhist
  endif

! find energy differences

  d21 = e2 - e1
  d31 = e3 - e1
  d41 = e4 - e1
  d32 = e3 - e2
  d42 = e4 - e2
  d43 = e4 - e3

! compute denominators

  dn1 = (d41*d42*d43) / vol
  dn2 = (d41*d42*d32) / vol
  dn3 = (d41*d32*d31) / vol
  dn4 = (d41*d31*d21) / vol

! loop over energylevels in histogram

  do i = nlow,nhigh

    ei = ehist(i)

    if (ei >= e1) then

      if (ei > e4) then

!       e > e4

        ce = vol
        chist(i) = chist(i) + ce

      else if (ei < e2) then

!       e1 <= e < e2

        de = (ei-e1)*(ei-e1) / dn4
        ce = (ei-e1)*de
        dhist(i) = dhist(i) + 3 * de
        chist(i) = chist(i) + ce

      else if (ei <= e3) then

!       e2 <= e <= e3

        if (e2 /= e3) then

          de = (e4-ei)*(ei-e2)/dn2 + (e3-ei)*(ei-e1)/dn3

          ce = (ei-e2)*(ei-e2) * (3*(e4-ei)+(ei-e2)) / (2*dn2)            &
             + (ei-e1)*(ei-e1) * (3*(e3-ei)+(ei-e1)) / (2*dn3)            &
             - d21*d21 * (3*d32+d21) / (2*dn3)
          if (dn4 > ZERO) ce = ce + d21*d21*d21/dn4
          dhist(i) = dhist(i) + 3*de
          chist(i) = chist(i) + ce
        else

!         e2 = e = e3

          de = ZERO

          if (e1 == e4) then
            write(6,"('  flat band at e = ',e12.4)") e1
          else
            de = vol/d41
          endif

          if (e1 == e2) de = de/2
          if (e3 == e4) de = de/2
          ce = ZERO
          if (dn4 > ZERO) ce = d21*d21*d21/dn4
          dhist(i) = dhist(i) + 3*de
          chist(i) = chist(i) + ce

        endif

      else if (ei <= e4) then

!       e3 < e <= e4

        de = (e4-ei)*(e4-ei) / dn1
        ce = vol - (e4-ei)* de
        dhist(i) = dhist(i) + 3*de
        chist(i) = chist(i) + ce

      endif

    endif                                                                !  ei > e1

  enddo                                                                  !  loop over histogram points

! points above range for uniform grids

  if(nhist > nhigh) then

    do i = nhigh+1,nhist
      ce = vol
      chist(i) = chist(i) + ce
    enddo

!   The previous loop can take 2/3 of the computing time
!    chlow(nhigh+1) = chlow(nhigh+1) + vol

  endif

  return
end subroutine dos_tet_idos
