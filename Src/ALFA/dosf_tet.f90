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

!>  calculates the contribution to the density of states  times a function
!>  from one tetrahedron
!>  implements tetra_QE.pdf (from quantm espresso)
!>
!>  \author       Jose Luis Martins, Carlos Loia Reis
!>  \version      4.99
!>  \date         2020-2021
!>  \copyright    GNU Public License v2

subroutine dosf_tet(vol,e,Fk,nhist,ehist,dhist,lunif)

! Adapted from dos_cub July 2020. CLR
! Modified, documentation, 19 September 2020. JLM
! copyright  J.L.Martins, Carlos Loia Reis, INESC-MN.

  implicit none

  integer, parameter  :: REAL64 = selected_real_kind(12)

! input:
  real(REAL64), intent(in)        ::  vol                           !<  weighted volume of tetrahedron in k-space
  real(REAL64), intent(in)        ::  e(4)                          !<  energy levels at corners of tetrahedron

  real(REAL64), intent(in)        ::  Fk(4)                         !<  values of the function at corners of tetrahedron

  integer, intent(in)             ::  nhist                         !<  number of points in histograms
  real(REAL64), intent(in)        ::  ehist(nhist)                  !<  energy levels in histograms
  logical, intent(in)             ::  lunif                         !<  true if grid is uniform

! indut and output

  real(REAL64), intent(inout)     ::  dhist(nhist)                  !<  accumulated density of states

! local variables

  real(REAL64)    ::  e1, e2, e3, e4
  real(REAL64)    ::  d21, d31, d41, d32, d42, d43
  real(REAL64)    ::  dn1, dn2, dn3, dn4
  real(REAL64)    ::  ei
  integer         ::  nlow, nhigh
  real(REAL64)    ::  dener

  real(REAL64)    ::  d12, d13, d14, d23, d24, d34

  REAL(REAL64)    ::  f12, f13, f14
  REAL(REAL64)    ::  f21, f23, f24
  REAL(REAL64)    ::  f31, f32, f34
  REAL(REAL64)    ::  f41, f42, f43

  REAL(REAL64)    ::  fk1,fk2,fk3,fk4, g, I1, I2, I3, I4

  integer, allocatable         ::   ind(:)

! counters

  integer    ::      i

! constants

  real(REAL64), parameter :: DELTA = 1.0E-9
  real(REAL64), parameter :: ZERO = 0.0_REAL64
  real(REAL64), parameter :: UM = 1.0_REAL64

! sort energylevels at corners in increasing order

  allocate(ind(4))

  call sort(4,e,ind)

  e1 = e(ind(1))
  e2 = e(ind(2))
  e3 = e(ind(3))
  e4 = e(ind(4))

  fk1 = Fk(ind(1))
  fk2 = Fk(ind(2))
  fk3 = Fk(ind(3))
  fk4 = Fk(ind(4))

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

  d12 = d21
  d13 = d31
  d14 = d41
  d23 = d32
  d24 = d42
  d34 = d43

! loop over energylevels in histogram

  do i = nlow,nhigh

    ei = ehist(i)

       if(abs(d14) < delta ) then
         f14 = UM
       else
         f14 = (ei - e4)/(e1-e4)
         f41 = (ei - e1)/(e4-e1)
       endif

       if(abs(d13) < delta ) then
         f13 = UM
       else
         f13 = (ei - e3)/(e1-e3)
         f31 = (ei - e1)/(e3-e1)
       endif

       if(abs(d23) < delta ) then
         f23 = UM
       else
         f23 = (ei - e3)/(e2-e3)
         f32 = (ei - e2)/(e3-e2)
       endif

       if(abs(d24) < delta ) then
         f24 = UM
       else
         f24 = (ei - e4)/(e2-e4)
         f42 = (ei - e2)/(e4-e2)
       endif

       if(abs(d12) < delta ) then
         f12 = UM
       else
         f12 = (ei - e2)/(e1-e2)
         f21 = (ei - e1)/(e2-e1)
       endif

       if(abs(d34) < delta ) then
         f34 = UM
       else
         f34 = (ei - e4)/(e3-e4)
         f43 = (ei - e3)/(e4-e3)
       endif


    if (ei < e4 .and. ei >= e3) then

         g = 3*f14*f24*f34/(e4-ei)

         I1 = f14/(3*UM)
         I2 = f24/(3*UM)
         I3 = f34/(3*UM)
         I4 = (f41+f42+f43)/(3*UM)

         dhist(i) = dhist(i) + g*vol*(fk1*I1 + fk2*I2 + fk3*I3 + fk4*I4)

!!
!         de = (e4-ei)*(e4-ei) / dn1
!         dhist(i) = dhist(i) + 3*de
!!

    elseif (ei < e3 .and. ei >= e2) then


         g = 3*(f23*f31+f32*f24)/(e4-e1)

         I1 = f14/(3*UM)
         I1 = I1 + f13*f31*f23/(g*(e4-e1))

         I2 = f23/(3*UM)
         I2 = I2 + f24*f24*f32/(g*(e4-e1))

         I3 = f32/(3*UM)
         I3 = I3 + f31*f31*f23/(g*(e4-e1))

         I4 = f41/(3*UM)
         I4 = I4 + f42*f24*f32/(g*(e4-e1))


         dhist(i) = dhist(i) + g*vol*(fk1*I1 + fk2*I2 + fk3*I3 + fk4*I4)

!!
!           de = (e4-ei)*(ei-e2)/dn2 + (e3-ei)*(ei-e1)/dn3
!           dhist(i) = dhist(i) + 3*de
!!

    elseif (ei < e2 .and. ei > e1) then


        g= 3*f21*f31*f41/(ei-e1)

        I1 = (f12 + f13+ f14)/(3*UM)
        I2 = f21/(3*UM)
        I3 = f31/(3*UM)
        I4 = f41/(3*UM)


        dhist(i) = dhist(i) + g*vol*(fk1*I1 + fk2*I2 + fk3*I3 + fk4*I4)

!!
!        de = (ei-e1)*(ei-e1) / dn4
!        dhist(i) = dhist(i) + 3 * de
!!


    endif

  enddo                                                             !  loop over histogram points


  return
  end subroutine dosf_tet
