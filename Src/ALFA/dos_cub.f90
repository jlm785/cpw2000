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

!>  computes the density of states from one cube with corners ec

subroutine dos_cub(vol,ec,nhist,ehist,dhist,chist,chlow,lunif,lidos)

! written november 9 1987. jlm
! modified November 17, 2013. jlm
! Modified, documentation, 19 September 2020. JLM
! Modified efficient idos, 19 October 2020. JLM
! copyright  J.L.Martins, INESC-MN.

! version 4.98 of cpw

  implicit none

  integer, parameter  :: REAL64 = selected_real_kind(12)

! input:

  real(REAL64), intent(in)        ::  vol                                !<  weighted volume of terahedron in k-space
  real(REAL64), intent(in)        ::  ec(8)                              !<  energy levels at corners of cube
! bottom  4 3   top  8 7
!         1 2        5 6
  integer, intent(in)             ::  nhist                              !<  number of points in histograms
  real(REAL64), intent(in)        ::  ehist(nhist)                       !<  energy levels in histograms
  logical, intent(in)             ::  lunif                              !<  true if grid is uniform
  logical, intent(in)             ::  lidos                              !<  true if integrated density of states is to be computed.

! indut and output

  real(REAL64), intent(inout)     ::  dhist(nhist)                       !<  accumulated density of states
  real(REAL64), intent(inout)     ::  chist(nhist)                       !<  accumulated integrated density of states
  real(REAL64), intent(inout)     ::  chlow(nhist)                       !<  add to chist for n >= j at the end 

! local variables

  real(REAL64)     ::    vol1, vol2, e(4)

  vol1 = vol/12
  vol2 = vol/6

! corner tetrahedrons (from corners 1,3,6,8)

  if(lidos) then

    e(1) = ec(1)
    e(2) = ec(2)
    e(3) = ec(4)
    e(4) = ec(5)
    call dos_tet_idos(vol1,e,nhist,ehist,dhist,chist,chlow,lunif)
    e(1) = ec(3)
    e(2) = ec(2)
    e(3) = ec(4)
    e(4) = ec(7)
    call dos_tet_idos(vol1,e,nhist,ehist,dhist,chist,chlow,lunif)
    e(1) = ec(6)
    e(2) = ec(2)
    e(3) = ec(5)
    e(4) = ec(7)
    call dos_tet_idos(vol1,e,nhist,ehist,dhist,chist,chlow,lunif)
    e(1) = ec(8)
    e(2) = ec(4)
    e(3) = ec(5)
    e(4) = ec(7)
    call dos_tet_idos(vol1,e,nhist,ehist,dhist,chist,chlow,lunif)

!   center tetrahedron

    e(1) = ec(2)
    e(2) = ec(4)
    e(3) = ec(5)
    e(4) = ec(7)
    call dos_tet_idos(vol2,e,nhist,ehist,dhist,chist,chlow,lunif)

!   do it all again

!   corner tetrahedrons (from corners 2,4,5,7)

    e(1) = ec(2)
    e(2) = ec(1)
    e(3) = ec(3)
    e(4) = ec(6)
    call dos_tet_idos(vol1,e,nhist,ehist,dhist,chist,chlow,lunif)
    e(1) = ec(4)
    e(2) = ec(1)
    e(3) = ec(3)
    e(4) = ec(8)
    call dos_tet_idos(vol1,e,nhist,ehist,dhist,chist,chlow,lunif)
    e(1) = ec(5)
    e(2) = ec(1)
    e(3) = ec(6)
    e(4) = ec(8)
    call dos_tet_idos(vol1,e,nhist,ehist,dhist,chist,chlow,lunif)
    e(1) = ec(7)
    e(2) = ec(3)
    e(3) = ec(6)
    e(4) = ec(8)
    call dos_tet_idos(vol1,e,nhist,ehist,dhist,chist,chlow,lunif)

!   center tetrahedron

    e(1) = ec(1)
    e(2) = ec(3)
    e(3) = ec(6)
    e(4) = ec(8)
    call dos_tet_idos(vol2,e,nhist,ehist,dhist,chist,chlow,lunif)

  else

    e(1) = ec(1)
    e(2) = ec(2)
    e(3) = ec(4)
    e(4) = ec(5)
    call dos_tet(vol1,e,nhist,ehist,dhist,lunif)
    e(1) = ec(3)
    e(2) = ec(2)
    e(3) = ec(4)
    e(4) = ec(7)
    call dos_tet(vol1,e,nhist,ehist,dhist,lunif)
    e(1) = ec(6)
    e(2) = ec(2)
    e(3) = ec(5)
    e(4) = ec(7)
    call dos_tet(vol1,e,nhist,ehist,dhist,lunif)
    e(1) = ec(8)
    e(2) = ec(4)
    e(3) = ec(5)
    e(4) = ec(7)
    call dos_tet(vol1,e,nhist,ehist,dhist,lunif)

!  center tetrahedron

    e(1) = ec(2)
    e(2) = ec(4)
    e(3) = ec(5)
    e(4) = ec(7)
    call dos_tet(vol2,e,nhist,ehist,dhist,lunif)

!   do it all again

!   corner tetrahedrons (from corners 2,4,5,7)

    e(1) = ec(2)
    e(2) = ec(1)
    e(3) = ec(3)
    e(4) = ec(6)
    call dos_tet(vol1,e,nhist,ehist,dhist,lunif)
    e(1) = ec(4)
    e(2) = ec(1)
    e(3) = ec(3)
    e(4) = ec(8)
    call dos_tet(vol1,e,nhist,ehist,dhist,lunif)
    e(1) = ec(5)
    e(2) = ec(1)
    e(3) = ec(6)
    e(4) = ec(8)
    call dos_tet(vol1,e,nhist,ehist,dhist,lunif)
    e(1) = ec(7)
    e(2) = ec(3)
    e(3) = ec(6)
    e(4) = ec(8)
    call dos_tet(vol1,e,nhist,ehist,dhist,lunif)

!   center tetrahedron

    e(1) = ec(1)
    e(2) = ec(3)
    e(3) = ec(6)
    e(4) = ec(8)
    call dos_tet(vol2,e,nhist,ehist,dhist,lunif)

  endif

  return
end subroutine dos_cub
