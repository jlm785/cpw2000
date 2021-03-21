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

!>  computes the density of states times a function from one cube with corners ec
!>
!>  \author       Jose Luis Martins, Carlos Loia Reis
!>  \version      4.99
!>  \date         2020-2021
!>  \copyright    GNU Public License v2

subroutine dosf_cub(vol,ec,Fkc, nhist,ehist,dhist,lunif)

! Adapted from dos_cub July 2020. CLR
! Modified, documentation, 19 September 2020. JLM
! copyright  J.L.Martins, Carlos Loia Reis, INESC-MN.

  implicit none

  integer, parameter  :: REAL64 = selected_real_kind(12)

! input:

  real(REAL64), intent(in)        ::  vol                                !<  weighted volume of terahedron in k-space
  real(REAL64), intent(in)        ::  ec(8)                              !<  energy levels at corners of cube
! bottom  4 3   top  8 7
!         1 2        5 6

  real(REAL64), intent(in)        ::  Fkc(8)                             !<  function at corners of cube

  integer, intent(in)             ::  nhist                              !<  number of points in histograms
  real(REAL64), intent(in)        ::  ehist(nhist)                       !<  energy levels in histograms
  logical, intent(in)             ::  lunif                              !<  true if grid is uniform

! input and output

  real(REAL64), intent(inout)     ::  dhist(nhist)                       !<  accumulated density of states

! local variables

  real(REAL64)                    ::    vol1, vol2
  real(REAL64), allocatable       ::    e(:), Fk(:)


  allocate(e(4), Fk(4))

  vol1 = vol/12
  vol2 = vol/6

! corner tetrahedrons (from corners 1,3,6,8)

  e(1) = ec(1)
  e(2) = ec(2)
  e(3) = ec(4)
  e(4) = ec(5)

  Fk(1) = Fkc(1)
  Fk(2) = Fkc(2)
  Fk(3) = Fkc(4)
  Fk(4) = Fkc(5)

  call dosf_tet(vol1,e,Fk,nhist,ehist,dhist,lunif)

  e(1) = ec(3)
  e(2) = ec(2)
  e(3) = ec(4)
  e(4) = ec(7)

  Fk(1) = Fkc(3)
  Fk(2) = Fkc(2)
  Fk(3) = Fkc(4)
  Fk(4) = Fkc(7)

  call dosf_tet(vol1,e,Fk,nhist,ehist,dhist,lunif)

  e(1) = ec(6)
  e(2) = ec(2)
  e(3) = ec(5)
  e(4) = ec(7)

  Fk(1) = Fkc(6)
  Fk(2) = Fkc(2)
  Fk(3) = Fkc(5)
  Fk(4) = Fkc(7)

  call dosf_tet(vol1,e,Fk,nhist,ehist,dhist,lunif)

  e(1) = ec(8)
  e(2) = ec(4)
  e(3) = ec(5)
  e(4) = ec(7)

  Fk(1) = Fkc(8)
  Fk(2) = Fkc(4)
  Fk(3) = Fkc(5)
  Fk(4) = Fkc(7)

  call dosf_tet(vol1,e,Fk,nhist,ehist,dhist,lunif)

! center tetrahedron

  e(1) = ec(2)
  e(2) = ec(4)
  e(3) = ec(5)
  e(4) = ec(7)

  Fk(1) = Fkc(2)
  Fk(2) = Fkc(4)
  Fk(3) = Fkc(5)
  Fk(4) = Fkc(7)

  call dosf_tet(vol2,e,Fk,nhist,ehist,dhist,lunif)

! do it all again

! corner tetrahedrons (from corners 2,4,5,7)

  e(1) = ec(2)
  e(2) = ec(1)
  e(3) = ec(3)
  e(4) = ec(6)

  Fk(1) = Fkc(2)
  Fk(2) = Fkc(1)
  Fk(3) = Fkc(3)
  Fk(4) = Fkc(6)

  call dosf_tet(vol1,e,Fk,nhist,ehist,dhist,lunif)

  e(1) = ec(4)
  e(2) = ec(1)
  e(3) = ec(3)
  e(4) = ec(8)

  Fk(1) = Fkc(4)
  Fk(2) = Fkc(1)
  Fk(3) = Fkc(3)
  Fk(4) = Fkc(8)

  call dosf_tet(vol1,e,Fk,nhist,ehist,dhist,lunif)

  e(1) = ec(5)
  e(2) = ec(1)
  e(3) = ec(6)
  e(4) = ec(8)

  Fk(1) = Fkc(5)
  Fk(2) = Fkc(1)
  Fk(3) = Fkc(6)
  Fk(4) = Fkc(8)

  call dosf_tet(vol1,e,Fk,nhist,ehist,dhist,lunif)

  e(1) = ec(7)
  e(2) = ec(3)
  e(3) = ec(6)
  e(4) = ec(8)

  Fk(1) = Fkc(7)
  Fk(2) = Fkc(3)
  Fk(3) = Fkc(6)
  Fk(4) = Fkc(8)

  call dosf_tet(vol1,e,Fk,nhist,ehist,dhist,lunif)

! center tetrahedron

  e(1) = ec(1)
  e(2) = ec(3)
  e(3) = ec(6)
  e(4) = ec(8)

  Fk(1) = Fkc(1)
  Fk(2) = Fkc(3)
  Fk(3) = Fkc(6)
  Fk(4) = Fkc(8)

  call dosf_tet(vol2,e,Fk,nhist,ehist,dhist,lunif)


  return
end subroutine dosf_cub
