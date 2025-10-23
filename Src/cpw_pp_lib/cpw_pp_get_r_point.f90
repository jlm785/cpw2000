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

!>  Asks the user to select a point or direction .
!>  User has a choice of coordinate system.  Output is in lattice
!>  and conventional cartesian coordinates.
!>
!>  \author       Jose Luis Martins
!>  \version      5.12
!>  \date         23 October 2025.
!>  \copyright    GNU Public License v2

subroutine cpw_pp_get_r_point(rpoint, rcar, adot, typeofr, ioreplay)

! Adapted from cpw_pp_get_k_vector

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input


  integer, intent(in)                ::  ioreplay                        !<  tape number for reproducing calculations

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in real space
  character(len=20), intent(in)      ::  typeofr                         !<  label for type of r-point

! output

  real(REAL64), intent(out)          ::  rpoint(3)                       !<  choice of r-point (lattice coordinates)
  real(REAL64), intent(out)          ::  rcar(3)                         !<  choice of r-point (cartesian coordinates)

! local variables

  integer           ::  icoor               !  coordinate system

  real(REAL64)      ::  avec(3,3)           !  primitive lattice vectors
  real(REAL64)      ::  bvec(3,3)           !  reciprocal primitive lattice vectors
  real(REAL64)      ::  aconv(3,3)          !  conventional lattice vectors of the Niggli cell
  real(REAL64)      ::  avecnig(3,3)        !  primitive lattice vectors of the Niggli cell

  real(REAL64)      ::  rin(3)

! parameters

  real(REAL64), parameter       ::  ZERO = 0.0_REAL64
  real(REAL64), parameter       ::  PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter       ::  BOHR = 0.5291772109_REAL64

! counters

  integer    ::  k, j, m


! gets the k-point, but first generates coordinate system

  call adot_to_avec_aconv(adot,avec,bvec,aconv,avecnig)

  write(6,*)
  write(6,*) '   Which coordinate system do you want to use?'
  write(6,*) '   1) Primitive lattice coordinates.'
  write(6,*) '   2) Cartesian coordinates.'
  write(6,*) '   3) Conventional lattice coordinates.'
  write(6,*)
  write(6,*) '   Enter your choice (1-3).'
  write(6,*)

  read(5,*) icoor
  write(ioreplay,'(2x,i8,"   coordinate choice for ",a20)') icoor, typeofr

  if(icoor <1 .or. icoor > 3) then
    write(6,*) '   Wrong choice, using primitive lattice'
    icoor = 1
  endif

! prints information about the coordinate system

  if(icoor == 1) then

    write(6,*)
    write(6,*) '   using PRIMITIVE LATTICE COORDINATES '
    write(6,*)
    write(6,*) '   The primitive lattice vectors are: '
    write(6,*)
    write(6,'("  a1 = (",f8.4,")   a2 = (",f8.4,")   a3 = (",f8.4,")")')   &
        avec(1,1),avec(1,2),avec(1,3)
     write(6,'("       (",f8.4,")        (",f8.4,")        (",f8.4,")")')  &
        avec(2,1),avec(2,2),avec(2,3)
     write(6,'("       (",f8.4,")        (",f8.4,")        (",f8.4,")")')  &
        avec(3,1),avec(3,2),avec(3,3)
     write(6,*)

  else

    if(icoor == 2) then

      write(6,*)
      write(6,*) '   using CARTESIAN COORDINATES (atomic units)'
      write(6,*)

    else

      write(6,*)
      write(6,*) '   using CONVENTIONAL LATTICE COORDINATES '
      write(6,*)

    endif

    write(6,*) '   The conventional lattice vectors are: '
    write(6,*)
    write(6,'("  a1 = (",f8.4,")   a2 = (",f8.4,")   a3 = (",f8.4,")")')   &
         aconv(1,1),aconv(1,2),aconv(1,3)
     write(6,'("       (",f8.4,")        (",f8.4,")        (",f8.4,")")')  &
         aconv(2,1),aconv(2,2),aconv(2,3)
     write(6,'("       (",f8.4,")        (",f8.4,")        (",f8.4,")")')  &
         aconv(3,1),aconv(3,2),aconv(3,3)
     write(6,*)

  endif

  write(6,*)
  write(6,'("  Enter ",a20)') typeofr
  write(6,*)

  read(5,*) rin(1),rin(2),rin(3)
  write(ioreplay,'(5x,3f20.10,5x,"k-point")') rin(1),rin(2),rin(3)

! transforms to primitive coordinates

  if(icoor == 1) then

    do j = 1,3
      rpoint(j) = rin(j)
    enddo

    do j = 1,3
      rcar(j) = ZERO
      do k = 1,3
        rcar(j) = rcar(j) + avec(j,k)*rpoint(k)
      enddo
    enddo

  elseif(icoor == 2) then

    do j = 1,3
      rcar(j) = rin(j)
    enddo

    do j = 1,3
      rpoint(j) = ZERO
      do k = 1,3
        rpoint(j) = rpoint(j) + rcar(k)*bvec(k,j)
      enddo
      rpoint(j) = rpoint(j) / (2*PI)
    enddo

  elseif(icoor == 3) then

    do j = 1,3
      rcar(j) = ZERO
      do k = 1,3
        rcar(j) = rcar(j) + aconv(j,k)*rin(k)
      enddo
    enddo

    do j = 1,3
      rpoint(j) = ZERO
      do k = 1,3
        rpoint(j) = rpoint(j) + rcar(k)*bvec(k,j)
      enddo
      rpoint(j) = rpoint(j) / (2*PI)
    enddo

  endif


  write(6,*)
  write(6,*) '    Coordinates of the chosen ', typeofr
  write(6,'(7x,"lattice coord.",17x," cartesian coord. (bohr)",8x," cartesian coord. (Ang)")')
  write(6,*)
  write(6,'(4x,3f9.4,5x,3f9.4,5x,3f9.4)') (rpoint(j),j=1,3), (rcar(j),j=1,3), (rcar(j)*BOHR,j=1,3)
  write(6,*)

  return

end subroutine cpw_pp_get_r_point
