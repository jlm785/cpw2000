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

!>  Asks the user to select a k-vector or direction.
!>  User has a choice of coordinate system.  Output is in lattice
!>  and conventional cartesian coordinates.
!>
!>  \author       Jose Luis Martins
!>  \version      5.12
!>  \date         6 November 2023. 24 September 2025.
!>  \copyright    GNU Public License v2

subroutine cpw_pp_get_k_vector(rkpt, rkcar, adot, typeofk, ioreplay)

! extracted from old out_effective_mass. 22 September 2023.
! Writes the coordinates. 24 September 2025. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input


  integer, intent(in)                ::  ioreplay                        !<  tape number for reproducing calculations

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in real space
  character(len=20), intent(in)      ::  typeofk                         !<  label for type of k-point

! output

  real(REAL64), intent(out)          ::  rkpt(3)                         !<  choice of k-vector (lattice coordinates)
  real(REAL64), intent(out)          ::  rkcar(3)                        !<  choice of k-vector (cartesian coordinates)

! local variables

  integer           ::  icoor               !  coordinate system

  real(REAL64)      ::  avec(3,3)           !  primitive lattice vectors
  real(REAL64)      ::  bvec(3,3)           !  reciprocal primitive lattice vectors
  real(REAL64)      ::  aconv(3,3)          !  conventional lattice vectors of the Niggli cell
  real(REAL64)      ::  avecnig(3,3)        !  primitive lattice vectors of the Niggli cell
  real(REAL64)      ::  bdot(3,3), vcell

  real(REAL64)      ::  adotconv(3,3)
  real(REAL64)      ::  bvecconv(3,3), tmp(3,3)

  real(REAL64)      ::  rkin(3)

! parameters

  real(REAL64), parameter       ::  ZERO = 0.0_REAL64
  real(REAL64), parameter       ::  PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter       ::  BOHR = 0.5291772109_REAL64

! counters

  integer    ::  k, j, m


! gets the k-point, but first generates coordinate system

  call adot_to_avec_aconv(adot,avec,bvec,aconv,avecnig)
  call adot_to_bdot(adot,vcell,bdot)

  write(6,*)
  write(6,*) '   Which coordinate system do you want to use?'
  write(6,*) '   1) Primitive lattice coordinates.'
  write(6,*) '   2) Cartesian coordinates.'
  write(6,*) '   3) Conventional lattice coordinates.'
  write(6,*)
  write(6,*) '   Enter your choice (1-3).'
  write(6,*)

  read(5,*) icoor
  write(ioreplay,'(2x,i8,"   coordinate choice for ",a20)') icoor, typeofk

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
  write(6,'("  Enter ",a20)') typeofk
  write(6,*)

  read(5,*) rkin(1),rkin(2),rkin(3)
  write(ioreplay,'(5x,3f20.10,5x,"k-point")') rkin(1),rkin(2),rkin(3)

! transforms to primitive coordinates

  if(icoor == 1) then

    do j = 1,3
      rkpt(j) = rkin(j)
    enddo

    do j = 1,3
      rkcar(j) = ZERO
      do k = 1,3
        rkcar(j) = rkcar(j) + bvec(j,k)*rkpt(k)
      enddo
    enddo

  elseif(icoor == 2) then

    do j = 1,3
      rkpt(j) = ZERO
      do k = 1,3
        rkpt(j) = rkpt(j) + rkin(k)*avec(k,j)
      enddo
      rkpt(j) = rkpt(j) / (2*PI)
    enddo

    do j = 1,3
      rkcar(j) = rkin(j)
    enddo

  elseif(icoor == 3) then

    do j = 1,3
    do k = 1,3
      adotconv(k,j) = ZERO
      do m = 1,3
        adotconv(k,j) = adotconv(k,j) + aconv(m,j)*aconv(m,k)
      enddo
    enddo
    enddo

    call adot_to_avec(adotconv,tmp,bvecconv)

    do j = 1,3
      rkcar(j) = ZERO
      do k = 1,3
        rkcar(j) = rkcar(j) + bvecconv(j,k)*rkin(k)
      enddo
    enddo

    do j = 1,3
      rkpt(j) = ZERO
      do k = 1,3
        rkpt(j) = rkpt(j) + rkcar(k)*avec(k,j)
      enddo
      rkpt(j) = rkpt(j) / (2*PI)
    enddo

  endif


  write(6,*)
  write(6,*) '    Coordinates of the chosen ', typeofk
  write(6,'(7x,"lattice coord.",17x," cartesian coord. (bohr-1)",6x," cartesian coord. (Ang-1)")')
  write(6,*)
  write(6,'(4x,3f9.4,5x,3f9.4,5x,3f9.4)') (rkpt(j),j=1,3), (rkcar(j),j=1,3), (rkcar(j)/BOHR,j=1,3)
  write(6,*)

  return

end subroutine cpw_pp_get_k_vector
