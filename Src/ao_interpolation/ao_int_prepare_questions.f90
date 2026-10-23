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

!>  Asks the choices for the calculation of the modified tight-binding
!>  (MTB) hamiltonian:  grid of k-points, spin-orbit, optical matrix
!>  elements and how the plane-wave energies are transferred.
!>
!>  \author       claude
!>  \version      5.13
!>  \date         23 October 2026.
!>  \copyright    GNU Public License v2

subroutine ao_int_prepare_questions(ioreplay, nk1, nk2, nk3,             &
    lso, loptical, imtb)

! Written 23 October 2026 from ao_interpolation_prepare. claude

  implicit none

! input

  integer, intent(in)                ::  ioreplay                        !<  tape number for reproducing calculations

! output

  integer, intent(out)               ::  nk1, nk2, nk3                   !<  grid of k-points
  logical, intent(out)               ::  lso                             !<  if true uses spin-orbit
  logical, intent(out)               ::  loptical                        !<  if true calculates optical matrix elements
  integer, intent(out)               ::  imtb                            !<  MTB method, 1 projection of plane-wave states, 2 matching

! local variables

  character(len=1)   ::  yesno
  integer            ::  ios


  write(6,*)
  write(6,*) '  please enter grid of k-points'
  write(6,*) '  recomended is 4 4 4 minimum'

  read(5,*) nk1,nk2,nk3
  write(ioreplay,*) nk1,nk2,nk3,'     k-point grid'

  write(6,'( "  interpolation in a mesh of", i3," x ",i3," x ",i3,       &
      &      " k-points"  )') nk1,nk2,nk3

  write(6,*)
  write(6,*) '  do you want a calculation including Spin-Orbit (y/n) ?'
  write(6,*)
  read(5,*) yesno
  write(ioreplay,*) yesno,'     spin-orbit'

  if(yesno == 'y' .or. yesno == 'Y') then
    lso = .TRUE.
  else
    lso = .FALSE.
  endif

  if(lso) then
    write(6,*)
    write(6,*) '  calculation will be done using SO perturbation.'
    write(6,*)
  endif

  write(6,*)
  write(6,*) '  do you want a calculation including Hamiltonian Derivatives (y/n) ?'
  write(6,*) '  this will be needed to compute the dielectric function '
  write(6,*) '  using interpolation (it is memory intensive...).'
  write(6,*)
  read(5,*) yesno
  write(ioreplay,*) yesno,'     optical'

  if(yesno == 'y' .or. yesno == 'Y') then
    loptical = .TRUE.
  else
    loptical = .FALSE.
  endif
  if(loptical) then
    write(6,*)
    write(6,*) '  calculation will also include dh0drk computation.'
    write(6,*)
  endif

  write(6,*)
  write(6,*) '  How should the plane-wave energies be transferred to the'
  write(6,*) '  atomic orbital hamiltonian?'
  write(6,*) '  1)  projection of the plane-wave states (recommended)'
  write(6,*) '  2)  matching of atomic orbital and plane-wave states (older method)'
  write(6,*)
  read(5,*,iostat=ios) imtb
  if(ios /= 0) imtb = 1
  if(imtb /= 2) imtb = 1
  write(ioreplay,*) imtb,'     MTB method'

  if(imtb == 1) then
    write(6,*)
    write(6,*) '  plane-wave states projected on the atomic orbitals'
    write(6,*)
  else
    write(6,*)
    write(6,*) '  atomic orbital and plane-wave states are matched'
    write(6,*)
  endif

  return

end subroutine ao_int_prepare_questions
