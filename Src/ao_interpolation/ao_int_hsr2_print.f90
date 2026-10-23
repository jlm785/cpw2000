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

!>  Prints a summary of the setup of the modified tight-binding
!>  data in hsr2.dat.
!>
!>  \author       claude
!>  \version      5.13
!>  \date         23 October 2026.
!>  \copyright    GNU Public License v2

subroutine ao_int_hsr2_print(title, subtitle, ztot, emax,                &
    ntype, natom, nameat, norbat, lorb,                                  &
    lso, loptical, imtb, mp_grid, nband,                                 &
    mxdtyp, mxdlao)

! Written 23 October 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdlao                          !<  array dimension of orbital per atom type

  character(len=50), intent(in)      ::  title                           !<  title for plots
  character(len=140), intent(in)     ::  subtitle                        !<  subtitle for plots
  real(REAL64), intent(in)           ::  ztot                            !<  total charge density (electrons/cell)
  real(REAL64), intent(in)           ::  emax                            !<  kinetic energy cutoff of plane wave expansion (hartree)
  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  character(len=2), intent(in)       ::  nameat(mxdtyp)                  !<  chemical symbol for the type i
  integer, intent(in)                ::  norbat(1,mxdtyp)                !<  number of atomic orbitals of atom type nt
  integer, intent(in)                ::  lorb(mxdlao,1,mxdtyp)           !<  angular momentum of orbital n of atom type nt
  integer, intent(in)                ::  lso                             !<  1 with spin-orbit, 0 without
  integer, intent(in)                ::  loptical                        !<  1 with optical matrix elements, 0 without
  integer, intent(in)                ::  imtb                            !<  MTB method, 1 projection, 2 matching
  integer, intent(in)                ::  mp_grid(3)                      !<  grid of k-points of the MTB calculation
  integer, intent(in)                ::  nband                           !<  number of basis functions (spinors with spin-orbit)

! local variables

  character(len=1)   ::  lname(0:3)                                      !  spectroscopic name of angular momentum

! counters

  integer   ::  nt, n

  data lname / 's', 'p', 'd', 'f' /


  write(6,*)
  write(6,*) '  The file hsr2.dat has the modified tight-binding data of'
  write(6,*)
  write(6,'(3x,a50)') title
  write(6,'(3x,a140)') subtitle
  write(6,*)
  write(6,'("   atoms:",10(3x,i3,1x,a2))') (natom(nt),nameat(nt),nt=1,ntype)
  write(6,'("   total charge: ",f10.3,"   cutoff: ",f10.3," hartree")')  &
       ztot, emax
  do nt = 1,ntype
    write(6,'("   basis for ",a2,":",i4," orbitals   ",40(1x,a1))')      &
         nameat(nt), norbat(1,nt),                                       &
         (lname(min(max(lorb(n,1,nt),0),3)),n=1,norbat(1,nt))
  enddo
  write(6,'("   number of basis functions:",i6)') nband
  write(6,'("   k-point grid:",3i4)') mp_grid
  if(lso == 1) then
    write(6,*) '  with spin-orbit'
  else
    write(6,*) '  without spin-orbit'
  endif
  if(loptical == 1) then
    write(6,*) '  with optical matrix elements'
  else
    write(6,*) '  without optical matrix elements'
  endif
  if(imtb == 2) then
    write(6,*) '  matching of atomic orbital and plane-wave states'
  else
    write(6,*) '  projection of the plane-wave states'
  endif
  write(6,*)

  return

end subroutine ao_int_hsr2_print
