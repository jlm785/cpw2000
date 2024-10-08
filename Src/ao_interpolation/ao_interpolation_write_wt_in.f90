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

!>  Writes the input file for WannierTools (wt.in)
!>
!>  \author       Carlos Loia Reis
!>  \version      5.11
!>  \date         2020, 7 October 2024.
!>  \copyright    GNU Public License v2



subroutine ao_interpolation_write_wt_in(mtb,                             &
    ztot, adot, ntype, natom, nameat, rat, norbat, lorb,                 &
    mxdtyp, mxdatm, mxdlao)

! Written by Carlos Lois Reis at an unknown date.
! Documentation, indentation, 7 October 2024. JLM

  use NonOrthoInterp

  implicit none

  integer, parameter            :: REAL64 = selected_real_kind(12)

  type (noiData_t) :: mtb

! input


  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type
  integer, intent(in)                ::  mxdlao                          !<  array dimension of orbital per atom type

  real(REAL64), intent(in)           ::  ztot                            !<  total charge density (electrons/cell)


  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  character(len=2)                   ::  nameat(mxdtyp)                  !<  chemical symbol for the type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i

  integer, intent(in)                ::  norbat(mxdtyp)                  !<  number of atomic orbitals for atom k
  integer, intent(in)                ::  lorb(mxdlao,mxdtyp)             !<  angular momentum of orbital n of atom k

! local variables

  real(REAL64)         ::  avec(3,3), bvec(3,3)
  integer              ::  ntotal_atoms, nproj, numoccupied
  integer              ::  itype, iatom, k, iorb, iline

  character(len=40)    ::  hrfile
  logical              ::  BulkBand_calc
  integer              ::  soc
  real(REAL64)         ::  e_fermi

  character(len= 6)    ::  orb_name(9)

  integer              ::  nlines,ndum
  real(REAL64)         ::  k_start(3), k_end(3)
  character(len=6)     ::  label_start, label_mid, label_end

! constants

  real(REAL64), parameter     ::  EV = 27.21138505_REAL64

! namelists   I do not like them...

  namelist /tb_file/ hrfile
  namelist /control/ BulkBand_calc
  namelist /system/ soc, e_fermi, numoccupied

!
  orb_name(1) = "s"
  orb_name(2) = "px"
  orb_name(3) = "py"
  orb_name(4) = "pz"
  orb_name(5) = "dxy"
  orb_name(6) = "dyz"
  orb_name(7) = "dzx"
  orb_name(8) = "dx2-y2"
  orb_name(9) = "dz2"

  do iorb = 1,9
    write(6,*) 'orb', iorb, orb_name(iorb)
  enddo

  open(unit=22, file="wt.in", form="formatted")

  hrfile  = "mtb_hr.dat"
  BulkBand_calc = .TRUE.

! input needed here

  if(mtb%lso == 1) then
    soc = 1
    numoccupied = nint(ztot)
  else
    soc = 0
    numoccupied = nint(ztot)/2
  endif

  e_fermi = mtb%eref * EV

  write(22,nml=tb_file)
  write(22,nml=control)
  write(22,nml=system)

  call adot_to_avec_sym(adot,avec,bvec)

  write(22,'("LATTICE")')
  write(22,'("Bohr")')

  write(22,'(3f14.8)') (avec(k,1),k=1,3)
  write(22,'(3f14.8)') (avec(k,2),k=1,3)
  write(22,'(3f14.8)') (avec(k,3),k=1,3)

  write(22,*)

  ntotal_atoms = 0
  do itype=1, ntype
  do iatom=1, natom(itype)
    ntotal_atoms = ntotal_atoms + 1
  enddo
  enddo

  write(22,'("ATOM_POSITIONS")')
  write(22,'(i5)') ntotal_atoms
  write(22,'("Direct")')
  do itype=1, ntype
  do iatom=1, natom(itype)
      write(22,'(a2, 3f14.8)') nameat(itype), rat(:,iatom,itype)
  enddo
  enddo

  write(22,*)

  write(22,'("PROJECTORS")')
!  write(22,'(50i5)') ((norbat(i), i=1, ntype),
  do itype=1, ntype
  do iatom=1, natom(itype)

    nproj = 0
    do k = 1,  norbat(itype)
      nproj = nproj + (2*lorb(k,itype) +1)
    enddo

    write(22,'(i5)', advance='no') nproj
  enddo
  enddo


  write(22,*)

  do itype=1, ntype
  do iatom=1, natom(itype)
    write(22,'(a2,a36)') nameat(itype), "  s px py pz dxy dyz dzx dx2-y2 dz2"
  enddo
  enddo


  label_start='      '
  label_mid  ='      '
  label_end  ='      '

! band lines
  open(unit=44, file="BAND_LINES.DAT", form="formatted")
  read(44,*) nlines, ndum

  write(22,*)
  write(22,'("KPATH_BULK")')
  write(22,*) nlines

  do iline=1, nlines
    read(44,*) k_start(1),k_start(2),k_start(3), k_end(1),k_end(2),k_end(3), ndum, label_start, label_mid, label_end
    write(22,'(a6,2x,3f14.8,2x,a6,2x, 3f14.8)') label_start, k_start, label_end, k_end
  enddo

  close(unit=44)

  write(22,*)
  write(22,'("SURFACE")')
  write(22,'("1 0 0 ")')
  write(22,'("0 1 0 ")')
  write(22,*)

  close(unit=22)

  return

end subroutine ao_interpolation_write_wt_in
