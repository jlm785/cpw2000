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

!>  Reads from the input file the choice of atomic basis set
!>  for each type of atom (block ChoiceOfDefaultAtomicBasisSet).
!>  Each line of the block has the species index (as in
!>  Chemical_Species_Label) and the chosen basis set.
!>  Species not listed, or no block, use the first basis set.
!>  The basis set index is not checked here, as the pseudopotentials
!>  have not yet been read.  An unavailable basis set is later
!>  replaced by the nearest available one (atomic_orbital_c16).
!>
!>  \author       claude
!>  \version      5.13
!>  \date         18 October 2026.
!>  \copyright    GNU Public License v2

subroutine read_esdf_atorb_choice(ipr, fname, ntype, nameat, ibset,      &
    mxdtyp)

! Written 18 October 2026. claude

  use esdf

  implicit none

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms

  integer, intent(in)                ::  ipr                             !<  print control: ipr < 1 no printing, except warnings.
  character(len=*), intent(in)       ::  fname                           !<  name of the input file (cpw.in)

  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  character(len=2), intent(in)       ::  nameat(mxdtyp)                  !<  chemical symbol for the type i

! output

  integer, intent(out)               ::  ibset(mxdtyp)                   !<  basis set chosen for atom type nt

! local variables

  integer          ::  nlines, ioerr
  integer          ::  nt, ib
  logical          ::  lblock

! counters

  integer          ::  i


  do nt = 1,mxdtyp
    ibset(nt) = 1
  enddo

  call esdf_init(fname)

  lblock = esdf_block('ChoiceOfDefaultAtomicBasisSet',nlines)

  if(lblock) then

    do i = 1,nlines

      read(block_data(i),*,iostat=ioerr) nt, ib

      if(ioerr /= 0) then
        write(6,*)
        write(6,*) '    Stopped in read_esdf_atorb_choice,'
        write(6,*) '    error reading ChoiceOfDefaultAtomicBasisSet,'
        write(6,*) '    line: ',trim(block_data(i))

        stop

      endif

      if(nt < 1 .or. nt > ntype) then
        write(6,*)
        write(6,*) '   WARNING   '
        write(6,*)
        write(6,'("   In ChoiceOfDefaultAtomicBasisSet species index ",  &
           &      i4," is outside 1...",i4,", line ignored")') nt, ntype
      else
        ibset(nt) = ib
      endif

    enddo

  endif

  call esdf_close

  if(lblock .and. ipr > 0) then
    write(6,*)
    write(6,'("  Choice of default atomic basis set (an unavailable ",   &
        &     "one is replaced by the nearest available):")')
    write(6,*)
    do nt = 1,ntype
      write(6,'(6x,i4,2x,a2,"   basis set ",i4)') nt, nameat(nt), ibset(nt)
    enddo
    write(6,*)
  endif

  return

end subroutine read_esdf_atorb_choice
