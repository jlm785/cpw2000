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

subroutine ao_interpolation_write_wt_in(mtb,ztot,adot,ntype,natom,nameat,rat, norbat, lorb, mxdtyp, mxdatm, mxdlao)
  use NonOrthoInterp

  implicit none

  integer, parameter            :: REAL64 = selected_real_kind(12)

  type (noiData_t) :: mtb
  real(REAL64), intent(in)      :: ztot

  real(REAL64), intent(in)      :: adot(3,3)
  integer, intent(in)           :: ntype
  integer, intent(in)           :: natom(mxdtyp)
  character(len=2)              :: nameat(mxdtyp)
  real(REAL64), intent(in)      :: rat(3,mxdatm,mxdtyp)
  integer, intent(in)           :: norbat(mxdtyp)
  integer, intent(in)           :: lorb(mxdlao, mxdtyp)
  integer, intent(in)           :: mxdtyp
  integer, intent(in)           :: mxdatm
  integer, intent(in)           :: mxdlao

  real(REAL64) :: avec(3,3), bvec(3,3)
  integer :: i,j,k, ntotal_atoms, nproj, numoccupied

  character(len=40) :: hrfile
  logical           :: BulkBand_calc
  integer           :: soc
  real(REAL64)      :: e_fermi

  character(len= 6) :: orb_name(9)


  namelist /tb_file/ hrfile
  namelist /control/ BulkBand_calc
  namelist /system/ soc, e_fermi, numoccupied

  integer nlines,ndum
  real(REAL64) k_start(3), k_end(3)
  character(len=6) label_start, label_mid, label_end

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

  do i=1, 9
    write(*,*) 'orb', i, orb_name(i)
  enddo


  open(unit=22, file="wt.in", form="formatted")

  hrfile  = "mtb_hr.dat"
  BulkBand_calc = .true.

! input needed here

  if(mtb%lso == 1) then
    soc = 1
    numoccupied = nint(ztot)
  else
    soc = 0
    numoccupied = nint(ztot)/2
  endif

  e_fermi = mtb %eref * 27.21138505_REAL64

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
  do i=1, ntype
  do j=1, natom(ntype)
    ntotal_atoms = ntotal_atoms + 1
  enddo
  enddo

  write(22,'("ATOM_POSITIONS")')
  write(22,'(i5)') ntotal_atoms
  write(22,'("Direct")')
  do i=1, ntype
  do j=1, natom(ntype)
      write(22,'(a2, 3f14.8)') nameat(i), rat(:,j,i)
  enddo
  enddo

  write(22,*)

  write(22,'("PROJECTORS")')
!  write(22,'(50i5)') ((norbat(i), i=1, ntype),
  do i=1, ntype
  do j=1, natom(ntype)

    nproj = 0
    do k = 1,  norbat(i)
      nproj = nproj + (2*lorb(k,i) +1)
    enddo

    write(22,'(i5)', advance='no') nproj
  enddo
  enddo


  write(22,*)

  do i=1, ntype
  do j=1, natom(ntype)
    write(22,'(a2,a36)') nameat(i), "  s px py pz dxy dyz dzx dx2-y2 dz2"
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

  do i=1, nlines
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
!end program
end subroutine
