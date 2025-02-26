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

!>  writes the pwscf.in file for use in QuantumEspresso
!>
!>  \author       Jose Luis Martins
!>  \version      5.11
!>  \date         15 October 2018, 25 February 2025.
!>  \copyright    GNU Public License v2

subroutine pw2o_qe_pwscf_in(meta_pwdat,                                  &
    adot, ntype, natom, nameat, rat, alatt,                              &
    emax, nbandin, nx,ny,nz, sx,sy,sz,                                   &
    mxdtyp, mxdatm)

! Written 15 October 2018.
! Modernized 12 February 2021. JLM
! Bravais lattice for phonon calculations. 23 january 2022. JLM
! Removed unused variables, 14 November 2024. JLM
! Extracted pw2o_qe_convert_crys_struct. 25 February 2025. JLM

  implicit none

  integer, parameter  :: REAL64 = selected_real_kind(12)

! input:

  integer, intent(in)                ::  mxdtyp                          !< array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !< array dimension of types of atoms

  character(len=250), intent(in)     ::  meta_pwdat                      !< metadata from cpw_in or PW.DAT

  real(REAL64), intent(in)           ::  adot(3,3)                       !< metric in direct space
  integer, intent(in)                ::  ntype                           !< number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !< number of atoms of type i
  character(len=2), intent(in)       ::  nameat(mxdtyp)                  !< chemical symbol for the type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !< k-th component (in lattice coordinates) of the position of the n-th atom of type i

  real(REAL64), intent(in)           ::  alatt                           !< lattice constant

  real(REAL64), intent(in)           ::  emax                            !< kinetic energy cutoff of plane wave expansion (Hartree).
  integer, intent(in)                ::  nbandin                         !< target for number of bands

  integer, intent(in)                ::  nx, ny, nz                      !< size of the integration mesh in k-space (nx*ny*nz)
  real(REAL64), intent(in)           ::  sx, sy, sz                      !< offset of the integration mesh (usually 0.5)

! allocatable array

  real(REAL64), allocatable          ::  rat_qe(:,:)

! metric variables

  integer                            ::  ibravais
  real(REAL64)                       ::  avec(3,3)
  real(REAL64)                       ::  avec_qe(3,3)

! local:

  integer               ::  nat, mxdnat
  integer               ::  io                                           !  tape number
  real(REAL64)          ::  atmass
  character(len=80)     ::  pseudofile

! QE stuff

  integer               ::  iqe
  real(REAL64)          ::  aa, bb, cc, cosbc, cosac, cosab

! counters

  integer       ::  nt, i, j, k


  mxdnat = 0
  do nt = 1,ntype
    mxdnat = mxdnat + natom(nt)
  enddo

  allocate(rat_qe(3,mxdnat))


! open file

  io = 10
  open(unit = io, file = 'pwscf.in',status='UNKNOWN', form='FORMATTED')

  call pw2o_qe_convert_crys_struct(nat, iqe, rat_qe, avec_qe,            &
      aa, bb, cc, cosbc, cosac, cosab, ibravais, avec,                   &
      adot, ntype, natom, rat, alatt,                                    &
      mxdtyp, mxdatm, mxdnat)

  write(io,'("&CONTROL")')
  write(io,'("  calculation = ''scf'',")')
  write(io,'("   pseudo_dir = ''.'',")')
  write(io,'("       prefix = ''cpw'',")')
  write(io,'("        title = ''",a50,"'' ,")') meta_pwdat(1:50)
  write(io,'("/")')

  write(io,'("&SYSTEM")')
  write(io,'("  ibrav     = ",i3," ,")') iqe
  write(io,'("  celldm(1) = ",f20.10," ,")') aa

  if(ibravais == 4) then
    write(io,'("  celldm(3) = ",f20.10," ,")') cc / aa
  elseif(ibravais > 4 .and. ibravais < 10) then
    write(io,'("  celldm(2) = ",f20.10," ,")') bb / aa
    write(io,'("  celldm(3) = ",f20.10," ,")') cc / aa
  elseif(ibravais == 10) then
    write(io,'("  celldm(3) = ",f20.10," ,")') cc / aa
  elseif(ibravais == 11) then
    write(io,'("  celldm(4) = ",f20.10," ,")') cosbc
  elseif(ibravais == 12 .or. ibravais == 13) then
    write(io,'("  celldm(2) = ",f20.10," ,")') bb / aa
    write(io,'("  celldm(3) = ",f20.10," ,")') cc / aa
    write(io,'("  celldm(5) = ",f20.10," ,")') cosac
  elseif(ibravais == 14) then
    write(io,'("  celldm(2) = ",f20.10," ,")') bb / aa
    write(io,'("  celldm(3) = ",f20.10," ,")') cc / aa
    write(io,'("  celldm(4) = ",f20.10," ,")') cosbc
    write(io,'("  celldm(5) = ",f20.10," ,")') cosac
    write(io,'("  celldm(6) = ",f20.10," ,")') cosab
  endif

  write(io,'("        nat = ",i5," ,")') nat
  write(io,'("       ntyp = ",i5," ,")') ntype
  write(io,'("       nbnd = ",i5," ,")') nbandin
  write(io,'("    ecutwfc = ",f14.3," ,")') 2*emax
  write(io,'("/")')

  write(io,'("&ELECTRONS")')
  write(io,'("  diagonalization = ''david'',")')
  write(io,'("/")')

  write(io,'("&IONS")')
  write(io,'("/")')


  if(iqe == 0) then
    write(io,'("&CELL")')
    write(io,'("/")')
    write(io,'("CELL_PARAMETERS {alat}")')
    write(io,'(3(3x,f16.8))') avec(1,1)/alatt,avec(2,1)/alatt,avec(3,1)/alatt
    write(io,'(3(3x,f16.8))') avec(1,2)/alatt,avec(2,2)/alatt,avec(3,2)/alatt
    write(io,'(3(3x,f16.8))') avec(1,3)/alatt,avec(2,3)/alatt,avec(3,3)/alatt
  endif

  write(io,'("ATOMIC_SPECIES")')
  do i = 1,ntype
    call p_tbl_mass(nameat(i),atmass)
    pseudofile = adjustl(trim(nameat(i)))//'_LDA_TM.UPF'
    write(io,'(a2,5x,f12.6,5x,a20)') nameat(i),atmass, trim(pseudofile)
  enddo

  write(io,'("ATOMIC_POSITIONS {crystal} ")')
  i = 0
  do nt = 1,ntype
  do j = 1,natom(nt)
    i = i+1
    write(io,'(3x,a2,5x,3f16.8)')  nameat(nt),(rat_qe(k,i),k=1,3)
  enddo
  enddo

  write(io,'("K_POINTS {automatic}")')
  write(io,'(6(2x,i6))') nx,ny,nz,nint(2*sx),nint(2*sy),nint(2*sz)

  close(unit = io)

  deallocate(rat_qe)

  return

end subroutine pw2o_qe_pwscf_in

