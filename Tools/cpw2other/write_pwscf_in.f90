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
!>  \version      5.04
!>  \date         15 October 2018, 23 january 2022.
!>  \copyright    GNU Public License v2

subroutine write_pwscf_in(meta_pwdat,                                    &
    adot, ntype, natom, nameat, rat, alatt,                              &
    emax, nbandin, nx,ny,nz, sx,sy,sz,                                   &
    mxdtyp, mxdatm)

! Written 15 October 2018.
! Modernized 12 February 2021. JLM
! Bravais lattice for phonon calculations. 23 january 2022. JLM

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

  real(REAL64), allocatable          ::  rat_car(:,:,:), rat_qe(:,:)

! metric variables

  real(REAL64)                       ::  adotnig(3,3), adotsym(3,3)
  integer                            ::  ibravais, mtotal(3,3)
  real(REAL64)                       ::  avec(3,3), bvec(3,3)
  real(REAL64)                       ::  aconv(3,3), avecnig(3,3)

  real(REAL64)                       ::  avecqe(3,3), bvecqe(3,3)
  real(REAL64)                       ::  vcellqe
  integer                            ::  iqe

! local:

  integer               ::  nat
  integer               ::  io                                           !  tape number
  real(REAL64)          ::  atmass
  character(len=80)     ::  pseudofile
  integer               ::  ipr

! QE stuff

  integer               ::  ibrqe
  real(REAL64)          ::  aa, bb, cc, cosbc, cosac, cosab

! parameter

  real(REAL64), parameter ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter ::  TOL = 1.0E-6_REAL64

! counters

  integer       ::  nt, i, j, k


! open file

  io = 10
  open(unit = io, file = 'pwscf.in',status='UNKNOWN', form='FORMATTED')


! lattice vectors

  call adot_to_avec_sym(adot,avec,bvec)
  call metric_ident(adot, adotnig, adotsym, ibravais, mtotal,            &
                      avec, aconv, avecnig, TOL, ipr)

! Finds total number of atoms, converts to cartesian

  nat = 0
  do nt=1,ntype
    nat = nat + natom(nt)
  enddo

  allocate(rat_car(3,mxdatm,mxdtyp))

  do nt=1,ntype
  do j = 1,natom(nt)
  do k = 1,3
    rat_car(k,j,nt) = avec(k,1)*rat(1,j,nt) + avec(k,2)*rat(2,j,nt)      &
                                            + avec(k,3)*rat(3,j,nt)
  enddo
  enddo
  enddo


! converts convention of bravais lattice

  iqe = 0
  aa = alatt
  if(ibravais == 1) then

    iqe = 1
    aa = aconv(1,1)

    avecqe(1,1) = aa
    avecqe(2,1) = ZERO
    avecqe(3,1) = ZERO

    avecqe(1,2) = ZERO
    avecqe(2,2) = aa
    avecqe(3,2) = ZERO

    avecqe(1,3) = ZERO
    avecqe(2,3) = ZERO
    avecqe(3,3) = aa

  elseif(ibravais == 2) then

    iqe = 2
    aa = aconv(1,1)

    avecqe(1,1) =-aa / 2
    avecqe(2,1) = ZERO
    avecqe(3,1) = aa / 2

    avecqe(1,2) = ZERO
    avecqe(2,2) = aa / 2
    avecqe(3,2) = aa / 2

    avecqe(1,3) =-aa / 2
    avecqe(2,3) = aa / 2
    avecqe(3,3) = ZERO

  elseif(ibravais == 3) then

    iqe = 3
    aa = aconv(1,1)

    avecqe(1,1) = aa / 2
    avecqe(2,1) = aa / 2
    avecqe(3,1) = aa / 2

    avecqe(1,2) =-aa / 2
    avecqe(2,2) = aa / 2
    avecqe(3,2) = aa / 2

    avecqe(1,3) =-aa / 2
    avecqe(2,3) =-aa / 2
    avecqe(3,3) = aa / 2

  elseif(ibravais == 4) then

    iqe = 6
    aa = aconv(1,1)
    cc = aconv(3,3)

    avecqe(1,1) = aa
    avecqe(2,1) = ZERO
    avecqe(3,1) = ZERO

    avecqe(1,2) = ZERO
    avecqe(2,2) = aa
    avecqe(3,2) = ZERO

    avecqe(1,3) = ZERO
    avecqe(2,3) = ZERO
    avecqe(3,3) = cc

  elseif(ibravais == 5) then

    iqe = 7
    aa = aconv(1,1)
    cc = aconv(3,3)

    avecqe(1,1) = aa / 2
    avecqe(2,1) =-aa / 2
    avecqe(3,1) = cc / 2

    avecqe(1,2) = aa / 2
    avecqe(2,2) = aa / 2
    avecqe(3,2) = cc / 2

    avecqe(1,3) =-aa / 2
    avecqe(2,3) =-aa / 2
    avecqe(3,3) = cc / 2

  elseif(ibravais == 6) then

    iqe = 8
    aa = aconv(1,1)
    bb = aconv(2,2)
    cc = aconv(3,3)

    avecqe(1,1) = aa
    avecqe(2,1) = ZERO
    avecqe(3,1) = ZERO

    avecqe(1,2) = ZERO
    avecqe(2,2) = bb
    avecqe(3,2) = ZERO

    avecqe(1,3) = ZERO
    avecqe(2,3) = ZERO
    avecqe(3,3) = cc

  elseif(ibravais == 7) then

    iqe = 9
    aa = aconv(1,1)
    bb = aconv(2,2)
    cc = aconv(3,3)

    avecqe(1,1) = aa / 2
    avecqe(2,1) = bb / 2
    avecqe(3,1) = ZERO

    avecqe(1,2) =-aa / 2
    avecqe(2,2) = bb / 2
    avecqe(3,2) = ZERO

    avecqe(1,3) = ZERO
    avecqe(2,3) = ZERO
    avecqe(3,3) = cc

  elseif(ibravais == 8) then

    iqe = 10
    aa = aconv(1,1)
    bb = aconv(2,2)
    cc = aconv(3,3)

    avecqe(1,1) = aa / 2
    avecqe(2,1) = ZERO
    avecqe(3,1) = cc / 2

    avecqe(1,2) = aa / 2
    avecqe(2,2) = bb / 2
    avecqe(3,2) = ZERO

    avecqe(1,3) = ZERO
    avecqe(2,3) = bb / 2
    avecqe(3,3) = cc / 2

  elseif(ibravais == 9) then

    iqe = 11
    aa = aconv(1,1)
    bb = aconv(2,2)
    cc = aconv(3,3)

    avecqe(1,1) = aa / 2
    avecqe(2,1) = bb / 2
    avecqe(3,1) = cc / 2

    avecqe(1,2) =-aa / 2
    avecqe(2,2) = bb / 2
    avecqe(3,2) = cc / 2

    avecqe(1,3) =-aa / 2
    avecqe(2,3) =-bb / 2
    avecqe(3,3) = cc / 2

  elseif(ibravais == 10) then

    iqe = 4
    aa = aconv(1,1)
    cc = aconv(3,3)

    avecqe(1,1) = aa
    avecqe(2,1) = ZERO
    avecqe(3,1) = ZERO

    avecqe(1,2) =-aa / 2
    avecqe(2,2) = aa*sqrt(3*UM) / 2
    avecqe(3,2) = ZERO

    avecqe(1,3) = ZERO
    avecqe(2,3) = ZERO
    avecqe(3,3) = cc

  elseif(ibravais == 11) then

    iqe = 5
    aa = sqrt((3*aconv(1,1)*aconv(1,1)+aconv(3,3)*aconv(3,3)) / (9*UM))
    cosbc = (2*aconv(3,3)*aconv(3,3)-3*aconv(1,1)*aconv(1,1)) / (18*aa*aa)

    avecqe(1,1) = aa*sqrt(UM - cosbc)/sqrt(2*UM)
    avecqe(2,1) =-aa*sqrt(UM - cosbc)/sqrt(6*UM)
    avecqe(3,1) = aa*sqrt(UM + 2*cosbc)/sqrt(3*UM)

    avecqe(1,2) = ZERO
    avecqe(2,2) = 2*aa*sqrt(UM - cosbc)/sqrt(6*UM)
    avecqe(3,2) = aa*sqrt(UM + 2*cosbc)/sqrt(3*UM)

    avecqe(1,3) =-aa*sqrt(UM - cosbc)/sqrt(2*UM)
    avecqe(2,3) =-aa*sqrt(UM - cosbc)/sqrt(6*UM)
    avecqe(3,3) = aa*sqrt(UM + 2*cosbc)/sqrt(3*UM)

  elseif(ibravais == 12) then

    iqe = -12
    aa = aconv(1,1)
    bb = aconv(2,2)
    cc = sqrt(aconv(1,3)*aconv(1,3)+aconv(3,3)*aconv(3,3))
    cosac = aconv(1,3) / cc

    avecqe(1,1) = aa
    avecqe(2,1) = ZERO
    avecqe(3,1) = ZERO

    avecqe(1,2) = ZERO
    avecqe(2,2) = bb
    avecqe(3,2) = ZERO

    avecqe(1,3) = aconv(1,3)
    avecqe(2,3) = ZERO
    avecqe(3,3) = aconv(3,3)

  elseif(ibravais == 13) then

    iqe = -13
    aa = aconv(1,1)
    bb = aconv(2,2)
    cc = sqrt(aconv(1,3)*aconv(1,3)+aconv(3,3)*aconv(3,3))
    cosac = aconv(1,3) / cc

    avecqe(1,1) = aa / 2
    avecqe(2,1) = bb / 2
    avecqe(3,1) = ZERO

    avecqe(1,2) =-aa / 2
    avecqe(2,2) = bb / 2
    avecqe(3,2) = ZERO

    avecqe(1,3) = aconv(1,3)
    avecqe(2,3) = ZERO
    avecqe(3,3) = aconv(3,3)

  elseif(ibravais == 14) then

    iqe = 14
    aa = aconv(1,1)
    bb = sqrt(aconv(1,2)*aconv(1,2)+aconv(2,2)*aconv(2,2))
    cc = sqrt(aconv(1,3)*aconv(1,3)+aconv(2,3)*aconv(2,3)+aconv(3,3)*aconv(3,3))
    cosbc = aconv(2,3) / cc
    cosac = aconv(1,3) / cc
    cosab = aconv(1,2) / bb

    avecqe(1,1) = aconv(1,1)
    avecqe(2,1) = ZERO
    avecqe(3,1) = ZERO

    avecqe(1,2) = aconv(1,2)
    avecqe(2,2) = aconv(2,2)
    avecqe(3,2) = ZERO

    avecqe(1,3) = aconv(1,3)
    avecqe(2,3) = aconv(2,3)
    avecqe(3,3) = aconv(3,3)

  endif

! converts cartesian atomic positions to QE lattice coordinates

  allocate(rat_qe(3,nat))

  bvecqe(1,1) = avecqe(2,2)*avecqe(3,3) - avecqe(3,2)*avecqe(2,3)
  bvecqe(2,1) = avecqe(3,2)*avecqe(1,3) - avecqe(1,2)*avecqe(3,3)
  bvecqe(3,1) = avecqe(1,2)*avecqe(2,3) - avecqe(2,2)*avecqe(1,3)
  bvecqe(1,2) = avecqe(2,3)*avecqe(3,1) - avecqe(3,3)*avecqe(2,1)
  bvecqe(2,2) = avecqe(3,3)*avecqe(1,1) - avecqe(1,3)*avecqe(3,1)
  bvecqe(3,2) = avecqe(1,3)*avecqe(2,1) - avecqe(2,3)*avecqe(1,1)
  bvecqe(1,3) = avecqe(2,1)*avecqe(3,2) - avecqe(3,1)*avecqe(2,2)
  bvecqe(2,3) = avecqe(3,1)*avecqe(1,2) - avecqe(1,1)*avecqe(3,2)
  bvecqe(3,3) = avecqe(1,1)*avecqe(2,2) - avecqe(2,1)*avecqe(1,2)

  vcellqe = bvecqe(1,1)*avecqe(1,1) + bvecqe(2,1)*avecqe(2,1) +          &
            bvecqe(3,1)*avecqe(3,1)

  do j = 1,3
  do k = 1,3
    bvecqe(k,j) = bvecqe(k,j) / vcellqe
  enddo
  enddo

  i = 0
  do nt=1,ntype
  do j = 1,natom(nt)
    i  = i+1
    do k = 1,3
      rat_qe(k,i) = bvecqe(1,k)*rat_car(1,j,nt) + bvecqe(2,k)*rat_car(2,j,nt)      &
                                                + bvecqe(3,k)*rat_car(3,j,nt)
    enddo
  enddo
  enddo

  write(io,'("&CONTROL")')
  write(io,'("  calculation = ''scf'',")')
  write(io,'("   pseudo_dir = ''.'',")')
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

  deallocate(rat_car,rat_qe)

  return
end subroutine write_pwscf_in

