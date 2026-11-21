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
!>
!>  \author       Jose Luis Martins
!>  \version      5.13
!>  \date         15 October 2018, 20 November 2026.
!>  \copyright    GNU Public License v2

subroutine pw2o_qe_pwscf_in(meta_pwdat,                                  &
    adot, ntype, natom, nameat, rat, alatt,                              &
    emax, nbandin, mpgrid, mpshift,                                      &
    mxdtyp, mxdatm)

! Written 15 October 2018.
! Modernized 12 February 2021. JLM
! Bravais lattice for phonon calculations. 23 january 2022. JLM
! Removed unused variables, 14 November 2024. JLM
! Extracted pw2o_qe_convert_crys_struct. 25 February 2025. JLM
! Documentation, one argument per declaration. 28 September 2026. JLM+claude
! mpgrid, mpshift, explicit list of k-points for a general mesh. 20 November 2026. JLM+claude WARNING NEW API
! K_POINTS automatic only if the QE and cpw lattice vectors are the same. 20 November 2026. JLM+claude


! A diagonal Monkhorst-Pack mesh with shifts 0 or 1/2 is written
! as K_POINTS automatic if the QE lattice vectors are the same as
! the cpw lattice vectors.  Other meshes are written as an explicit
! list (K_POINTS crystal, in the QE lattice vectors) of the points
! in the irreducible wedge of the lattice point group (holohedry),
! as QE expects, with weights normalized to one.



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

  integer, intent(in)                ::  mpgrid(3,3)                     !< supercell of the Monkhorst-Pack mesh,  A_i = sum_j mpgrid(i,j) a_j  (line i of kgrid_Monkhorst_Pack)
  real(REAL64), intent(in)           ::  mpshift(3)                      !< shift of the mesh in units of the reciprocal vectors of the supercell

! allocatable array

  real(REAL64), allocatable          ::  rat_qe(:,:)

  real(REAL64), allocatable          ::  rk(:,:)                         !  k-points in the irreducible wedge of the holohedry (cpw lattice coordinates)
  real(REAL64), allocatable          ::  wgk(:)                          !  weights of the k-points
  integer, allocatable               ::  nband(:)                        !  not used
  integer, allocatable               ::  indk(:,:)                       !  not used
  integer, allocatable               ::  kmap(:,:,:,:)                   !  not used

! k-point mesh

  integer                            ::  nx, ny, nz                      !  box of the mesh (see kgrid_box)
  real(REAL64)                       ::  sx, sy, sz                      !  shifts of the box
  integer                            ::  kmat(3,3)                       !  rk = kmat (((i-1)+sx)/nx, ...)
  logical                            ::  lauto                           !  mesh can be written as K_POINTS automatic
  integer                            ::  mxdnrk, nrk
  integer                            ::  norbit                          !  sum of the sizes of the orbits of the irreducible points
  real(REAL64)                       ::  sumw
  real(REAL64)                       ::  tqe(3,3)                        !  QE lattice vectors in cpw lattice coordinates
  real(REAL64)                       ::  avinv(3,3)                      !  inverse of avec
  real(REAL64)                       ::  vcell
  real(REAL64)                       ::  xkqe(3)

! lattice point group (holohedry)

  integer                            ::  ntrh                            !  number of operations of the holohedry
  integer                            ::  mtrxh(3,3,48)                   !  operations in reciprocal lattice coordinates
  real(REAL64)                       ::  tnph(3,48)                      !  zero
  real(REAL64)                       ::  adoth(3,3)                      !  copy of adot (sym_identify symmetrizes it)
  integer                            ::  nath(1)                         !  one atom ...
  real(REAL64)                       ::  rath(3,1,1)                     !  ... at the origin

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

  integer       ::  nt, i, j, k, n

! parameters

  real(REAL64), parameter ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter ::  EPS = 1.0E-8_REAL64
  real(REAL64), parameter ::  TOLSYM = 1.0E-5_REAL64
  real(REAL64), parameter ::  TOLT = 1.0E-4_REAL64


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
  write(io,'("       prefix = ''cpw2qe'',")')
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

  write(io,*)


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

! k-points

  call kgrid_box(mpgrid, mpshift, nx, ny, nz, sx, sy, sz, kmat)

! QE lattice vectors in cpw lattice coordinates, tqe = avec^-1 avec_qe

  avinv(1,1) = avec(2,2)*avec(3,3) - avec(3,2)*avec(2,3)
  avinv(1,2) = avec(3,2)*avec(1,3) - avec(1,2)*avec(3,3)
  avinv(1,3) = avec(1,2)*avec(2,3) - avec(2,2)*avec(1,3)
  avinv(2,1) = avec(2,3)*avec(3,1) - avec(3,3)*avec(2,1)
  avinv(2,2) = avec(3,3)*avec(1,1) - avec(1,3)*avec(3,1)
  avinv(2,3) = avec(1,3)*avec(2,1) - avec(2,3)*avec(1,1)
  avinv(3,1) = avec(2,1)*avec(3,2) - avec(3,1)*avec(2,2)
  avinv(3,2) = avec(3,1)*avec(1,2) - avec(1,1)*avec(3,2)
  avinv(3,3) = avec(1,1)*avec(2,2) - avec(2,1)*avec(1,2)

  vcell = avec(1,1)*avinv(1,1) + avec(2,1)*avinv(1,2) + avec(3,1)*avinv(1,3)

  do j = 1,3
  do i = 1,3
    tqe(i,j) = ( avinv(i,1)*avec_qe(1,j) + avinv(i,2)*avec_qe(2,j) +    &
                 avinv(i,3)*avec_qe(3,j) ) / vcell
  enddo
  enddo

! K_POINTS automatic is along the QE lattice vectors.  It is used
! if they are the cpw lattice vectors and the mesh is a usual grid.

  lauto = .TRUE.
  do j = 1,3
  do i = 1,3
    if(i /= j .and. abs(tqe(i,j)) > TOLT) lauto = .FALSE.
    if(i == j .and. abs(tqe(i,j) - UM) > TOLT) lauto = .FALSE.
  enddo
  enddo
  do j = 1,3
  do i = 1,3
    if(i /= j .and. kmat(i,j) /= 0) lauto = .FALSE.
    if(i == j .and. kmat(i,j) /= 1) lauto = .FALSE.
  enddo
  enddo
  if(abs(sx) > EPS .and. abs(sx - UM/2) > EPS) lauto = .FALSE.
  if(abs(sy) > EPS .and. abs(sy - UM/2) > EPS) lauto = .FALSE.
  if(abs(sz) > EPS .and. abs(sz - UM/2) > EPS) lauto = .FALSE.

  if(lauto) then

    write(io,'("K_POINTS {automatic}")')
    write(io,'(6(2x,i6))') nx,ny,nz,nint(2*sx),nint(2*sy),nint(2*sz)

  else

!   holohedry, lattice with one atom at the origin

    adoth = adot
    nath(1) = 1
    rath(1,1,1) = ZERO
    rath(2,1,1) = ZERO
    rath(3,1,1) = ZERO

    call sym_identify(1, 0, TOLSYM, ntrh, mtrxh, tnph,                   &
        adoth, 1, nath, rath, 1, 1)

!   irreducible k-points of the holohedry

    call size_mxdnrk(nx,ny,nz, sx,sy,sz, kmat, adot, ntrh, mtrxh,        &
        mxdnrk)

    allocate(rk(3,mxdnrk))
    allocate(wgk(mxdnrk))
    allocate(nband(mxdnrk))
    allocate(indk(6,mxdnrk))
    allocate(kmap(3,nx,ny,nz))

    call int_pnt(1, nx,ny,nz, sx,sy,sz, kmat, 0, adot, ntrh, mtrxh,      &
        nrk, rk, wgk, nband, indk, kmap, mxdnrk, 1)

!   is the mesh invariant under the holohedry?

    norbit = 0
    do n = 1,nrk
      call orbit_size(rk(:,n), ntrh, mtrxh, k)
      norbit = norbit + k
    enddo

    if(norbit /= nx*ny*nz) then
      write(6,*)
      write(6,*) '   WARNING:  the k-point mesh is not invariant under the'
      write(6,*) '   point group of the lattice.  QE will use the mesh'
      write(6,*) '   symmetrized by the lattice point group, with ',norbit
      write(6,*) '   instead of ',nx*ny*nz,' points.  cpw used the same'
      write(6,*) '   mesh only if the crystal has the symmetry of the lattice.'
      write(6,*)
    endif

    sumw = ZERO
    do n = 1,nrk
      sumw = sumw + wgk(n)
    enddo

    write(io,'("K_POINTS {crystal}")')
    write(io,'(2x,i8)') nrk
    do n = 1,nrk
      do j = 1,3
        xkqe(j) = tqe(1,j)*rk(1,n) + tqe(2,j)*rk(2,n) + tqe(3,j)*rk(3,n)
      enddo
      write(io,'(4(2x,f16.12))') (xkqe(j),j=1,3), wgk(n)/sumw
    enddo

    deallocate(rk)
    deallocate(wgk)
    deallocate(nband)
    deallocate(indk)
    deallocate(kmap)

  endif

  close(unit = io)

  deallocate(rat_qe)

  return

  contains

! number of distinct points in the orbit of k under the operations and inversion (time reversal)

  subroutine orbit_size(xk, nop, mop, norb)

    real(REAL64), intent(in)   ::  xk(3)
    integer, intent(in)        ::  nop
    integer, intent(in)        ::  mop(3,3,48)
    integer, intent(out)       ::  norb

    real(REAL64)   ::  xorb(3,96), xt(3), d
    logical        ::  lnew
    integer        ::  iop, isg, m, l

    norb = 0
    do iop = 1,nop
    do isg = 1,-1,-2
      do l = 1,3
        xt(l) = isg*(mop(l,1,iop)*xk(1) + mop(l,2,iop)*xk(2) + mop(l,3,iop)*xk(3))
      enddo
      lnew = .TRUE.
      do m = 1,norb
        d = ZERO
        do l = 1,3
          d = max(d, abs(xt(l) - xorb(l,m) - nint(xt(l) - xorb(l,m))))
        enddo
        if(d < EPS) then
          lnew = .FALSE.

          exit

        endif
      enddo
      if(lnew) then
        norb = norb + 1
        xorb(:,norb) = xt(:)
      endif
    enddo
    enddo

  end subroutine orbit_size

end subroutine pw2o_qe_pwscf_in

