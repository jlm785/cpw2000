!> writes the pwscf.in file

subroutine write_pwscf_in(meta_pwdat,                                    &
    adot, ntype, natom, nameat, rat, alatt,                              &
    emax, nbandin, nx,ny,nz, sx,sy,sz,                                   &
    mxdtyp, mxdatm)

! Written 15 October 2018.
! Modernized 12 February 2021. JLM
! copyright  J.L.Martins, INESC-MN.

! version 4.99

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

! local:
  
  real(REAL64)          ::  avec(3,3),bvec(3,3)
  integer               ::  nat
  integer               ::  io                                           !  tape number
  real(REAL64)          ::  atmass
  character(len=80)     ::  pseudofile

! counters

  integer       ::  nt, i, j


! open file

  io = 10 
  open(unit = io, file = 'pwscf.in',status='UNKNOWN', form='FORMATTED')

! Finds total number of atoms

  nat = 0
  do i=1,ntype
    nat = nat + natom(i)
  enddo

! lattice vectors
  
  call adot_to_avec_sym(adot,avec,bvec)

  write(io,'("&CONTROL")')
  write(io,'("  calculation = ''scf'',")')
  write(io,'("   pseudo_dir = ''.'',")')
  write(io,'("        title = ''",a50,"'' ,")') meta_pwdat(1:50)
  write(io,'("/")')

  write(io,'("&SYSTEM")')
  write(io,'("  ibrav     = 0,")')
  write(io,'("  celldm(1) = ",f20.10," ,")') alatt
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

  write(io,'("&CELL")')
  write(io,'("/")')
  
  write(io,'("CELL_PARAMETERS {alat}")')
  write(io,'(3(3x,f16.8))') avec(1,1)/alatt,avec(2,1)/alatt,avec(3,1)/alatt
  write(io,'(3(3x,f16.8))') avec(1,2)/alatt,avec(2,2)/alatt,avec(3,2)/alatt
  write(io,'(3(3x,f16.8))') avec(1,3)/alatt,avec(2,3)/alatt,avec(3,3)/alatt

  write(io,'("ATOMIC_SPECIES")')
  do i = 1,ntype
    call p_tbl_mass(nameat(i),atmass)
    pseudofile = adjustl(trim(nameat(i)))//'_LDA_ncpp.UPF'
    write(io,'(a2,5x,f12.6,5x,a20)') nameat(i),atmass, trim(pseudofile)
  enddo
  
  write(io,'("ATOMIC_POSITIONS {alat} ")')
  do nt = 1,ntype
  do i = 1,natom(nt)
     write(io,'(3x,a2,5x,3f16.8)')  nameat(nt),(rat(j,i,nt),j=1,3)
  enddo
  enddo

  write(io,'("K_POINTS {automatic}")')
  write(io,'(6(2x,i6))') nx,ny,nz,nint(2*sx),nint(2*sy),nint(2*sz)
  
  close(unit = io)

  return
end subroutine write_pwscf_in

