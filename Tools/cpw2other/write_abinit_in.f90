!>  writes the abinit.in file

subroutine write_abinit_in(meta_pwdat,                                   &
  adot,ntype,natom,nameat,rat,alatt,                                     &
  emax,nbandin,nx,ny,nz,sx,sy,sz,                                        &
  mxdtyp,mxdatm)

! Written 30 October 2018.
! Documentation, 12 february 2021. JLM
! copyright  J.L.Martins, INESC-MN.

! version 4.99

  implicit none

  integer, parameter  :: REAL64 = selected_real_kind(12)

! input:

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of types of atoms

  character(len=250), intent(in)     ::  meta_pwdat                      !<  metadata from cpw_in or PW.DAT

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  character(len=2), intent(in)       ::  nameat(mxdtyp)                  !<  chemical symbol for the type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i

  real(REAL64), intent(in)           ::  alatt                           !<  lattice constant

  real(REAL64), intent(in)           ::  emax                            !<  kinetic energy cutoff of plane wave expansion (Hartree).
  integer, intent(in)                ::  nbandin                         !<  target for number of bands      

  integer, intent(in)                ::  nx, ny, nz                      !<  size of the integration mesh in k-space (nx*ny*nz)
  real(REAL64), intent(in)           ::  sx, sy, sz                      !<  offset of the integration mesh (usually 0.5)

! local:
  
  real(REAL64)          ::  avec(3,3),bvec(3,3)
  integer               ::  nat
  integer               ::  io                                           !  tape number
  real(REAL64)          ::  atmass
  character(len=80)     ::  pseudofile

  integer, allocatable  ::  iznuc(:)
  

! counters

  integer       ::  nt, i, j


! open file

  io = 10 
  open(unit = io, file = 'abinit.in',status='UNKNOWN', form='FORMATTED')

! Finds total number of atoms

  nat = 0
  do i=1,ntype
    nat = nat + natom(i)
  enddo

! lattice vectors
  
  call adot_to_avec_sym(adot,avec,bvec)

! nuclear charge

  allocate(iznuc(mxdtyp))
  do nt = 1,ntype
    call p_tbl_charge(nameat(nt),iznuc(nt))
  enddo

  write(io,*) '#Definition of the unit cell'
  write(io,'("acell",3(4x,f16.10))') alatt,alatt,alatt
  write(io,*) 'rprim'
  write(io,'(3(3x,f16.8))') avec(1,1)/alatt,avec(2,1)/alatt,avec(3,1)/alatt
  write(io,'(3(3x,f16.8))') avec(1,2)/alatt,avec(2,2)/alatt,avec(3,2)/alatt
  write(io,'(3(3x,f16.8))') avec(1,3)/alatt,avec(2,3)/alatt,avec(3,3)/alatt
  write(io,*)
  write(io,*)

  write(io,*) '#Definition of the atom types'
  write(io,'("ntypat",4x,i5)') ntype
  write(io,'("znucl",20(5x,i5))') (iznuc(nt),nt=1,ntype)
  write(io,*)
  write(io,*)

  write(io,*) '#Definition of the atoms'
  write(io,'("natom",4x,i5)') nat
!  write(io,'("typat",(200(5x,i5)))') ((nt,j=1,natom(nt)),nt=1,ntype)
  write(io,'("typat")')
  do nt = 1,ntype
    write(io,'(20(5x,i5))') (nt,j=1,natom(nt))
  enddo
  write(io,'("xred")')
  do nt = 1,ntype
  do i = 1,natom(nt)
    write(io,'(3(5x,f15.10))') (rat(j,i,nt),j=1,3)
  enddo
  enddo
  write(io,*)
  write(io,*)

  write(io,*) '#Definition of the planewave basis set'
  write(io,'("ecut",5x,f20.15)') emax
  write(io,*)
  write(io,*)

  write(io,*) '#Definition of the k-point grid'
  write(io,'("kptopt",4x,i5)') 1
  write(io,'("chksymbreak",2x,i2)') 0
  write(io,'("ngkpt",3(5x,i5))')  nx,ny,nz
  write(io,'("nshiftk",3x,i5)') 1
  write(io,*)
  write(io,*)

  write(io,*) '#Definition of the SCF procedure'
  write(io,'("nstep",5x,i5)') 50
  write(io,'("toldfe",5x,e15.5)') 1.0e-6
  write(io,'("diemac",5x,f12.5)') 1.0

  close(unit = io)

  deallocate(iznuc)

  return
end subroutine write_abinit_in

