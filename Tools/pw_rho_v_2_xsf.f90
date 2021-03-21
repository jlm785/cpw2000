program pw_rho_v_2_xsf

! This program reads the PW_RHO_V.DAT file and writes the corresponding xsf file

! writen December 5, 2016.jlm
! Modified, write_xsf, 31 January 2021
! copyright  Jose Luis Martins/INESC-MN

! version 4.99

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! dimensions

  integer                            ::  mxdtyp                     !  array dimension of types of atoms
  integer                            ::  mxdatm                     !  array dimension of number of atoms of a given type

! main variables

  real(REAL64)                       ::  adot(3,3)                  !  metric in real space
  integer                            ::  ntype                      !  number of types of atoms
  integer,allocatable                ::  natom(:)                   !  number of atoms of type i
  character(len=2),allocatable       ::  nameat(:)                  !  chemical symbol for the type i
  real(REAL64),allocatable           ::  rat(:,:,:)                 !  k-th component (in lattice coordinates) of the position of the n-th atom of type i

  integer      ::  ng,ns
  integer      ::  iotape

! counters

  integer      :: i, j, k

! constants

  real(REAL64), parameter  :: BOHR = 0.52917721_REAL64

! reads data

  open(unit=10,file='PW_RHO_V.DAT',status='old',form='unformatted')

  read(10) ntype,ng,ns ! ,mxdl

  mxdtyp = ntype
  allocate(natom(mxdtyp))
  allocate(nameat(mxdtyp))

  read(10) (natom(i),i=1,ntype)
  read(10) ! bdate,btime
  read(10) ! author,flgscf,flgdal
  read(10) ! emaxin,teleck,nx,ny,nz,sx,sy,sz
  read(10) ! io,iostat = ioerr) line140


  read(10) ((adot(i,j),i=1,3),j=1,3)
  read(10) (nameat(i),i=1,ntype)

  mxdatm = natom(1)
  do i=1,ntype
    if(mxdatm < natom(i)) mxdatm = natom(i)
  enddo
  allocate(rat(3,mxdatm,mxdtyp))

  do i=1,ntype
    read(10) ((rat(j,k,i),j=1,3),k=1,natom(i))
  enddo

  close(unit=10)


! are you old enough to remember they were actual tapes?

  iotape = 26
  open(unit = iotape,file = 'vesta_relaxed.xsf',status='UNKNOWN',        &
                   form='FORMATTED')

  call plot_xsf_crys(iotape,.TRUE.,                                      &
     adot,ntype,natom,nameat,rat,                                        &
     mxdtyp,mxdatm)

! closes

  close(unit = iotape)

  deallocate(natom)
  deallocate(nameat)

  stop
end program pw_rho_v_2_xsf

