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

!>  Reads the PW_RHO_V.DAT file and writes the corresponding xsf file
!>
!>  \author       Jose Luis Martins
!>  \version      5.13
!>  \date         5 December 2016, 3 October 2026.
!>  \copyright    GNU Public License v2

program pw_rho_v_2_xsf

! writen December 5, 2016.jlm
! Modified, write_xsf, 31 January 2021
! Constants updated to CODATA 2022 (BOHR). 3 October 2026. JLM+claude

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

  real(REAL64), parameter  :: BOHR = 0.5291772105_REAL64

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

