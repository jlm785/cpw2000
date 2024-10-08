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

!>  Calculates the atomic orbital interpolation
!>  for the density of states
!>
!>  \author       Carlos Loia Reis
!>  \version      5.11
!>  \date         July 2014, 8 October 2024.
!>  \copyright    GNU Public License v2

subroutine ao_interpolation_out_dos(noiData, ztot, adot, ntrans, mtrx)

! Adapted to use ao interpolatiion package July, 2014. clr
! modified 25 November 2015. JLM
! Modified documentation May 2020. JLM
! Modfied, new dos file format June 2021. CLR
! Modified, e_of_k allocations, 14 Janury 2024. JLM
! Modified, indentation, ao_int_, 8 October 2024. JLM

  use NonOrthoInterp

  implicit none

  integer, parameter                 :: REAL64 = selected_real_kind(12)

!      input

  type(noiData_t)                    ::  noiData                         !<  see NonOrthoInterp

  real(REAL64), intent(in)           ::  ztot                            !<  total charge density (electrons/cell)

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  integer, intent(in)                ::  ntrans                          !<  number of symmetry operations in the factor group
  integer, intent(in)                ::  mtrx(3,3,48)                    !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group

!      local allocatable arrays

  integer, allocatable               ::  kmap2(:,:,:,:)                  !  kmap(1,...) corresponding k-point, kmap(2,...) = 1 additional inversion, kmap(3,...) symmetry operation
  integer, allocatable               ::  nband2(:)                       !  number of bands for each k-points
  integer, allocatable               ::  indk2(:,:)                      !  index of the six k-points neighbouring k-point i
  real(REAL64),allocatable           ::  rk2(:,:)                        !  component in lattice coordinates of the k-point in the mesh
  real(REAL64),allocatable           ::  w2(:)                           !  weight in the integration of k-point

  real(REAL64), allocatable          ::  e_of_k(:,:)                     !  band energies of k-point in plot
  real(REAL64), allocatable          ::  e_of_k_so(:,:)                  !  spin-orbit band energies of k-point in plot
  real(REAL64), allocatable          ::  basxpsi_of_k(:,:,:)
  real(REAL64), allocatable          ::  ev_interp(:)

  complex(REAL64), allocatable      ::   psi(:,:)

!      local variables

  integer                            ::  mxdbnd                          !  array dimension for the number of bands
  integer                            ::  mxdpnt                          !  dimensions for dos k-points

  integer                            ::  neig                            !  number of eigenvectors required (maybe modified on output)
  real(REAL64)                       ::  rkpt(3)                         !  j-th component in lattice coordinates of the k-point

  integer                            ::  irk
  integer                            ::  ipr,nrk2
  logical                            ::  lfile
  integer                            ::  nbandi2,nx2,ny2,nz2
  real(REAL64)                       ::  sx2,sy2,sz2
  real(REAL64)                       ::  vcell, bdot(3,3)

  integer                            ::  io                              !  tape numbers
  character(len=128)                 ::  syscom

  logical :: lscl, lso
  integer :: identif
  character(len=12)                  ::  filedos

!      counters
  integer                            ::  iband, jband


  character(len=50)                  ::  title                           !<  title for plots
  character(len=140)                 ::  subtitle                        !<  subtitle for plots


  io = 12

  ipr = 2

  write(title,*) "title"
  write(subtitle,*) "subtitle"

  lfile = .false.
  open(unit=11,file='DOS_MESH.DAT',status='old',err=10, form = 'formatted')
  lfile = .true.
  10    continue

  if(lfile) then
    read(11,*) nbandi2,nx2,ny2,nz2,sx2,sy2,sz2
    close(unit=11)
    mxdbnd = nbandi2
  else
    mxdbnd = nint(ztot) + 4
    nbandi2 = mxdbnd
    nx2 = 8
    ny2 = 8
    nz2 = 8
    sx2 = 0.0
    sy2 = 0.0
    sz2 = 0.0
  endif

  mxdpnt = nx2*ny2*nz2

  allocate(kmap2(3,nx2,ny2,nz2))
  allocate(nband2(mxdpnt))
  allocate(indk2(6,mxdpnt))
  allocate(rk2(3,mxdpnt))
  allocate(w2(mxdpnt))

  call int_pnt(nbandi2, nx2, ny2, nz2, sx2, sy2, sz2, ipr,               &
      adot,                                                              &
      ntrans, mtrx,                                                      &
      nrk2, rk2, w2, nband2, indk2, kmap2,                               &
      mxdpnt, mxdbnd)

  neig = nbandi2

  !      writes in a file for later processing


  call adot_to_bdot(adot,vcell,bdot)

  allocate(ev_interp(noiData%nband))

  !      loop over k-points

  if(noiData%lso==1) then
    allocate(e_of_k_so(2*neig,nrk2))
    allocate(e_of_k(1,1))
    allocate(basxpsi_of_k(noiData%nband,noiData%nband,nrk2))
    lscl = .false.
    lso  = .true.
  else
    allocate(e_of_k_so(1,1))
    allocate(e_of_k(neig,nrk2))
    allocate(basxpsi_of_k(noiData%nband,noiData%nband,nrk2))
    lscl = .true.
    lso  = .false.
  endif

  identif = 0

  filedos = 'dos_file.dat'

  allocate (psi(noiData%nband,noiData%nband))

  do irk=1,nrk2

    rkpt(1) = rk2(1,irk)
    rkpt(2) = rk2(2,irk)
    rkpt(3) = rk2(3,irk)

    call fi_hamiltonian_get_hk(noiData%fiData, rkpt, noiData%Hao_tr, noiData%nband)
    call fi_hamiltonian_get_sk(noiData%fiData, rkpt, noiData%Uao, noiData%nband)

    call ao_int_DiagByLowdin(noiData%nband, noiData%Hao_tr, noiData%Uao, ev_interp, psi)

!    call NonOrthoInterpRun(noiData,rkpt,ev_interp)

    if(noiData%lso==1) then
      do iband = 1, 2*neig
        e_of_k_so(iband,irk) = ev_interp(iband)
      do jband = 1, noiData%nband
        basxpsi_of_k(jband,iband,irk) = real(psi(jband,iband)*conjg(psi(jband,iband)),REAL64)
      enddo
      enddo
    else
      do iband = 1, neig
        e_of_k(iband,irk) = ev_interp(iband)
      do jband = 1, noiData%nband
        basxpsi_of_k(jband,iband,irk) = real(psi(jband,iband)*conjg(psi(jband,iband)),REAL64)
      enddo
      enddo
    endif

    call progress(irk,nrk2)
    call flush(6)

  enddo

  open(unit=190,file="dos_file_proj.dat", form="unformatted")
  write(190) nrk2
  write(190) neig
  write(190) noiData%nband
  do irk=1,nrk2
    write(190) basxpsi_of_k(:,:,irk)
  enddo

  close(unit=190)

  call out_dos_write(filedos, io, title, subtitle,                       &
      lscl, lso, identif,                                                &
      nrk2, nx2, ny2, nz2, ztot, adot, ntrans, mtrx,                     &
      nband2, rk2, w2, indk2, kmap2, e_of_k, e_of_k_so,                  &
      mxdpnt, mxdbnd)


  if(noiData%lso==1) then
    syscom = "rm PW_DOS.DAT"
  else
    syscom = "rm PW_DOS_SO.DAT"
  endif
  call execute_command_line(trim(syscom))

  !      Cleanup

  deallocate(basxpsi_of_k)

  deallocate(e_of_k_so)
  deallocate(e_of_k)

  deallocate(kmap2)
  deallocate(nband2)
  deallocate(indk2)
  deallocate(rk2)
  deallocate(w2)

  deallocate(ev_interp)
  deallocate(psi)

  return

end subroutine ao_interpolation_out_dos








