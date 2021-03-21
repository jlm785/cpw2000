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

!>  Expands the optical data from the irreducible grid to the full grid
!>  Uses files to avoid exceeding RAM.

subroutine opt_grid_unsym(el, filedhdrk, io_dhdrk, fileunsym, io_unsym,      &
  kmap, mtrx, neig, nx,ny,nz, nrk, mxdbnd)

! Extracted from previous code. 23 October 2020. JLM
! Modified to use files to avoid exceeding RAM.  11 December 2020. JLM
! copyright  Carlos Loia Reis/INESC-MN

! version 4.99

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)
  integer, parameter          :: REAL32 = selected_real_kind(6)

! input

  integer, intent(in)                ::  mxdbnd                          !<  array dimension for the number of bands
  integer, intent(in)                ::  nrk                             !<  number of irreducible k-points
  integer, intent(in)                ::  nx,ny,nz                        !<  grid size

  real(REAL64), intent(in)           ::  el(mxdbnd,nrk)                  !<  eigenvalues in Hartree

  character(len=16), intent(in)      ::  filedhdrk                       !<  file with ev, dhdrk, symmetric
  character(len=16), intent(in)      ::  fileunsym                       !<  file with ev, dhdrk, unsymmetrized

  integer, intent(in)                ::  io_dhdrk                        !<  tape number for symmetrized results
  integer, intent(in)                ::  io_unsym                        !<  tape number for unsymmetrized results


  integer, intent(in)                ::  kmap(3,nx,ny,nz)                !<  kmap(1,...) corresponding k-point, kmap(2,...) = 1 additional inversion, kmap(3,...) symmetry operation
  integer , intent(in)               ::  mtrx(3,3,48)                    !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group

  integer, intent(in)                ::  neig                            !<  number of eigenvectors

! allocatable arrays

  complex(REAL32), allocatable       ::  dhdrk_32(:,:,:)                 !  d h / d rk

  real(REAL32), allocatable          ::  ev_grid_32(:)                   !  eigenvalue in grid
  complex(REAL32), allocatable       ::  dhdrk_grid_32(:,:,:)            !  d h / d rk in grid


! local variables

  integer                      ::  irec_err, irk_rd, ir_size
  integer                      ::  nxyz

! counters

  integer                      ::  i, j, k, iorb, kk, jj


  allocate(dhdrk_32(neig,neig,3))

  allocate(ev_grid_32(neig))
  allocate(dhdrk_grid_32(neig,neig,3))

  inquire(iolength = ir_size) irk_rd, dhdrk_32(:,:,:)
 
  open(unit = io_dhdrk, file = trim(filedhdrk), access="direct", recl=ir_size)

  inquire(iolength = ir_size) ev_grid_32(:), dhdrk_grid_32(:,:,:)
 
  open(unit = io_unsym, file = trim(fileunsym), access="direct", recl=ir_size)

  do k = 1,nz
  do j = 1,ny
  do i = 1,nx
    kk = kmap(1,i,j,k)
    jj = kmap(3,i,j,k)
    nxyz = (k-1)*ny*nx + (j-1)*nx + i 


    read(io_dhdrk, rec=kk, iostat=irec_err) irk_rd, dhdrk_32(:,:,:)

    do iorb=1,neig
      ev_grid_32(iorb) = el(iorb,kk)
    enddo

    dhdrk_grid_32(:,:,1) = mtrx(1,1,jj)*dhdrk_32(:,:,1) +                &    
                           mtrx(2,1,jj)*dhdrk_32(:,:,2) +                &
                           mtrx(3,1,jj)*dhdrk_32(:,:,3)
    dhdrk_grid_32(:,:,2) = mtrx(1,2,jj)*dhdrk_32(:,:,1) +                &
                           mtrx(2,2,jj)*dhdrk_32(:,:,2) +                &
                           mtrx(3,2,jj)*dhdrk_32(:,:,3)
    dhdrk_grid_32(:,:,3) = mtrx(1,3,jj)*dhdrk_32(:,:,1) +                &
                           mtrx(2,3,jj)*dhdrk_32(:,:,2) +                &
                           mtrx(3,3,jj)*dhdrk_32(:,:,3)

    write(unit = io_unsym, rec = nxyz, iostat=irec_err) ev_grid_32(:), dhdrk_grid_32(:,:,:)

  enddo
  enddo
  enddo

  deallocate(dhdrk_32)

  deallocate(ev_grid_32)
  deallocate(dhdrk_grid_32)
 
  close(unit = io_dhdrk)
  close(unit = io_unsym)

  return
end subroutine opt_grid_unsym
 
