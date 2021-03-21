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

!>  Fills the e (eigenvalue difference) and fk (dipole) in grid


subroutine opt_set_in_grid(ix,iy, adot, fileunsym, io_unsym, filetmp, io_tmp, filegrid, io_grid, &
                           neig, nval, ncond, nx, ny, nz)

! Written by Carlos Loia Reis. July 2020
! Modified, documentation, 20 September 2020. JLM
! Modified, merged two subroutines. 20 October 2020. JLM
! Modified, use of files to avoid exceeding RAM. 11 december 2020. JLM
! copyright  Carlos Loia Reis/INESC-MN

! version 4.99

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)
  integer, parameter          :: REAL32 = selected_real_kind(6)

! input
  
  integer, intent(in)                ::  ix, iy                          !<  tensor component(s) to 
  
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in real space
  
  character(len=16), intent(in)      ::  fileunsym                       !<  file with ev, dhdrk, unsymmetrized
  integer, intent(in)                ::  io_unsym                        !<  tape number for unsymmetrized results
  character(len=16), intent(in)      ::  filetmp                         !<  temporary file
  integer, intent(in)                ::  io_tmp                          !<  tape number for temp file
  character(len=16), intent(in)      ::  filegrid                        !<  file with grid data
  integer, intent(in)                ::  io_grid                         !<  tape number for grid

  integer, intent(in)                ::  neig                            !<  number of eigenvectors
  integer, intent(in)                ::  nval                            !<  number of valence bands
  integer, intent(in)                ::  ncond                           !<  number of conduction

  integer, intent(in)                ::  nx,ny,nz                        !<  grid size
 
! local allocatable arrays
    
  integer, allocatable               ::  indx(:)

  real(REAL64), allocatable          ::  efmei(:)
  real(REAL32), allocatable          ::  efmei_32(:)
  complex(REAL32), allocatable       ::  ematif_32(:)
  
  real(REAL32), allocatable          ::  ev_grid_32(:)                   !  eigenvalue in one k-point
  complex(REAL32), allocatable       ::  dhdrk_grid_32(:,:,:)            !  d h / d rk in one k-point

  real(REAL32), allocatable          ::  egrid_exc_32(:)                 !  excitation values in one k-point
  real(REAL32), allocatable          ::  fkgrid_exc_32(:)                !  dipole matrix elements in one k-point

  real(REAL32), allocatable          ::  egrid_part(:,:,:,:)             !  excitation values for several k-points
  real(REAL32), allocatable          ::  fkgrid_part(:,:,:,:)            !  dipole matrix elements for several k-points

  real(REAL32), allocatable          ::  egrid_32(:,:,:)                 !  excitation values in one "joint band"
  real(REAL32), allocatable          ::  fkgrid_32(:,:,:)                !  dipole matrix elements in one "joint band"

! local variables
  
  real(REAL64)                 ::  avec(3,3), bvec(3,3)
  real(REAL32)                 ::  adot_32(3,3), avec_32(3,3)

  integer                      ::  irec_err, ir_size
  integer                      ::  nxyz

! constants

  real(REAL64), parameter      ::  PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter      ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter   ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)

! counters

  integer                      ::  i, j, k, ij, iorb, jorb, n, m

! begin
  
  allocate(efmei(nval*ncond))
  allocate(efmei_32(nval*ncond))
  allocate(ematif_32(nval*ncond))

  allocate(indx(nval*ncond))
  
  call adot_to_avec_sym(adot,avec,bvec)

  adot_32(:,:) = adot(:,:)
  avec_32(:,:) = avec(:,:)

  allocate(ev_grid_32(neig))
  allocate(dhdrk_grid_32(neig,neig,3))

  allocate(egrid_exc_32(nval*ncond))
  allocate(fkgrid_exc_32(nval*ncond))

  inquire(iolength = ir_size) ev_grid_32(:), dhdrk_grid_32(:,:,:)
 
  open(unit = io_unsym, file = trim(fileunsym), access="direct", recl=ir_size)
  
  inquire(iolength = ir_size) egrid_exc_32(:), fkgrid_exc_32(:)
 
  open(unit = io_tmp, file = trim(filetmp), access="direct", recl=ir_size)
  
  
  do k = 1,nz
  do j = 1,ny
  do i = 1,nx
   nxyz = (k-1)*ny*nx + (j-1)*nx + i 

   read(unit = io_unsym, rec = nxyz, iostat=irec_err) ev_grid_32(:), dhdrk_grid_32(:,:,:)
    
    do iorb = 1,nval
    do jorb = 1,ncond
      ij = (iorb-1)*ncond + jorb
      efmei_32(ij) = ev_grid_32(nval+jorb) - ev_grid_32(iorb)
      ematif_32(ij) = C_ZERO

      if(ix == 0 .and. iy == 0) then

        do n = 1,3
        do m = 1,3
          ematif_32(ij) = ematif_32(ij) + dhdrk_grid_32(nval+jorb,iorb,n)   &
                         * adot_32(n,m) * dhdrk_grid_32(iorb,nval+jorb,m)

       enddo
       enddo
                   
       ematif_32(ij) = ematif_32(ij) / (4*PI*PI) /(3*UM)
       
      else

        ematif_32(ij) = ematif_32(ij) +                                  &
        (avec_32(ix,1)*dhdrk_grid_32(nval+jorb,iorb,1)  +                &
         avec_32(ix,2)*dhdrk_grid_32(nval+jorb,iorb,2)  +                &
         avec_32(ix,3)*dhdrk_grid_32(nval+jorb,iorb,3)) *                &
        (avec_32(iy,1)*dhdrk_grid_32(iorb,nval+jorb,1)  +                &
         avec_32(iy,2)*dhdrk_grid_32(iorb,nval+jorb,2)  +                &
         avec_32(iy,3)*dhdrk_grid_32(iorb,nval+jorb,3)) 
                          
         ematif_32(ij) = ematif_32(ij) / (4*PI*PI)

      endif

    enddo
    enddo

    indx(:) = 0
    efmei(:) = efmei_32(:)
    
    call sort(nval*ncond,efmei,indx)

    do ij = 1,nval*ncond
      egrid_exc_32(ij)  = efmei_32(indx(ij))
      fkgrid_exc_32(ij) = real(ematif_32(indx(ij)))
    enddo

    write(unit = io_tmp, rec = nxyz, iostat=irec_err) egrid_exc_32(:), fkgrid_exc_32(:)
  
  enddo
  enddo
  enddo
   
  close(unit = io_unsym)

  deallocate(ev_grid_32)
  deallocate(dhdrk_grid_32)

  deallocate(indx)
  deallocate(efmei)
  deallocate(efmei_32)
  deallocate(ematif_32)
 
  allocate(egrid_part(nx,ny,nz,ncond))
  allocate(fkgrid_part(nx,ny,nz,ncond))

  allocate(egrid_32(nx,ny,nz))
  allocate(fkgrid_32(nx,ny,nz))

  inquire(iolength = ir_size) egrid_32(:,:,:), fkgrid_32(:,:,:)
 
  open(unit = io_grid, file = trim(filegrid), access="direct", recl=ir_size)

  do iorb = 1,nval

    do k = 1,nz
    do j = 1,ny
    do i = 1,nx
      nxyz = (k-1)*ny*nx + (j-1)*nx + i 

      read(unit = io_tmp, rec = nxyz, iostat=irec_err) egrid_exc_32(:), fkgrid_exc_32(:)

      do jorb = 1,ncond
        ij = (iorb-1)*ncond + jorb
        egrid_part(i,j,k,jorb) = egrid_exc_32(ij)
        fkgrid_part(i,j,k,jorb) = fkgrid_exc_32(ij)
      enddo

    enddo
    enddo
    enddo
    
    do jorb = 1,ncond
      ij = (iorb-1)*ncond + jorb

      egrid_32(:,:,:) = egrid_part(:,:,:,jorb)
      fkgrid_32(:,:,:) = fkgrid_part(:,:,:,jorb)

      write(unit = io_grid, rec = ij, iostat=irec_err) egrid_32(:,:,:), fkgrid_32(:,:,:)

    enddo
  enddo
  
  close(unit = io_grid)
  
  close(unit = io_tmp)
 
  deallocate(egrid_part)
  deallocate(fkgrid_part)

  deallocate(egrid_exc_32)
  deallocate(fkgrid_exc_32)

  deallocate(egrid_32)
  deallocate(fkgrid_32)
  
  return
end subroutine opt_set_in_grid
