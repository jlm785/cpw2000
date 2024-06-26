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

subroutine ao_interpolation_write_chk(ortho, nband)
  use NonOrthoInterp  

  implicit none

  integer, parameter            :: REAL64 = selected_real_kind(12)

  type (fiData_t)               :: ortho 
  integer                       :: nband
  
!  real(REAL64), intent(in)      :: adot(3,3)
  
  real(REAL64)                  :: avec(3,3), bvec(3,3)
    
  integer                         :: i,j,k,l, ikpt


  integer, parameter              :: chk_unit = 17
  
  character(len=33)               :: header
!  integer                         :: num_bands
!  integer                         :: num_exclude_bands
!  integer, allocatable            :: exclude_bands(:)
  
!  real(REAL64)                    :: real_lattice(3,3)
!  real(REAL64)                    :: recip_lattice(3,3)
  
  integer                         :: num_kpts
!  integer                         :: mp_grid(3)
      
!  real(REAL64),allocatable        :: kpt_latt(:,:)
  
  integer                         :: nntot
  integer                         :: num_wann
  
  character(len=20)               :: chkpt1
!  logical                         :: have_disentangled
  
!  real(REAL64)                    :: omega_invariant
  
!  logical,  allocatable           :: lwindow(:, :)
!  integer,  allocatable           :: ndimwin(:)
  
!  complex(REAL64), allocatable    :: u_matrix_opt(:,:,:)
  
  complex(REAL64), allocatable    :: u_matrix(:,:,:)
  complex(REAL64), allocatable    :: m_matrix(:,:,:,:)

  real(REAL64),  allocatable      :: wannier_centres(:, :)
  real(REAL64),   allocatable     :: wannier_spreads(:)
  
  real(REAL64), parameter         :: bohr_to_ang= 0.52917720859
  
  integer, allocatable            :: nnlist(:, :)                          !*
  integer, parameter              :: num_nnmax = 12
  
      
  num_wann =  nband
  num_kpts =  ortho%num_kpts
  nntot  = 8                                                             ! this must be calculated for a given lattice (to steal ..) 
            
  open(unit=chk_unit, file="wan.chk", form="unformatted")
  
  call adot_to_avec_sym(ortho%real_metric,avec,bvec)
  
  
  write(header,*) " auto generated by mtb"
    
  write(chk_unit) header
  write(chk_unit) nband
  write(chk_unit) 0                                                      ! num_exclude_bands

  write(chk_unit)                                                        ! num_exclude_bands
  
  write (chk_unit) ((avec(j, i)*bohr_to_ang, i=1, 3), j=1, 3)                        ! Real lattice
  write (chk_unit) ((bvec(j, i)/bohr_to_ang, i=1, 3), j=1, 3)                        ! Reciprocal lattice
  write (chk_unit)  num_kpts                                             ! Number of k-points
  write (chk_unit) (ortho%mp_grid(i), i=1, 3)                            ! M-P grid
  
  write (chk_unit) ((ortho%kpt_latt(i, k), i=1, 3), k=1,  num_kpts)      ! K-points  
  write (chk_unit) nntot                                                     ! Number of nearest k-point neighbours nntot
  write (chk_unit) nband                                                 ! Number of wannier functions
  


  write (888,*) ((avec(j, i)*bohr_to_ang, i=1, 3), j=1, 3)                        ! Real lattice
  write (888,*) ((bvec(j, i)/bohr_to_ang, i=1, 3), j=1, 3)                        ! Reciprocal lattice
  write (888,*)  num_kpts                                             ! Number of k-points
  write (888,*) (ortho%mp_grid(i), i=1, 3)                            ! M-P grid
  
  write (888,*) ((ortho%kpt_latt(i, k), i=1, 3), k=1,  num_kpts)      ! K-points  
  write (888,*) nntot                                                     ! Number of nearest k-point neighbours nntot
  write (888,*) nband                                                 ! Number of wannier functions



  
  write(chkpt1,*)  "postwann"
  
  write (chk_unit) chkpt1                                                ! Position of checkpoint
  write (chk_unit) .false.                                               ! Whether a disentanglement has been performed
  
  
  allocate(u_matrix(num_wann,num_wann,num_kpts))
  allocate(m_matrix(num_wann,num_wann,nntot,num_kpts))

  allocate(wannier_centres(3,num_wann))
  allocate(wannier_spreads(num_wann))
 
  open(unit=253,file="amn.dat",  form="unformatted")  

  open(unit = 99, file = "mmnk.dat", form="unformatted")

  allocate(nnlist(num_kpts,num_nnmax))
  
  read(99) nntot
  read(99) nnlist
    
  write(*,*) 'nntot is', nntot
     
   do ikpt = 1, num_kpts
      write(*,'(12i5)') ikpt,( nnlist(ikpt, i), i=1,nntot)
   enddo
  
!  stop "stoped in write_chk debug purposes"
    
  do k=1, num_kpts
!  do i=1, num_wann
!  do j=1, num_wann
!    u_matrix(i,j,k) = cmplx(1.0,0.0,REAL64)                               ! dont forget to change this one
    read(253) u_matrix(:,:,k)
    do i=1,nntot
!      read(99) m_matrix(:,:,i,k)    
    enddo
    
    
!  enddo
!  enddo
  enddo
  
  close(253)
  close(99)

  do l=1,nntot
  do k=1, num_kpts
  do i=1, num_wann
  do j=1, num_wann
!    m_matrix(i,j,l,k) = cmplx(1.0,0.0,REAL64)                             ! this one is way more tricky..
  enddo
  enddo
  enddo
  enddo
  
  do j=1, num_wann
    wannier_centres(1,j) = 0.0_REAL64
    wannier_centres(2,j) = 0.0_REAL64
    wannier_centres(3,j) = 0.0_REAL64
    
    wannier_spreads(j) = 0.0_REAL64
  
  enddo

    
  write (chk_unit) (((u_matrix(j, i, k), i=1, num_wann), j=1, num_wann), k=1, num_kpts)                  ! U_matrix    


  write (chk_unit) ((((m_matrix(i, j, k, l), i=1, num_wann), j=1, num_wann), k=1, nntot), l=1, num_kpts) ! M_matrix

  write (chk_unit) ((wannier_centres(i, j), i=1, 3), j=1, num_wann)
  write (chk_unit) (wannier_spreads(i), i=1, num_wann)
  
  
  deallocate(u_matrix)
  deallocate(m_matrix)

  deallocate(wannier_centres)
  deallocate(wannier_spreads)
  
  
  
  close(unit=chk_unit)


end subroutine
