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

module FourierInterpolation
! Non Orthogonal Fourier Interpolation package.
! adapted from Wannier90 w90_hamiltonian
! Copyright Carlos Loia Reis
! June 2014.

  implicit none
  integer, parameter,private          :: REAL64 = selected_real_kind(12)

  type fiData_t
    complex(REAL64), allocatable  :: ham_k(:,:,:)
    complex(REAL64), allocatable  :: S_k(:,:,:)
    complex(REAL64), allocatable  :: ham_r(:,:,:)
    complex(REAL64), allocatable  :: S_r(:,:,:)

    integer                       :: mp_grid(3)
    integer                       :: ws_search_size(3)
    real(REAL64)                  :: real_metric(3,3)

    integer,  allocatable         :: irvec(:,:)
    integer,  allocatable         :: ndegen(:)


    integer                       :: nrpts
    integer                       :: num_kpts
    real(REAL64),    allocatable  :: kpt_latt(:,:)
        
    complex(REAL64), allocatable  ::  dh0drk_GridK(:,:,:,:)     ! d <Psi|H|Psi> d k
    complex(REAL64), allocatable  ::  dh0drk_GridR(:,:,:,:)     ! d <Psi|H|Psi> d k
    
    
  end type
  
!  integer, parameter, public          :: dp = selected_real_kind(15,300)
  
  ! Module variables
!  logical, save :: ham_have_setup=.false.
!  logical, save :: have_ham_r=.false.
!  logical, save :: have_ham_k=.false.
!  logical, save :: have_dh0drk=.false.

contains

  subroutine fi_hamiltonian_get_hk(this,rkpt,ham, nband)
    integer, parameter          :: REAL64 = selected_real_kind(12)        
    type(fiData_t):: this
    real(REAL64), parameter     :: pi=3.141592653589793238462643383279_REAL64
    real(REAL64), parameter     :: twopi=2.0_REAL64*pi
    complex(REAL64), parameter  :: cmplx_i = (0.0_REAL64,1.0_REAL64)
    complex(REAL64), parameter  :: cmplx_0 = (0.0_REAL64,0.0_REAL64)

    integer nband
    real(REAL64) rkpt(3)    
    complex(REAL64) ham(nband,nband)

    real(REAL64) rdotk
    complex(REAL64) fac
    
    integer loop_rpt  
!    integer loop_kpt 
    
    ham=cmplx_0
    do loop_rpt=1,this%nrpts
        rdotk=twopi*dot_product(rkpt,this%irvec(:,loop_rpt))
        fac=exp(cmplx_i*rdotk)/real(this%ndegen(loop_rpt),REAL64)
        ham=ham+fac*this%ham_r(:,:,loop_rpt)
    end do
  end subroutine fi_hamiltonian_get_hk


  subroutine fi_hamiltonian_get_sk(this,rkpt,S, nband)  
    implicit none
    integer, parameter          :: REAL64 = selected_real_kind(12)        
    type(fiData_t):: this
    real(REAL64), parameter     :: pi=3.141592653589793238462643383279_REAL64
    real(REAL64), parameter     :: twopi=2.0_REAL64*pi
    complex(REAL64), parameter  :: cmplx_i = (0.0_REAL64,1.0_REAL64)
    complex(REAL64), parameter  :: cmplx_0 = (0.0_REAL64,0.0_REAL64)

    integer nband
    real(REAL64) rkpt(3)    
    complex(REAL64) S(nband,nband)

    real(REAL64) rdotk
    complex(REAL64) fac
    
    integer loop_rpt  
!    integer loop_kpt
    
    S=cmplx_0
    do loop_rpt=1,this%nrpts
        rdotk=twopi*dot_product(rkpt,this%irvec(:,loop_rpt))
        fac=exp(cmplx_i*rdotk)/real(this%ndegen(loop_rpt),REAL64)
        S=S+fac*this%S_r(:,:,loop_rpt)
    end do
  end subroutine fi_hamiltonian_get_sk
          
  subroutine fi_hamiltonian_get_hr(this)
    implicit none
    integer, parameter          :: REAL64 = selected_real_kind(12)        
    type(fiData_t):: this
    real(REAL64), parameter     :: pi=3.141592653589793238462643383279_REAL64
    real(REAL64), parameter     :: twopi=2.0_REAL64*pi
    complex(REAL64), parameter  :: cmplx_i = (0.0_REAL64,1.0_REAL64)
    complex(REAL64), parameter  :: cmplx_0 = (0.0_REAL64,0.0_REAL64)

    integer loop_kpt,loop_rpt
    real(REAL64) rdotk
    complex(REAL64) fac
    
    this%ham_r = cmplx_0
    this%S_r = cmplx_0
    
    do loop_rpt=1,this%nrpts
      do loop_kpt=1,this%num_kpts
          rdotk=twopi*dot_product(this%kpt_latt(:,loop_kpt),real(this%irvec(:,loop_rpt),REAL64))
          fac=exp(-cmplx_i*rdotk)/real(this%num_kpts,REAL64)
          this%ham_r(:,:,loop_rpt)=this%ham_r(:,:,loop_rpt)+fac*this%ham_k(:,:,loop_kpt)
          this%S_r(:,:,loop_rpt)=this%S_r(:,:,loop_rpt)+fac*this%S_k(:,:,loop_kpt)
        enddo
      call progress(loop_rpt,this%nrpts)
      call flush(6)        
    enddo

  end subroutine fi_hamiltonian_get_hr


  subroutine fi_hamiltonian_get_Op_k(this,rkpt,Op_R,Op_k, nband)
    implicit none
    integer, parameter          :: REAL64 = selected_real_kind(12)        
    type(fiData_t):: this
    real(REAL64), parameter     :: pi=3.141592653589793238462643383279_REAL64
    real(REAL64), parameter     :: twopi=2.0_REAL64*pi
    complex(REAL64), parameter  :: cmplx_i = (0.0_REAL64,1.0_REAL64)
    complex(REAL64), parameter  :: cmplx_0 = (0.0_REAL64,0.0_REAL64)
    integer nband
    real(REAL64) rkpt(3)    

    complex(REAL64) Op_k(nband,nband,3)
    complex(REAL64) Op_R(nband,nband,3,this%nrpts)

    real(REAL64) rdotk
    complex(REAL64) fac
    
    integer loop_rpt  
        
    Op_k(:,:,:)=cmplx_0
        
    do loop_rpt=1,this%nrpts
        rdotk=twopi*dot_product(rkpt,this%irvec(:,loop_rpt))
        fac=exp(cmplx_i*rdotk)/real(this%ndegen(loop_rpt),REAL64)
        Op_k(:,:,1)=Op_k(:,:,1)+fac*Op_r(:,:,1,loop_rpt)
        Op_k(:,:,2)=Op_k(:,:,2)+fac*Op_r(:,:,2,loop_rpt)
        Op_k(:,:,3)=Op_k(:,:,3)+fac*Op_r(:,:,3,loop_rpt)
    end do
    
  end subroutine fi_hamiltonian_get_Op_k

  subroutine fi_hamiltonian_set_in_Rgrid(this,Op_k,Op_R,nband)
    implicit none
    integer, parameter          :: REAL64 = selected_real_kind(12)    
    type(fiData_t):: this  
    real(REAL64), parameter     :: pi=3.141592653589793238462643383279_REAL64
    real(REAL64), parameter     :: twopi=2.0_REAL64*pi
    complex(REAL64), parameter  :: cmplx_i = (0.0_REAL64,1.0_REAL64)
    complex(REAL64), parameter  :: cmplx_0 = (0.0_REAL64,0.0_REAL64)

    integer loop_kpt,loop_rpt
    real(REAL64) rdotk
    complex(REAL64) fac
    
    integer nband
    
    complex(REAL64) Op_k(nband,nband,3,this%num_kpts)
    complex(REAL64) Op_R(nband,nband,3,this%nrpts)

    Op_R(:,:,:,:) = cmplx_0
    
    do loop_rpt=1,this%nrpts
      do loop_kpt=1,this%num_kpts
          rdotk=twopi*dot_product(this%kpt_latt(:,loop_kpt),real(this%irvec(:,loop_rpt),REAL64))
          fac=exp(-cmplx_i*rdotk)/real(this%num_kpts,REAL64)
          Op_R(:,:,1,loop_rpt)=Op_R(:,:,1,loop_rpt)+fac*Op_k(:,:,1,loop_kpt)
          Op_R(:,:,2,loop_rpt)=Op_R(:,:,2,loop_rpt)+fac*Op_k(:,:,2,loop_kpt)
          Op_R(:,:,3,loop_rpt)=Op_R(:,:,3,loop_rpt)+fac*Op_k(:,:,3,loop_kpt)
        enddo
        
!        write(6,FMT="(3x,A1,A,t21,F6.2,A)",ADVANCE="NO") achar(13), &
!        & " Percent Complete: ", (real(loop_rpt)/real(this%nrpts))*100.0, "%"
      call progress(loop_rpt,this%nrpts)
      call flush(6)
        
    enddo

  end subroutine fi_hamiltonian_set_in_Rgrid


  subroutine fi_hamiltonian_setup(this,nband,n1,n2,n3,ws_n1,ws_n2,ws_n3,adot)
    implicit none
    type(fiData_t):: this
    
    integer :: ierr

    integer nband    
    integer n1,n2,n3
    integer ws_n1, ws_n2,ws_n3
    integer i,j,k,ikpt
    real(REAL64) adot(3,3)

    complex(REAL64), parameter  :: cmplx_0 = (0.0_REAL64,0.0_REAL64)
        
    this%num_kpts = n1*n2*n3
    allocate(this%kpt_latt(3,this%num_kpts))
    
    this%mp_grid(1) = n1
    this%mp_grid(2) = n2
    this%mp_grid(3) = n3
    
    this%real_metric(:,:) = adot(:,:)

    write(*,'("in setup nband, nk_grid, ws_search_size", 7i5)') nband,n1,n2,n3, ws_n1, ws_n2, ws_n3
    
    this%ws_search_size(1) = ws_n1
    this%ws_search_size(2) = ws_n2
    this%ws_search_size(3) = ws_n3
    

!    if (ham_have_setup) return

    !
    !
    ! Set up Wigner-Seitz vectors
    !
    call fi_hamiltonian_wigner_seitz(this,count_pts=.true.,iprint=3)
    !
    allocate(this%irvec(3,this%nrpts),stat=ierr)
    this%irvec=0
    !
    allocate(this%ndegen(this%nrpts),stat=ierr)
    this%ndegen=0
    !
    allocate(this%ham_r(nband,nband,this%nrpts),stat=ierr)
    allocate(this%S_r(nband,nband,this%nrpts),stat=ierr)
!    if (ierr/=0) call io_error('Error in allocating ham_r in hamiltonian_setup')
    this%ham_r(:,:,:)= cmplx_0
    this%S_r(:,:,:)  = cmplx_0
    !
    allocate(this%ham_k(nband,nband,this%num_kpts),stat=ierr)
    allocate(this%S_k(nband,nband,this%num_kpts),stat=ierr)
!    if (ierr/=0) call io_error('Error in allocating ham_k in hamiltonian_setup')
    this%ham_k(:,:,:)= cmplx_0
    this%S_k(:,:,:)  = cmplx_0
    !
    ! Set up the wigner_seitz vectors
    !
    call fi_hamiltonian_wigner_seitz(this,count_pts=.false.,iprint=3)
    !

!    ham_have_setup = .true.
    
    ikpt=1
    do i=1,n1
    do j=1,n2
    do k=1,n3
      this%kpt_latt(1,ikpt) = (i-1)*1.0D0/n1
      this%kpt_latt(2,ikpt) = (j-1)*1.0D0/n2
      this%kpt_latt(3,ikpt) = (k-1)*1.0D0/n3    
      ikpt=ikpt+1
    enddo
    enddo
    enddo
    
!    write(*,*) 'kpt lattice'
!    do ikpt=1,this%num_kpts
!      write(*,'(3f8.5)') this%kpt_latt(:,ikpt)
!    enddo
        
    return
  end subroutine fi_hamiltonian_setup
  
  
  subroutine fi_hamiltonian_wigner_seitz_uniform(this,count_pts, iprint)
    implicit none
    
    ! uniform conjugate grid of r-points irvec. CLR

    type(fiData_t):: this
    
    integer i, n1, n2, n3, iprint
    logical count_pts
    
    integer, parameter :: stdout = 6
    
    
    if (count_pts) then 
      this %nrpts = this%mp_grid(1)*this%mp_grid(2)*this%mp_grid(3)
      return    
    endif
        
    i=1    
    do n1=-this%mp_grid(1)/2 + 1, this%mp_grid(1)/2 
    do n2=-this%mp_grid(2)/2 + 1, this%mp_grid(2)/2 
    do n3=-this%mp_grid(3)/2 + 1, this%mp_grid(3)/2     
          this%irvec(1, i) = n1
          this%irvec(2, i) = n2
          this%irvec(3, i) = n3
          this%ndegen(i) = 1.0_REAL64
          i=i+1
    enddo
    enddo
    enddo
    
    write(*,*) 'i=',i

    if (iprint >= 3 ) then
      write (stdout, '(1x,i4,a,/)') this%nrpts, ' lattice points in Wigner-Seitz supercell:'
      do i = 1, this%nrpts
        write (stdout, '(4x,a,3(i3,1x),a,i2)') '  vector ', this%irvec(1, i), this%irvec(2, i), &
          this%irvec(3, i), '  degeneracy: ', this%ndegen(i)
      enddo
      
    endif



  end subroutine fi_hamiltonian_wigner_seitz_uniform
  
  !================================================================================!
  subroutine fi_hamiltonian_wigner_seitz(this,count_pts, iprint)
    implicit none

    ! <stolen> from wannier90 3.0.0 and adapted to our needs CLR
    !================================================================================!
    !! Calculates a grid of points that fall inside of (and eventually on the
    !! surface of) the Wigner-Seitz supercell centered on the origin of the B
    !! lattice with primitive translations nmonkh(1)*a_1+nmonkh(2)*a_2+nmonkh(3)*a_3
    !================================================================================!

    type(fiData_t):: this
    
    integer, parameter :: stdout = 6

    real(REAL64), parameter:: eps7=1.0D-7
    real(REAL64), parameter:: eps8=1.0D-8

    real(REAL64), parameter:: ws_distance_tol=1.0D-5

    ! irvec(i,irpt)     The irpt-th Wigner-Seitz grid point has components
    !                   irvec(1:3,irpt) in the basis of the lattice vectors
    ! ndegen(irpt)      Weight of the irpt-th point is 1/ndegen(irpt)
    ! nrpts             number of Wigner-Seitz grid points


    logical, intent(in) :: count_pts
    !! Only count points and return
    integer       :: ndiff(3)
    real(REAL64) :: tot, dist_min
    real(REAL64), allocatable :: dist(:)
    integer       ::iprint, n1, n2, n3, i1, i2, i3, icnt, i, j, ierr, dist_dim


    dist_dim = 1
    do i = 1, 3
      dist_dim = dist_dim*((this%ws_search_size(i) + 1)*2 + 1)
    end do
    
    allocate (dist(dist_dim), stat=ierr)

    ! The Wannier functions live in a supercell of the real space unit cell
    ! this supercell is mp_grid unit cells long in each direction
    !
    ! We loop over grid points r on a unit cell that is (2*ws_search_size+1)**3 times
    ! larger than this primitive supercell.
    !
    ! One of these points is in the W-S cell if it is closer to R=0 than any of the
    ! other points, R (where R are the translation vectors of the supercell)

    ! In the end nrpts contains the total number of grid
    ! points that have been found in the Wigner-Seitz cell

    this%nrpts = 0
    ! Loop over the lattice vectors of the primitive cell
    ! that live in a supercell which is (2*ws_search_size+1)**2
    ! larger than the Born-von Karman supercell.
    ! We need to find which among these live in the Wigner-Seitz cell
    do n1 = -this%ws_search_size(1)*this%mp_grid(1), this%ws_search_size(1)*this%mp_grid(1)
      do n2 = -this%ws_search_size(2)*this%mp_grid(2), this%ws_search_size(2)*this%mp_grid(2)
        do n3 = -this%ws_search_size(3)*this%mp_grid(3), this%ws_search_size(3)*this%mp_grid(3)
          ! Loop over the lattice vectors R of the Born-von Karman supercell
          ! that contains all the points of the previous loop.
          ! There are (2*(ws_search_size+1)+1)**3 points R. R=0 corresponds to
          ! i1=i2=i3=0, or icnt=((2*(ws_search_size+1)+1)**3 + 1)/2
          icnt = 0
          do i1 = -this%ws_search_size(1) - 1, this%ws_search_size(1) + 1
            do i2 = -this%ws_search_size(2) - 1, this%ws_search_size(2) + 1
              do i3 = -this%ws_search_size(3) - 1, this%ws_search_size(3) + 1
                icnt = icnt + 1
                ! Calculate distance squared |r-R|^2
                ndiff(1) = n1 - i1*this%mp_grid(1)
                ndiff(2) = n2 - i2*this%mp_grid(2)
                ndiff(3) = n3 - i3*this%mp_grid(3)
                dist(icnt) = 0.0_REAL64
                do i = 1, 3
                  do j = 1, 3
                    dist(icnt) = dist(icnt) + real(ndiff(i), REAL64)*this%real_metric(i, j) &
                                 *real(ndiff(j), REAL64)
                  enddo
                enddo
              enddo
            enddo
          enddo
          ! AAM: On first pass, we reference unallocated variables (ndegen,irvec)
          dist_min = minval(dist)
          if (abs(dist((dist_dim + 1)/2) - dist_min) .lt. ws_distance_tol**2) then
            this%nrpts = this%nrpts + 1
            if (.not. count_pts) then
              this%ndegen(this%nrpts) = 0
              do i = 1, dist_dim
                if (abs(dist(i) - dist_min) .lt. ws_distance_tol**2) this%ndegen(this%nrpts) = this%ndegen(this%nrpts) + 1
              end do
              this%irvec(1, this%nrpts) = n1
              this%irvec(2, this%nrpts) = n2
              this%irvec(3, this%nrpts) = n3
              !
              ! Record index of r=0
!              if (n1 == 0 .and. n2 == 0 .and. n3 == 0) rpt_origin = this%nrpts
            endif
          end if

          !n3
        enddo
        !n2
      enddo
      !n1
    enddo
    !
    deallocate (dist, stat=ierr)
    if (ierr /= 0) stop 'Error in deallocating dist hamiltonian_wigner_seitz'
    if (count_pts) return

    ! Check the "sum rule"
    tot = 0.0_REAL64
    do i = 1, this%nrpts
      tot = tot + 1.0_REAL64/real(this%ndegen(i), REAL64)
    enddo

    if (iprint >= 3 ) then
      write (stdout, '(1x,i4,a,/)') this%nrpts, ' lattice points in Wigner-Seitz supercell:'
      do i = 1, this%nrpts
        write (stdout, '(4x,a,3(i3,1x),a,i2)') '  vector ', this%irvec(1, i), this%irvec(2, i), &
          this%irvec(3, i), '  degeneracy: ', this%ndegen(i)
      enddo
      write (stdout, '(1x,a,f12.3)') ' tot = ', tot
      write (stdout, '(1x,a,i12)') ' mp_grid product = ', this%mp_grid(1)*this%mp_grid(2)*this%mp_grid(3)
    endif
    if (abs(tot - real(this%mp_grid(1)*this%mp_grid(2)*this%mp_grid(3), REAL64)) > eps8) then
      write(*,*) 'ERROR in hamiltonian_wigner_seitz: error in finding Wigner-Seitz points'
      write(*,*) 'You may try to increase Wigner-Seitz search size... '
      stop
    endif


    return

  end subroutine fi_hamiltonian_wigner_seitz

  !============================================!
  
  
  
  

end module
