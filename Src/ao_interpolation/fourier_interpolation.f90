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
! Debug print in fi_hamiltonian_setup removed. 4 October 2026. JLM+claude
! Wigner-Seitz search replaced by the nearest images of each hopping between
! atomic sites (near_dist_image_all), irvec are the grid classes,
! ndegen and ws_search_size removed. 23 October 2026. JLM+claude

  implicit none
  integer, parameter,private          :: REAL64 = selected_real_kind(12)

  type fiData_t
    complex(REAL64), allocatable  :: ham_k(:,:,:)
    complex(REAL64), allocatable  :: S_k(:,:,:)
    complex(REAL64), allocatable  :: ham_r(:,:,:)
    complex(REAL64), allocatable  :: S_r(:,:,:)

    integer                       :: mp_grid(3)
    real(REAL64)                  :: real_metric(3,3)

    integer,  allocatable         :: irvec(:,:)                ! R of the grid classes, 0 <= R_i < mp_grid(i)

    integer                       :: natsite                   ! number of atomic sites
    real(REAL64), allocatable     :: tau(:,:)                  ! lattice coordinates of the sites
    integer,  allocatable         :: isite(:)                  ! site of the basis function
    integer,  allocatable         :: nimg(:,:,:)               ! number of nearest images of R + tau_b - tau_a
    integer,  allocatable         :: irimg(:,:,:,:,:)          ! lattice vectors of the nearest images (R + T)


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

  subroutine fi_hamiltonian_phase(this,rkpt,ph)

!   phase factors sum_j exp(i 2 pi k.(R+T_j)) / nimg  for each pair of
!   sites and grid class, the images T_j are the nearest ones of the
!   hopping R + tau_b - tau_a.  23 October 2026. JLM+claude

    implicit none
    integer, parameter          :: REAL64 = selected_real_kind(12)
    type(fiData_t):: this
    real(REAL64), parameter     :: pi=3.141592653589793238462643383279_REAL64
    real(REAL64), parameter     :: twopi=2.0_REAL64*pi
    complex(REAL64), parameter  :: cmplx_0 = (0.0_REAL64,0.0_REAL64)

    real(REAL64) rkpt(3)
    complex(REAL64) ph(this%natsite,this%natsite,this%nrpts)

    real(REAL64) rdotk
    integer ir, ia, ib, j

    do ir=1,this%nrpts
      do ib=1,this%natsite
      do ia=1,this%natsite
        ph(ia,ib,ir) = cmplx_0
        do j=1,this%nimg(ia,ib,ir)
          rdotk=twopi*(rkpt(1)*this%irimg(1,j,ia,ib,ir) +                &
                       rkpt(2)*this%irimg(2,j,ia,ib,ir) +                &
                       rkpt(3)*this%irimg(3,j,ia,ib,ir))
          ph(ia,ib,ir) = ph(ia,ib,ir) + cmplx(cos(rdotk),sin(rdotk),REAL64)
        enddo
        ph(ia,ib,ir) = ph(ia,ib,ir) / this%nimg(ia,ib,ir)
      enddo
      enddo
    enddo

  end subroutine fi_hamiltonian_phase


  subroutine fi_hamiltonian_get_hk(this,rkpt,ham, nband)
    implicit none
    integer, parameter          :: REAL64 = selected_real_kind(12)
    type(fiData_t):: this
    complex(REAL64), parameter  :: cmplx_0 = (0.0_REAL64,0.0_REAL64)

    integer nband
    real(REAL64) rkpt(3)
    complex(REAL64) ham(nband,nband)

    complex(REAL64), allocatable :: ph(:,:,:)

    integer loop_rpt, i, j

    allocate(ph(this%natsite,this%natsite,this%nrpts))
    call fi_hamiltonian_phase(this,rkpt,ph)

    ham=cmplx_0
    do loop_rpt=1,this%nrpts
      do j=1,nband
      do i=1,nband
        ham(i,j)=ham(i,j)+ph(this%isite(i),this%isite(j),loop_rpt)*this%ham_r(i,j,loop_rpt)
      enddo
      enddo
    end do

    deallocate(ph)

  end subroutine fi_hamiltonian_get_hk


  subroutine fi_hamiltonian_get_sk(this,rkpt,S, nband)
    implicit none
    integer, parameter          :: REAL64 = selected_real_kind(12)
    type(fiData_t):: this
    complex(REAL64), parameter  :: cmplx_0 = (0.0_REAL64,0.0_REAL64)

    integer nband
    real(REAL64) rkpt(3)
    complex(REAL64) S(nband,nband)

    complex(REAL64), allocatable :: ph(:,:,:)

    integer loop_rpt, i, j

    allocate(ph(this%natsite,this%natsite,this%nrpts))
    call fi_hamiltonian_phase(this,rkpt,ph)

    S=cmplx_0
    do loop_rpt=1,this%nrpts
      do j=1,nband
      do i=1,nband
        S(i,j)=S(i,j)+ph(this%isite(i),this%isite(j),loop_rpt)*this%S_r(i,j,loop_rpt)
      enddo
      enddo
    end do

    deallocate(ph)

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
    complex(REAL64), parameter  :: cmplx_0 = (0.0_REAL64,0.0_REAL64)
    integer nband
    real(REAL64) rkpt(3)

    complex(REAL64) Op_k(nband,nband,3)
    complex(REAL64) Op_R(nband,nband,3,this%nrpts)

    complex(REAL64), allocatable :: ph(:,:,:)
    complex(REAL64) fac

    integer loop_rpt, i, j

    allocate(ph(this%natsite,this%natsite,this%nrpts))
    call fi_hamiltonian_phase(this,rkpt,ph)

    Op_k(:,:,:)=cmplx_0

    do loop_rpt=1,this%nrpts
      do j=1,nband
      do i=1,nband
        fac=ph(this%isite(i),this%isite(j),loop_rpt)
        Op_k(i,j,1)=Op_k(i,j,1)+fac*Op_r(i,j,1,loop_rpt)
        Op_k(i,j,2)=Op_k(i,j,2)+fac*Op_r(i,j,2,loop_rpt)
        Op_k(i,j,3)=Op_k(i,j,3)+fac*Op_r(i,j,3,loop_rpt)
      enddo
      enddo
    end do

    deallocate(ph)

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

      call progress(loop_rpt,this%nrpts)
      call flush(6)

    enddo

  end subroutine fi_hamiltonian_set_in_Rgrid


  subroutine fi_hamiltonian_setup(this,nband,n1,n2,n3,adot,natsite,tau,isite)

!   grid, nearest images and space for H and S.  The images replace the
!   Wigner-Seitz search.  23 October 2026. JLM+claude

    implicit none
    type(fiData_t):: this

    integer nband
    integer n1,n2,n3
    real(REAL64) adot(3,3)
    integer natsite
    real(REAL64) tau(3,natsite)
    integer isite(nband)

    integer i,j,k,ikpt

    complex(REAL64), parameter  :: cmplx_0 = (0.0_REAL64,0.0_REAL64)

    this%num_kpts = n1*n2*n3
    allocate(this%kpt_latt(3,this%num_kpts))

    this%mp_grid(1) = n1
    this%mp_grid(2) = n2
    this%mp_grid(3) = n3

    this%real_metric(:,:) = adot(:,:)

    call fi_hamiltonian_images(this,nband,natsite,tau,isite)

    allocate(this%ham_r(nband,nband,this%nrpts))
    allocate(this%S_r(nband,nband,this%nrpts))
    this%ham_r(:,:,:)= cmplx_0
    this%S_r(:,:,:)  = cmplx_0

    allocate(this%ham_k(nband,nband,this%num_kpts))
    allocate(this%S_k(nband,nband,this%num_kpts))
    this%ham_k(:,:,:)= cmplx_0
    this%S_k(:,:,:)  = cmplx_0

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

    return
  end subroutine fi_hamiltonian_setup


  subroutine fi_hamiltonian_images(this,nband,natsite,tau,isite)

!   The inverse Fourier transform on the grid gives the hopping only
!   modulo the lattice of the Born-von Karman supercell.  For each
!   grid class R and pair of sites a, b the hopping R + tau_b - tau_a
!   is assigned to its nearest images in that supercell, all of them
!   with equal weights if several are at the same distance.  This
!   keeps the symmetry of the interpolated hamiltonian.
!   Needs mp_grid and real_metric.  23 October 2026. JLM+claude

    implicit none
    type(fiData_t):: this

    integer nband
    integer natsite
    real(REAL64) tau(3,natsite)
    integer isite(nband)

    integer, parameter          :: MXDIMG = 8
    real(REAL64), parameter     :: TOL = 1.0E-4_REAL64            ! tolerance for equal distances (bohr)

    real(REAL64)   ::  adots(3,3), adotb(3,3), dr(3), distmin
    integer        ::  mtb(3,3), itrans(3,MXDIMG), nimg
    integer        ::  i, j, k, ir, ia, ib, n

    this%nrpts = this%mp_grid(1)*this%mp_grid(2)*this%mp_grid(3)

    allocate(this%irvec(3,this%nrpts))
    ir = 0
    do i=0,this%mp_grid(1)-1
    do j=0,this%mp_grid(2)-1
    do k=0,this%mp_grid(3)-1
      ir = ir + 1
      this%irvec(1,ir) = i
      this%irvec(2,ir) = j
      this%irvec(3,ir) = k
    enddo
    enddo
    enddo

    this%natsite = natsite
    allocate(this%tau(3,natsite))
    allocate(this%isite(nband))
    this%tau(:,:) = tau(:,:)
    this%isite(:) = isite(:)

!   reduced cell of the supercell

    do j=1,3
    do i=1,3
      adots(i,j) = this%mp_grid(i)*this%real_metric(i,j)*this%mp_grid(j)
    enddo
    enddo

    call metric_buerger(adots, adotb, mtb)

    allocate(this%nimg(natsite,natsite,this%nrpts))
    allocate(this%irimg(3,MXDIMG,natsite,natsite,this%nrpts))

    do ir=1,this%nrpts
      do ib=1,natsite
      do ia=1,natsite
        do i=1,3
          dr(i) = (this%irvec(i,ir) + tau(i,ib) - tau(i,ia)) / this%mp_grid(i)
        enddo

        call near_dist_image_all(nimg, itrans, distmin,                  &
            adotb, mtb, dr, TOL,                                         &
            MXDIMG)

        this%nimg(ia,ib,ir) = nimg
        do n=1,nimg
          do i=1,3
            this%irimg(i,n,ia,ib,ir) = this%irvec(i,ir) + this%mp_grid(i)*itrans(i,n)
          enddo
        enddo
      enddo
      enddo
    enddo

    return
  end subroutine fi_hamiltonian_images


  subroutine fi_hamiltonian_expand(this,nband,Op_r,nexp,irexp,Opexp)

!   Operator in the list of all distinct lattice vectors of the
!   nearest images, as needed for files with a single list of
!   lattice vectors (Wannier90 _hr format, ndegen = 1).
!   23 October 2026. JLM+claude

    implicit none
    type(fiData_t):: this

    integer nband
    complex(REAL64) Op_r(nband,nband,this%nrpts)
    integer nexp
    integer, allocatable :: irexp(:,:)
    complex(REAL64), allocatable :: Opexp(:,:,:)

    complex(REAL64), parameter  :: cmplx_0 = (0.0_REAL64,0.0_REAL64)

    integer, allocatable :: indx(:,:,:)
    integer   ::  nmin(3), nmax(3), iv(3)
    integer   ::  i, j, ir, ia, ib, n, ie, i1, i2, i3

!   box containing all the images

    nmin(:) = 0
    nmax(:) = 0
    do ir=1,this%nrpts
      do ib=1,this%natsite
      do ia=1,this%natsite
        do n=1,this%nimg(ia,ib,ir)
          do i=1,3
            nmin(i) = min(nmin(i),this%irimg(i,n,ia,ib,ir))
            nmax(i) = max(nmax(i),this%irimg(i,n,ia,ib,ir))
          enddo
        enddo
      enddo
      enddo
    enddo

    allocate(indx(nmin(1):nmax(1),nmin(2):nmax(2),nmin(3):nmax(3)))
    indx = 0

    nexp = 0
    do ir=1,this%nrpts
      do ib=1,this%natsite
      do ia=1,this%natsite
        do n=1,this%nimg(ia,ib,ir)
          iv(:) = this%irimg(:,n,ia,ib,ir)
          if(indx(iv(1),iv(2),iv(3)) == 0) then
            nexp = nexp + 1
            indx(iv(1),iv(2),iv(3)) = nexp
          endif
        enddo
      enddo
      enddo
    enddo

    allocate(irexp(3,nexp))
    allocate(Opexp(nband,nband,nexp))
    Opexp(:,:,:) = cmplx_0

    do i1=nmin(1),nmax(1)
    do i2=nmin(2),nmax(2)
    do i3=nmin(3),nmax(3)
      ie = indx(i1,i2,i3)
      if(ie > 0) then
        irexp(1,ie) = i1
        irexp(2,ie) = i2
        irexp(3,ie) = i3
      endif
    enddo
    enddo
    enddo

    do ir=1,this%nrpts
      do j=1,nband
      do i=1,nband
        ia = this%isite(i)
        ib = this%isite(j)
        do n=1,this%nimg(ia,ib,ir)
          iv(:) = this%irimg(:,n,ia,ib,ir)
          ie = indx(iv(1),iv(2),iv(3))
          Opexp(i,j,ie) = Opexp(i,j,ie) + Op_r(i,j,ir) / this%nimg(ia,ib,ir)
        enddo
      enddo
      enddo
    enddo

    deallocate(indx)

    return
  end subroutine fi_hamiltonian_expand


end module
