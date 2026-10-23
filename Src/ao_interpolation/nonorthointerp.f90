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

module NonOrthoInterp
! Atomic sites (tau, isite) replace the Wigner-Seitz search size, new hsr.dat
! content. 23 October 2026. JLM+claude
! Write and read on an open unit (file hsr2.dat with a header), FileInfo
! removed. 23 October 2026. JLM+claude
  use FourierInterpolation
  implicit none
  integer, parameter,private          :: REAL64 = selected_real_kind(12)

  type noiData_t
    type (fiData_t):: fiData
    integer nband,nkpt,nequal
    integer lso
    integer loptical

    real(REAL64)     :: eref

    complex(REAL64), allocatable:: Hao_tr(:,:)
    complex(REAL64), allocatable:: Uao(:,:)
    complex(REAL64), allocatable:: wrk(:,:),S12(:,:),S12_inv(:,:)
    real(REAL64),    allocatable:: ev_ao(:),ev_s(:)
  end type

  contains

  subroutine NonOrthoInterpFinish(this)
    implicit none
    type(noiData_t) :: this

    deallocate(this%Hao_tr)

    deallocate(this%S12)
    deallocate(this%S12_inv)
    deallocate(this%Uao)

    deallocate(this%wrk)

    deallocate(this%ev_ao)
    deallocate(this%ev_s)

    deallocate(this %fiData %irvec)
    deallocate(this %fiData %tau)
    deallocate(this %fiData %isite)
    deallocate(this %fiData %nimg)
    deallocate(this %fiData %irimg)
    deallocate(this %fiData %ham_r)
    deallocate(this %fiData %S_r)
    if (allocated(this %fiData %kpt_latt)) then
      deallocate(this %fiData %kpt_latt)
    endif

  end subroutine NonOrthoInterpFinish

  subroutine NonOrthoInterpInternalInit(this,nband)
    implicit none
    type(noiData_t) :: this
    integer nband

    allocate(this%Hao_tr(nband,nband))

    allocate(this%S12(nband,nband))
    allocate(this%S12_inv(nband,nband))
    allocate(this%Uao(nband,nband))

    allocate(this%wrk(nband,nband))

    allocate(this%ev_ao(nband))
    allocate(this%ev_s(nband))

    this %eref = 0.0

  end subroutine NonOrthoInterpInternalInit

  subroutine NonOrthoInterpInit(this,lso,loptical,nk1,nk2,nk3,adot,natsite,tau,isite,nkpt,nband,nequal)
    implicit none
    type(noiData_t) :: this
    logical, intent(in)         ::  lso                                  !  if true use spin orbit
    logical, intent(in)         ::  loptical                             !  if true calculate optical matrix elements
    integer nk1,nk2,nk3,nkpt
    integer natsite                                                      ! number of atomic sites
!    integer nk1,nk2,nk3,nkpt
    real(REAL64) adot(3,3)
    integer nband, nequal
    real(REAL64) tau(3,natsite)                                          ! lattice coordinates of the sites
    integer isite(nband)                                                 ! site of the basis function

    if(lso) then
      this%lso = 1
    else
      this%lso = 0
    endif
    if(loptical) then
      this%loptical   = 1
    else
      this%loptical   = 0
    endif
    this%nkpt   = nkpt
    this%nband  = nband
    this%nequal = nequal

    call NonOrthoInterpInternalInit(this,nband)

    call fi_hamiltonian_setup(this%fiData,nband,nk1,nk2,nk3,adot,natsite,tau,isite)

    if (this%loptical == 1) then
      write(6,*) 'allocating space for optical properties computation'

      allocate( this%fiData%dh0drk_GridK(nband,nband,3,this %fiData %num_kpts))
      allocate( this%fiData%dh0drk_GridR(nband,nband,3,this %fiData %nrpts))

    endif

  end subroutine NonOrthoInterpInit

  subroutine NonOrthoInterpSetGridData(this,ikpt,Hao,S,ev_pw)
    implicit none
    type(noiData_t) :: this
    integer ikpt

    complex(REAL64) Hao(this%nband,this%nband)
    complex(REAL64) S(this%nband,this%nband)
    real(REAL64) ev_pw(this%nband)

    this % fiData % S_k(:,:,ikpt) = S(:,:)

      call   ao_int_GetHpw(this%nband,this%nequal,Hao,S,ev_pw,                     &
    &               this%fiData%Ham_k(:,:,ikpt),this%S12,this%ev_ao,               &
    &               this%S12_inv,this%Uao,this%Hao_tr, this%wrk, this%ev_s)


  end subroutine NonOrthoInterpSetGridData


    subroutine NonOrthoInterpWriteToFile(this,io)

!   writes to the open unformatted unit io (after the header of hsr2.dat)

    implicit none
    type(noiData_t) :: this
    integer io
    integer ir, idir

!! experimental
!    integer, parameter           :: REAL64 = selected_real_kind(12)
!    complex(REAL64), allocatable :: S12(:,:), S12_inv(:,:) , wrk(:,:)
!    real(REAL64), allocatable    :: ev_wrk(:)
!!

    write(io) this %nband
    write(io) this %nequal
    write(io) this %lso, this %loptical

    write(io) this %fiData %real_metric(:,:)

    write(io) this %fiData %mp_grid(1), this %fiData %mp_grid(2), this %fiData %mp_grid(3)

    write(io) this %eref

    write(io) this %fiData %natsite
    write(io) this %fiData %tau(:,:)
    write(io) this %fiData %isite(:)

    do ir=1,this %fiData %nrpts
!      do iband=1,this %nband
!      do jband=1,this %nband
        write(io) this %fiData %ham_r(:,:,ir)
!      enddo
!      enddo
      call progress(ir,2*(this %fiData %nrpts))
      call flush(6)
    enddo

    do ir=1,this %fiData %nrpts
!      do iband=1,this %nband
!      do jband=1,this %nband
        write(io) this %fiData %S_r(:,:,ir)
!      enddo
!      enddo
      call progress(this% fiData %nrpts+ir,(2*this %fiData %nrpts))
      call flush(6)
    enddo

    if ( this %loptical == 1) then
      do ir=1,this %fiData %nrpts
        do idir=1,3
            write(io) this %fiData % dh0drk_GridR(:,:,idir,ir)
        enddo
      enddo
    endif


!!   experimental section

!    subroutine GetS12(S,S12,S12_inv,wrk,ev_wrk,nband)

!    allocate(S12(this%nband,this%nband))
!    allocate(S12_inv(this%nband,this%nband))
!    allocate(wrk(this%nband,this%nband))
!    allocate(ev_wrk(this%nband))

!    do ir=1,this %fiData %nrpts
!      write(*,*) 'ir=', ir
!      call GetS12( this %fiData %S_r(:,:,ir),S12,S12_inv,wrk,ev_wrk,this %nband)
!    enddo


    end subroutine NonOrthoInterpWriteToFile





    subroutine NonOrthoInterpReadFromFile(this,io)

!   reads from the open unformatted unit io (after the header of hsr2.dat)

    implicit none
    type(noiData_t) :: this
    integer io
    integer ir, nrpts, nband, idir
    integer natsite
    real(REAL64), allocatable :: tau(:,:)
    integer, allocatable :: isite(:)

!    write(*,*) 'deb1: '

    read(io) this %nband
    read(io) this %nequal
    read(io) this %lso, this %loptical

    read(io) this %fiData %real_metric(:,:)

    read(io) this %fiData %mp_grid(1), this %fiData %mp_grid(2), this %fiData %mp_grid(3)

    read(io) this %eref

    read(io) natsite
    allocate(tau(3,natsite))
    allocate(isite(this %nband))
    read(io) tau(:,:)
    read(io) isite(:)

    nband = this %nband

!   grid classes and nearest images

    call fi_hamiltonian_images(this %fiData,nband,natsite,tau,isite)

    deallocate(tau)
    deallocate(isite)

    nrpts = this %fiData %nrpts

!    write(*,*) 'deb2: ', nband,nrpts

    call NonOrthoInterpInternalInit(this,nband)

!    write(*,*) 'deb3: '

    allocate(this %fiData %ham_r(nband,nband,nrpts))
    allocate(this %fiData %S_r(nband,nband,nrpts))

!    write(*,*) 'deb4: '

!    write(*,*) 'deb6: '

    do ir=1,this %fiData %nrpts
!      do iband=1,this %nband
!      do jband=1,this %nband
        read(io) this %fiData %ham_r(:,:,ir)
!      enddo
!      enddo
    enddo

!    write(*,*) 'deb7: '

    do ir=1,this %fiData %nrpts
!      do iband=1,this %nband
!      do jband=1,this %nband
        read(io) this %fiData %S_r(:,:,ir)
!      enddo
!      enddo
    enddo

!    write(*,*) 'deb8: '

    if ( this %loptical == 1) then
      allocate( this%fiData%dh0drk_GridR(this %nband,this %nband,3,this %fiData %nrpts))
      do ir=1,this %fiData %nrpts
        do idir=1,3
            read(io) this %fiData % dh0drk_GridR(:,:,idir,ir)
        enddo
      enddo
    endif





    end subroutine NonOrthoInterpReadFromFile


   subroutine NonOrthoInterpFinishSetup(this)
    implicit none
    type(noiData_t) :: this
    write(*,*) '   Transforming H(k),S(k) to H(R) and S(R)'
    call fi_hamiltonian_get_hr(this%fiData)
    write(*,*) '   Done'
    if ( this %loptical == 1) then
    write(*,*) '   Transforming dH0drk(k) to dH0drk(R)'
      call fi_hamiltonian_set_in_Rgrid(this %fiData,this %fiData %dh0drk_GridK,this %fiData %dh0drk_GridR, this%nband)
    write(*,*) '   Done'
    endif

!   free some memory dealocating Hk and Sk
!
    write(*,*) '   Freeing Space'
    deallocate(this %fiData %Ham_k)
    deallocate(this %fiData %S_k)
    if ( this %loptical == 1) then
      deallocate(this %fiData %dh0drk_GridK)
    endif

  end subroutine NonOrthoInterpFinishSetup

   subroutine NonOrthoInterpRun(this,rkpt,ev)
    implicit none
    type(noiData_t) :: this
    real(REAL64) rkpt(3)
    real(REAL64) ev(this%nband)

    integer  ::  info

    call fi_hamiltonian_get_hk(this%fiData,rkpt,this%Hao_tr, this%nband)
    call fi_hamiltonian_get_sk(this%fiData,rkpt,this%Uao, this%nband)
    call diag_c16_gen(this%nband,this%Hao_tr,this%Uao,ev,this%wrk,this%nband,info)


    if(info /= 0) stop


    end subroutine NonOrthoInterpRun


end module

