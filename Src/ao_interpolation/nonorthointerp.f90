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
    deallocate(this %fiData %ndegen)
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

  subroutine NonOrthoInterpInit(this,lso,loptical,nk1,nk2,nk3,ws_n1,ws_n2,ws_n3,adot,nkpt,nband,nequal)
    implicit none
    type(noiData_t) :: this
    integer lso, loptical
    integer nk1,nk2,nk3,nkpt
    integer ws_n1, ws_n2, ws_n3
!    integer nk1,nk2,nk3,nkpt
    real(REAL64) adot(3,3)
    integer nband, nequal

    this%lso        = lso
    this%loptical   = loptical
    this%nkpt   = nkpt
    this%nband  = nband
    this%nequal = nequal

    call NonOrthoInterpInternalInit(this,nband)

    call fi_hamiltonian_setup(this%fiData,nband,nk1,nk2,nk3,ws_n1,ws_n2,ws_n3,adot)

    if (loptical==1) then
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

      call   GetHpw(this%nband,this%nequal,Hao,S,ev_pw,                            &
    &               this%fiData%Ham_k(:,:,ikpt),this%S12,this%ev_ao,               &
    &               this%S12_inv,this%Uao,this%Hao_tr, this%wrk, this%ev_s)


  end subroutine NonOrthoInterpSetGridData

   subroutine NonOrthoInterpRun2(this,Hao,S,rkpt,ev)
    implicit none
    type(noiData_t) :: this
    real(REAL64) rkpt(3)
    real(REAL64) ev(this%nband)

    integer i,j

    complex(REAL64) Hao(this%nband,this%nband)
    complex(REAL64) S(this%nband,this%nband)

    call fi_hamiltonian_get_hk(this%fiData,rkpt,Hao, this%nband)
    call fi_hamiltonian_get_sk(this%fiData,rkpt,S,   this%nband)

    do i=1,this%nband
    do j=i+1,this%nband
      Hao(j,i) = conjg(Hao(i,j))
      S(j,i)   = conjg(S(i,j))
    enddo
    enddo



    call DiagByLowdin(this%nband,Hao,S,this%S12,ev,this%S12_inv,this%Uao,this%Hao_tr, this%wrk, this%ev_s)

  end subroutine NonOrthoInterpRun2

    subroutine NonOrthoInterpWriteToFile(this)
    implicit none
    type(noiData_t) :: this
    integer ir, idir

!! experimental
!    integer, parameter           :: REAL64 = selected_real_kind(12)
!    complex(REAL64), allocatable :: S12(:,:), S12_inv(:,:) , wrk(:,:)
!    real(REAL64), allocatable    :: ev_wrk(:)
!!

    open(unit=777,file='hsr.dat',form='unformatted')

    write(777) this %nband
    write(777) this %nequal
    write(777) this %lso, this %loptical

    write(777) this %fiData %real_metric(:,:)

    write(777) this %fiData %mp_grid(1), this %fiData %mp_grid(2), this %fiData %mp_grid(3)

    write(777) this %eref

    write(777) this %fiData %ws_search_size(1), this %fiData %ws_search_size(2), this %fiData %ws_search_size(3)

    write(777) this %fiData %nrpts

    do ir=1,this %fiData %nrpts
      write(777) this %fiData %irvec(1,ir),this %fiData %irvec(2,ir),this %fiData %irvec(3,ir)
    enddo

    do ir=1,this %fiData %nrpts
      write(777) this %fiData %ndegen(ir)
    enddo

    do ir=1,this %fiData %nrpts
!      do iband=1,this %nband
!      do jband=1,this %nband
        write(777) this %fiData %ham_r(:,:,ir)
!      enddo
!      enddo
      call progress(ir,2*(this %fiData %nrpts))
      call flush(6)
    enddo

    do ir=1,this %fiData %nrpts
!      do iband=1,this %nband
!      do jband=1,this %nband
        write(777) this %fiData %S_r(:,:,ir)
!      enddo
!      enddo
      call progress(this% fiData %nrpts+ir,(2*this %fiData %nrpts))
      call flush(6)
    enddo

    if ( this %loptical == 1) then
      do ir=1,this %fiData %nrpts
        do idir=1,3
            write(777) this %fiData % dh0drk_GridR(:,:,idir,ir)
        enddo
      enddo
    endif

    close(unit=777)

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


    subroutine NonOrthoInterpFileInfo(this)
      implicit none
      type(noiData_t) :: this
      open(unit=777,file='hsr.dat',form='unformatted')
      read(777) this %nband
      read(777) this %nequal
      read(777) this %lso, this %loptical
      read(777) this %fiData % real_metric(:,:)
      read(777) this %fiData %mp_grid(1), this %fiData %mp_grid(2), this %fiData %mp_grid(3)

      write(*,*) '  hsr.dat file was prepared with the following setup:'

       write(6,*)
       write(6,*) '  real-space metric'
       write(6,*)
       write(6,'(3x ,3(1x,f16.8),5x," metric  g11,g12,g13 ")')           &
     &        this%fiData%real_metric(1,1),this%fiData%real_metric(1,2),this%fiData%real_metric(1,3)
       write(6,'(20x,2(1x,f16.8),5x," metric      g22,g23 ")')           &
     &                  this%fiData%real_metric(2,2),this%fiData%real_metric(2,3)
       write(6,'(38x,     f16.8 ,5x," metric          g33 ")')           &
     &                            this%fiData%real_metric(3,3)


      write(*,*) '  k-grid:'
      write(*,*) this %fiData %mp_grid(:)
      if(this %lso ==1) then
        write(*,*) '  Spin-Orbit: yes'
      else
        write(*,*) '  Spin-Orbit: no'
      endif
      if(this %loptical ==1) then
        write(*,*) '  dh0drk: yes'
      else
        write(*,*) '  dh0drk: no'
      endif

      close(unit=777)

    end subroutine NonOrthoInterpFileInfo



    subroutine NonOrthoInterpReadFromFile(this)
    implicit none
    type(noiData_t) :: this
    integer ir, nrpts, nband, idir

    open(unit=777,file='hsr.dat',form='unformatted')

!    write(*,*) 'deb1: '

    read(777) this %nband
    read(777) this %nequal
    read(777) this %lso, this %loptical

    read(777) this %fiData %real_metric(:,:)

    read(777) this %fiData %mp_grid(1), this %fiData %mp_grid(2), this %fiData %mp_grid(3)

    read(777) this %eref

    read(777) this %fiData %ws_search_size(1), this %fiData %ws_search_size(2), this %fiData %ws_search_size(3)

    read(777) this %fiData %nrpts

    nrpts = this %fiData %nrpts
    nband = this %nband

!    write(*,*) 'deb2: ', nband,nrpts

    call NonOrthoInterpInternalInit(this,nband)

!    write(*,*) 'deb3: '

    allocate(this %fiData %irvec(3,nrpts))
    allocate(this %fiData %ndegen(nrpts))
    allocate(this %fiData %ham_r(nband,nband,nrpts))
    allocate(this %fiData %S_r(nband,nband,nrpts))

!    write(*,*) 'deb4: '

    do ir=1,this %fiData %nrpts
      read(777) this %fiData %irvec(1,ir),this %fiData %irvec(2,ir),this %fiData %irvec(3,ir)
    enddo

!    write(*,*) 'deb5: '

    do ir=1,this %fiData %nrpts
      read(777) this %fiData %ndegen(ir)
    enddo

!    write(*,*) 'deb6: '

    do ir=1,this %fiData %nrpts
!      do iband=1,this %nband
!      do jband=1,this %nband
        read(777) this %fiData %ham_r(:,:,ir)
!      enddo
!      enddo
    enddo

!    write(*,*) 'deb7: '

    do ir=1,this %fiData %nrpts
!      do iband=1,this %nband
!      do jband=1,this %nband
        read(777) this %fiData %S_r(:,:,ir)
!      enddo
!      enddo
    enddo

!    write(*,*) 'deb8: '

    if ( this %loptical == 1) then
      allocate( this%fiData%dh0drk_GridR(this %nband,this %nband,3,this %fiData %nrpts))
      do ir=1,this %fiData %nrpts
        do idir=1,3
            read(777) this %fiData % dh0drk_GridR(:,:,idir,ir)
        enddo
      enddo
    endif




    close(unit=777)

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


subroutine PrintMatrix(iunit,lab,M,n,nmax,prtreal)
  implicit none
  integer, parameter          :: REAL64 = selected_real_kind(12)
  integer            :: iunit, n, nmax
  integer            :: i,j
  complex(REAL64)    :: M(n,n)
  logical            :: prtreal
  character(len=4)   :: lab

  return

  write(iunit,*) '*m--------------------------------------------------------------------'
  write(6,*) '*m--------------------------------------------------------------------'

  do i=1, nmax
    if(prtreal) then
      write(iunit,'(" *m",a4,i3,100f8.4)') lab,i, (ABS(M(i,j)), j=1,nmax)
      write(6    ,'(" *m",a4,i3,100f8.4)') lab,i, (ABS(M(i,j)), j=1,nmax)
!
      write(iunit,'(" *m",a4,i3,100f8.4)') lab,i, (IMAG(M(i,j)), j=1,nmax)
      write(6    ,'(" *m",a4,i3,100f8.4)') lab,i, (IMAG(M(i,j)), j=1,nmax)
      write(6,*) '*m-----'
!
    else
      write(iunit,'(" *m",a4,i3,100f8.4)') lab,i, (M(i,j), j=1,nmax)
      write(6    ,'(" *m",a4,i3,100f8.4)') lab,i, (M(i,j), j=1,nmax)
    endif
  enddo

  write(iunit,*) '*m--------------------------------------------------------------------'
  write(6,*) '*m--------------------------------------------------------------------'


end subroutine PrintMatrix
