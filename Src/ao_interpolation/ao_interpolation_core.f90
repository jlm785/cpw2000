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

subroutine GetS12(S, S12, S12_inv, wrk, ev_wrk, nband)

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  nband                           !<  number of bands
  complex(REAL64), intent(in)        ::  S(nband,nband)                  !<  Overlap natrix

! output

  complex(REAL64), intent(out)       ::  S12(nband,nband)                !<  S12 is S^(-1/2)
  complex(REAL64), intent(out)       ::  S12_inv(nband,nband)            !<  S12_inv is S^(1/2) = S x S^(-1/2)

! work arrays

  complex(REAL64), intent(out)       ::  wrk(nband,nband)
  real(REAL64), intent(out)          ::  ev_wrk(nband)

! local variables

  integer     ::  info

! counters

  integer     ::  i

!  parameters

  real(REAL64)   , parameter         :: ZERO = 0.0_REAL64
  real(REAL64)   , parameter         :: UM = 1.0_REAL64
  complex(REAL64), parameter         :: C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter         :: C_UM = cmplx(UM,ZERO,REAL64)


!  Lowdin symmetric otrthogonalization
!  S12 is S^(-1/2)
!  S12_inv is S^(1/2) = S x S^(-1/2)

! S is positive definite, ev_wrk > 0

  call diag_c16(nband, S, ev_wrk, wrk, nband, info)


  if(info /= 0) then
    write(6,*)
    write(6,*) '    STOPPED in GetS12'
    write(6,*) '    diagonalization of S failed, info = ', info
     write(6,*)

    stop

  endif

! reuse S12

  S12(:,:) = C_ZERO

  do i = 1,nband

    if(ev_wrk(i) <= ZERO) then
      write(6,*)
      write(6,*) '    STOPPED in GetS12'
      write(6,*) '    overlap matrix is not positive definite!'
      write(6,*) '    negative eigenvalue', i, ev_wrk(i)
      write(6,*)

      stop

    else
      S12(i,i) = cmplx(UM/dsqrt(ev_wrk(i)),ZERO,REAL64)
    endif

  enddo

!  wrk = matmul(wrk,matmul(S12,transpose(wrk)))

! S12_inv is work array here

!  call ZMul('N','C',S12,wrk,S12_inv,nband)
  call zgemm('N','C', nband,nband,nband, C_UM,S12,nband, wrk,nband,      &
             C_ZERO,S12_inv,nband)

!  call ZMul('N','N',wrk,S12_inv,S12,nband)
  call zgemm('N','N', nband,nband,nband, C_UM,wrk,nband, S12_inv,nband,  &
              C_ZERO,S12,nband)

! S12_inv is calculated now

!  call ZMul('N','N',S,S12,S12_inv,nband)
  call zgemm('N','N', nband,nband,nband, C_UM,S,nband, S12,nband,        &
                C_ZERO,S12_inv,nband)

  return

end subroutine GetS12



!   THIS SUBROUNTINE IS NOT USED   NonOrthoInterpRun2   IS NEVER CALLED

!   subroutine DiagByLowdin(nband,Hao,S,S12,ev_ao,S12_inv,Uao,Hao_tr, wrk, ev_s)


subroutine GetHpw(nband, nequal, Hao, S, ev_pw,                          &
                  Hpw, S12, ev_ao,                                       &
                  S12_inv, Uao, Hao_tr, wrk, ev_s)

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! For given k-point input Hamiltonian Hao=<ao|H|ao> and
! Overlap matrix: S = <ao|ao> in an atomic orbital basis {|ao>}
! and the eigenvalues obtained in the plane wave calculatian ev_pw,
! this subroutine outputs the Hamiltonian in the atomic basis
! Hpw=<ao|Hpw|ao> which when diagonalized gives the first nequal
! eigenvalues equal to ev_pw.

! input

  integer, intent(in)                ::  nband                           !<  number of bands
  integer, intent(in)                ::  nequal                          !<  number of eigenvalues that should be equal

  complex(REAL64), intent(in)        ::  Hao(nband,nband)                !<  Hamiltonian in the atomic orbital basis
  complex(REAL64), intent(in)        ::  S(nband,nband)                  !<  Overlap matrix

  real(REAL64), intent(in)           ::  ev_pw(nband)                    !<  Eigenvalues of the plane-waves

! output

  complex(REAL64), intent(out)       ::  Hpw(nband,nband)                !<  MTB hamiltonian
  complex(REAL64), intent(out)       ::  S12(nband,nband)                !<  S12 is S^(-1/2)
  real(REAL64), intent(out)          ::  ev_ao(nband)                    !<  Eigenvalues of the atomic orbital Hamiltonian


! work arrays, but may have physical content

  complex(REAL64), intent(out)       ::  S12_inv(nband,nband)            !<  S12_inv is S^(1/2) = S x S^(-1/2) (work)
  complex(REAL64), intent(out)       ::  Uao(nband,nband)                !<  Wave-functions in the orthogonalized atomic basis (work)
  complex(REAL64), intent(out)       ::  Hao_tr(nband,nband)             !<  Hamiltonian in orthogonalized atomic basis (work)
  complex(REAL64), intent(out)       ::  wrk(nband,nband)                !<  Work array

  real(REAL64), intent(out)          ::  ev_s(nband)                     !<  Eigenvalues of matrix S (work?)

! local variables

  integer       ::  info

!  counters

  integer       ::  i,j

! parameters

  real(REAL64)   , parameter         :: ZERO = 0.0_REAL64
  real(REAL64)   , parameter         :: UM = 1.0_REAL64
  complex(REAL64), parameter         :: C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter         :: C_UM = cmplx(UM,ZERO,REAL64)


  call GetS12(S, S12, S12_inv, wrk, ev_s, nband)

!    Hao_tr = matmul(S12,matmul(Hao,S12))
!    call ZMul(transA,transB,A,B,C,n)

!    call ZMul('N','N',Hao,S12,wrk,nband)
  call zgemm('N','N', nband,nband,nband, C_UM,Hao,nband, S12,nband,        &
                C_ZERO,wrk,nband)
!    call ZMul('N','N',S12,wrk,Hao_tr,nband)
  call zgemm('N','N', nband,nband,nband, C_UM,S12,nband, wrk,nband,        &
                C_ZERO,Hao_tr,nband)

  call diag_c16(nband, Hao_tr, ev_ao, Uao, nband, info)


  if(info /= 0) then
    write(6,*)
    write(6,*) '    STOPPED in GetHpw'
    write(6,*) '    diagonalization of S failed, info = ', info
    write(6,*)

    stop

  endif

! Now reconstirtute the original ham from S*S^-1/2 ev_wrk and U

  do i=1,nband
  do j=1,nband
    if(i==j) then
!      wrk(i,j) = cmplx(ev_wrk(i),0.0D0,REAL64)
      Hpw(i,j) = cmplx(ev_pw(i),ZERO,REAL64)
    else
      Hpw(i,j) = C_ZERO
    endif
  enddo
  enddo

  do j = nequal+1,nband
      Hpw(j,j) = cmplx(ev_ao(j),ZERO,REAL64)
  enddo

!    wrk = matmul(Uao,matmul(Hpw,transpose(Uao)))

!    call ZMul('N','C',Hpw,Uao,wrk,nband)
  call zgemm('N','C', nband,nband,nband, C_UM,Hpw,nband, Uao,nband,        &
                C_ZERO,wrk,nband)
!    call ZMul('N','N',Uao,wrk,Hpw,nband)
  call zgemm('N','N', nband,nband,nband, C_UM,Uao,nband, wrk,nband,        &
                C_ZERO,Hpw,nband)

!    Hpw = matmul(S12_inv,matmul(Hpw,S12_inv))

!    call ZMul('N','N',Hpw,S12_inv,wrk,nband)
  call zgemm('N','N', nband,nband,nband, C_UM,Hpw,nband, S12_inv,nband,        &
                C_ZERO,wrk,nband)
!    call ZMul('N','N',S12_inv,wrk,Hpw,nband)
  call zgemm('N','N', nband,nband,nband, C_UM,S12_inv,nband, wrk,nband,        &
                C_ZERO,Hpw,nband)


end subroutine GetHpw



subroutine OrthoH(H,S,Hw,n)

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer n
  complex(REAL64) :: H(n,n), S(n,n), Hw(n,n)

  complex(REAL64), allocatable :: S12(:,:), S12_inv(:,:) , wrk(:,:)
  real(REAL64), allocatable    :: ev_wrk(:)

  real(REAL64)   , parameter         :: ZERO = 0.0_REAL64
  real(REAL64)   , parameter         :: UM = 1.0_REAL64
  complex(REAL64), parameter         :: C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter         :: C_UM = cmplx(UM,ZERO,REAL64)

  allocate(S12(n,n))
  allocate(S12_inv(n,n))
  allocate(wrk(n,n))
  allocate(ev_wrk(n))

  call GetS12(S,S12,S12_inv,wrk,ev_wrk,n)

!  call ZMul('N','N',H,S12,wrk,n)
  call zgemm('N','N', n,n,n, C_UM,H,n, S12,n,        &
                C_ZERO,wrk,n)
!  call ZMul('N','N',S12,wrk,Hw,n)
  call zgemm('N','N', n,n,n, C_UM,S12,n, wrk,n,      &
                C_ZERO,Hw,n)

  deallocate(S12,S12_inv,wrk,ev_wrk)

  return

end subroutine OrthoH




subroutine DiagByLowdin2(nband, Hao, S, ev, psi)

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  nband                           !<  number of bands

  complex(REAL64), intent(in)        ::  Hao(nband,nband)                !<  Hamiltonian in atomic basis
  complex(REAL64), intent(in)        ::  S(nband,nband)                  !<  Overlap natrix

! output

  real(REAL64), intent(out)          ::  ev(nband)                       !<  Eigenvalues
  complex(REAL64) , intent(out)      ::  psi(nband,nband)                !<  Wave-functions in atomic basis

! allocatable work arrays

  complex(REAL64), allocatable       ::  S12(:,:)

  complex(REAL64), allocatable       :: S12_inv(:,:)
  complex(REAL64), allocatable       :: Hao_tr(:,:)
  complex(REAL64), allocatable       :: wrk(:,:)

  real(REAL64), allocatable          :: ev_s(:)

  integer  ::  info

  real(REAL64)   , parameter         :: ZERO = 0.0_REAL64
  real(REAL64)   , parameter         :: UM = 1.0_REAL64
  complex(REAL64), parameter         :: C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter         :: C_UM = cmplx(UM,ZERO,REAL64)

! begin

  allocate(S12(nband,nband))

  allocate(S12_inv(nband,nband))
  allocate(Hao_tr(nband,nband))
  allocate(wrk(nband,nband))

  allocate(ev_s(nband))


  call GetS12(S, S12, S12_inv, wrk, ev_s, nband)

!    call ZMul('N','N',Hao,S12,wrk,nband)
  call zgemm('N','N', nband,nband,nband, C_UM,Hao,nband, S12,nband,        &
                C_ZERO,wrk,nband)
!    call ZMul('N','N',S12,wrk,Hao_tr,nband)
  call zgemm('N','N', nband,nband,nband, C_UM,S12,nband, wrk,nband,        &
                C_ZERO,Hao_tr,nband)

  call diag_c16(nband, Hao_tr, ev, psi, nband, info)

  if(info /= 0) then

    write(6,*)
    write(6,*) '     STOPPED IN DiagByLowdin2'
    write(6,*) '     Hao_tr diagonalization failed, info = ', info
    write(6,*)

    stop

  endif

  deallocate(S12,  S12_inv, Hao_tr, wrk, ev_s)

  return

end subroutine DiagByLowdin2



subroutine ZMul(transA,transB,A,B,C,n)
  implicit none
  integer, parameter          :: REAL64 = selected_real_kind(12)
  integer n
  character(len=1)     ::  transA,transB
  complex(REAL64) A(n,n),B(n,n),C(n,n)
  complex(REAL64), parameter:: zum   = (1.0D0,0.0D0)
  complex(REAL64), parameter:: zzero = (0.0D0,0.0D0)
  if(n>2) stop
end subroutine


!   THIS SUBROUNTINE IS NOT USED   NonOrthoInterpRun2   IS NEVER CALLED

!   subroutine DiagByLowdin(nband,Hao,S,S12,ev_ao,S12_inv,Uao,Hao_tr, wrk, ev_s)

subroutine DiagByLowdin(nband,Hao,S,S12,ev_ao,S12_inv,Uao,Hao_tr, wrk, ev_s)
  implicit none
  integer, parameter          :: REAL64 = selected_real_kind(12)


!      input
  integer nband

  complex(REAL64)                    :: Hao(nband,nband)
  complex(REAL64)                    :: S(nband,nband)

!      output
  complex(REAL64)                    :: S12(nband,nband)
  real(REAL64)                       :: ev_ao(nband)


!      wrk
  complex(REAL64)                    :: S12_inv(nband,nband)
  complex(REAL64)                    :: Uao(nband,nband)
  complex(REAL64)                    :: Hao_tr(nband,nband)
  complex(REAL64)                    :: wrk(nband,nband)

  real(REAL64)                       :: ev_s(nband)

  integer  ::  info

!  misc

!  integer i,j




    if(nband > 2) stop



end subroutine
