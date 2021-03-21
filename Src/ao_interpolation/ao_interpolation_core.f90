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

subroutine ZMul(transA,transB,A,B,C,n)
  implicit none
  integer, parameter          :: REAL64 = selected_real_kind(12)
  integer n
  character*1 transA,transB
  complex(REAL64) A(n,n),B(n,n),C(n,n)
  complex(REAL64), parameter:: zum   = (1.0D0,0.0D0)
  complex(REAL64), parameter:: zzero = (0.0D0,0.0D0)
  call zgemm(transA,transB,n,n,n,zum,A,n,B,n,zzero,C,n)
end subroutine

subroutine GetS12(S,S12,S12_inv,wrk,ev_wrk,nband)
  implicit none
  integer, parameter          :: REAL64 = selected_real_kind(12)  
  integer nband
  complex(REAL64) S(nband,nband)
  complex(REAL64) S12(nband,nband)
  complex(REAL64) S12_inv(nband,nband)
  complex(REAL64) wrk(nband,nband)
  real(REAL64) ev_wrk(nband)
  integer i
!  complex(REAL64), parameter:: zum   = (1.0D0,0.0D0)
  complex(REAL64), parameter:: zzero = (0.0D0,0.0D0)

  integer  ::  info

!  Lowdin symmetric otrthogonalization  
!  S12 is S^(-1/2)
!  S12_inv is S^(1/2) = S x S^(-1/2)
  
! S is positive definite, ev_wrk > 0
  call diag_c16(nband,S,ev_wrk,wrk,nband,info)

  
  if(info /= 0) stop


! reuse S12 
  
  S12(:,:) = zzero

  do i=1,nband
      if(ev_wrk(i) <= 0.0D0) then
        write(*,*) 'matrix not positive definite!'
        write(*,*) 'negative eigenvalue', i, ev_wrk(i)
        stop
      else
        S12(i,i) = cmplx(1.0D0/dsqrt(ev_wrk(i)),0.0D0,REAL64)      
      endif
  enddo
  
!  wrk = matmul(wrk,matmul(S12,transpose(wrk)))

! S12_inv is work array  
  call ZMul('N','C',S12,wrk,S12_inv,nband)
  call ZMul('N','N',wrk,S12_inv,S12,nband)
  
! S12_inv is calculated    
  call ZMul('N','N',S,S12,S12_inv,nband)

  
end subroutine

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

    call GetS12(S,S12,S12_inv,wrk,ev_s,nband)

    call ZMul('N','N',Hao,S12,wrk,nband)
    call ZMul('N','N',S12,wrk,Hao_tr,nband)
    
    call diag_c16(nband,Hao_tr,ev_ao,Uao,nband,info)

  
    if(info /= 0) stop



end subroutine
 
subroutine GetHpw(nband,nequal,Hao,S,ev_pw,                             &
                      Hpw,S12,ev_ao,                                    & 
&                     S12_inv,Uao,Hao_tr, wrk, ev_s)                
  implicit none
  integer, parameter          :: REAL64 = selected_real_kind(12)
  
  
! For given k-point input Hamiltonian Hao=<ao|H|ao> and 
! Overlap matrix: S = <ao|ao> in an atomic orbital basis {|ao>} 
! and the eigenvalues obtained in the plane wave calculatian ev_pw, 
! this subroutine outputs the Hamiltonian in the atomic basis 
! Hpw=<ao|Hpw|ao> which when diagonalized gives the first nequal 
! eigenvalues equal to ev_pw.

!      input
  integer nband, nequal
  
  complex(REAL64)                    :: Hao(nband,nband)
  complex(REAL64)                    :: S(nband,nband)
  real(REAL64)                       :: ev_pw(nband)

!      output
  complex(REAL64)                    :: Hpw(nband,nband)
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
  
  integer i,j

    call GetS12(S,S12,S12_inv,wrk,ev_s,nband)
    
!    Hao_tr = matmul(S12,matmul(Hao,S12))    
!    call ZMul(transA,transB,A,B,C,n)    

    call ZMul('N','N',Hao,S12,wrk,nband)
    call ZMul('N','N',S12,wrk,Hao_tr,nband)
    
    call diag_c16(nband,Hao_tr,ev_ao,Uao,nband,info)

  
    if(info /= 0) stop


        
    ! Now reconstirtute the original ham from S*S^-1/2 ev_wrk and U
    
    do i=1,nband
    do j=1,nband
      if(i==j) then
!        wrk(i,j) = cmplx(ev_wrk(i),0.0D0,REAL64)
        Hpw(i,j) = cmplx(ev_pw(i),0.0D0,REAL64)
      else
        Hpw(i,j) = cmplx(0.0D0,0.0D0,REAL64)
      endif
    enddo
    enddo

    do j=nequal+1,nband
        Hpw(j,j) = cmplx(ev_ao(j),0.0D0,REAL64)
    enddo
    
!    wrk = matmul(Uao,matmul(Hpw,transpose(Uao)))    

    call ZMul('N','C',Hpw,Uao,wrk,nband)
    call ZMul('N','N',Uao,wrk,Hpw,nband)
    
!    Hpw = matmul(S12_inv,matmul(Hpw,S12_inv))
    
    call ZMul('N','N',Hpw,S12_inv,wrk,nband)
    call ZMul('N','N',S12_inv,wrk,Hpw,nband)

  
end subroutine

subroutine OrthoH(H,S,Hw,n)
  implicit none
  integer, parameter          :: REAL64 = selected_real_kind(12)
  integer n
  complex(REAL64) :: H(n,n), S(n,n), Hw(n,n)
  
  complex(REAL64), allocatable :: S12(:,:), S12_inv(:,:) , wrk(:,:)
  real(REAL64), allocatable    :: ev_wrk(:)
    
  allocate(S12(n,n))
  allocate(S12_inv(n,n))
  allocate(wrk(n,n))
  allocate(ev_wrk(n))
    
  call GetS12(S,S12,S12_inv,wrk,ev_wrk,n)

  call ZMul('N','N',H,S12,wrk,n)
  call ZMul('N','N',S12,wrk,Hw,n)

  deallocate(S12,S12_inv,wrk,ev_wrk)

end subroutine



