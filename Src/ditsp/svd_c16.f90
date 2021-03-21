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

!>     Interface to LAPACK single value decomposition to orthogonalize the vectors psi

       subroutine svd_c16(psi,singval,mtxd,neig,info,mxddim,mxdbnd)

!      Written February 2014. clr, jlm
!      Modified, documentation, January 2020. JLM
!      Modified, outputs the singular values, 19 August 2020. JLM
!      copyright INESC-MN/Jose Luis Martins and Carlos Loia Reis

!      version 4.98

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxddim                     !<  array dimension of plane-waves
       integer, intent(in)                ::  mxdbnd                     !<  array dimension for number of bands

       integer, intent(in)                ::  mtxd                       !<  dimension of eigenvectors
       integer, intent(in)                ::  neig                       !<  number of eigenvectors (requested on input, modified by degeneracies on output)

!      input and output

       complex(REAL64), intent(inout)     ::  psi(mxddim,mxdbnd)         !<  component j of reference eigenvector

!      output

       real(REAL64), intent(out)          ::  singval(mxdbnd)            !<  singular values
       integer, intent(out)               ::  info                       !<  info = 0 if successful
       
!      local allocatable arrays

       complex(REAL64), allocatable       ::  u(:,:)                     !  left singular vectors (not used)
       complex(REAL64), allocatable       ::  vt(:,:)                    !  hermitian of v

       complex(REAL64), allocatable       ::  work(:)
       real(REAL64), allocatable          ::  rwork(:)
       integer, allocatable               ::  iwork(:)  
       
!      local variables

       integer  ::  lwork, lrwork


!       lwork = 2*neig + mtxd

       allocate(u(mxddim,neig))
       allocate(vt(neig,neig))
       allocate(work(1))
       allocate(rwork(1))
       allocate(iwork(8*neig))

       call zgesdd('O',mtxd,neig,psi,mxddim,singval,u,mxddim,            &
     &        vt,neig,work,-1,rwork,iwork,info)

       lwork = nint(real(work(1))) 
       deallocate(work)
       allocate(work(lwork))
       lrwork = 5*neig*neig + 7*neig
       deallocate(rwork)
       allocate(rwork(lrwork)) 

       call zgesdd('O',mtxd,neig,psi,mxddim,singval,u,mxddim,            &
     &        vt,neig,work,lwork,rwork,iwork,info)

       deallocate(u)
       deallocate(vt)
       deallocate(work)
       deallocate(iwork)
       deallocate(rwork)

       return
       end subroutine svd_c16



