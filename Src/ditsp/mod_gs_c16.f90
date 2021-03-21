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

!>     Performs a modified gram-schmidt orthogonalization step for
!>     nvec xvec vectors.
!>     If also=.T. performs the corresponding transformation in hxvec

       subroutine mod_gs_c16(also,xvec,hxvec,mtxd,nvec,irow,xscl,        &
     &                     mxddim)

!      Steps are repeated if needed according to the "twice is enough" principle.

!      It also indicates linearly dependent vectors based on xscl.

!      Written January 17 2008. jlm
!      Modified (f90) February 11 2014.
!      Modified real(zdotc(...),REAL64) 9 Agosto 2015. jlm
!      Modified, documentation, January 2020. JLM
!      copyright  Jose Luis Martins/INESC-MN

!      version 4.94

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)


!      input

       integer, intent(in)                ::  mxddim                     !<  array dimension for the hamiltonian

       logical, intent(in)                ::  also                       !<  indicates if hxvec is also transformed
       integer, intent(in)                ::  mtxd                       !<  dimension of the hamiltonian
       integer, intent(in)                ::  nvec                       !<  number of vectors
       real(REAL64), intent(in)           ::  xscl(nvec)                 !<  natural scale (norm) of vector i

!      input and output

       complex(REAL64), intent(inout)     ::  xvec(mxddim,nvec)          !<  component j of vector i
       complex(REAL64), intent(inout)     ::  hxvec(mxddim,nvec)         !<  component j of companion vector i
       integer, intent(inout)             ::  irow(nvec)                 !<  if irow(i)=0 the vector is linearly dependent on the others and the corresponding xvec is zero

!      local variables

       real(REAL64)     ::  suma,xni,xnp,peq
       complex(REAL64)  ::  prod

!      constants

       real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
       real(REAL64), parameter  :: EPS = 1.0E-20

!      counters

       integer    ::  i, j

!      functions

       complex(REAL64)  ::  zdotc
       external  zdotc


       peq = UM/EPS
       peq = peq*tiny(peq)

       suma = real(zdotc(mtxd,xvec(1,1),1,xvec(1,1),1),REAL64)
       
       if(xscl(1) > peq) then
         if( suma/xscl(1) < EPS) then
           suma = ZERO
           irow(1) = 0
         else
           suma = UM/sqrt(suma)
         endif
       else
         suma = ZERO
         irow(1) = 0
       endif

       call zdscal(mtxd,suma,xvec(1,1),1)
       if(also) then
         call zdscal(mtxd,suma,hxvec(1,1),1)
       endif

       if(nvec > 1) then

         do i=2,nvec

           xni = real(zdotc(mtxd,xvec(1,i),1,xvec(1,i),1),REAL64)
           xnp = ZERO

           do j=1,i-1
             prod = zdotc(mtxd,xvec(1,j),1,xvec(1,i),1)
             xnp = xnp + real(prod*conjg(prod),REAL64)
             call zaxpy(mtxd,-prod,xvec(1,j),1,xvec(1,i),1)
             if(also) then
               call zaxpy(mtxd,-prod,hxvec(1,j),1,hxvec(1,i),1)
             endif
           enddo

!          repeat if vector is almost linearly dependent

           if( xnp/xni > 0.9999) then
             do j=1,i-1
               prod = zdotc(mtxd,xvec(1,j),1,xvec(1,i),1)
               call zaxpy(mtxd,-prod,xvec(1,j),1,xvec(1,i),1)
               if(also) then
                 call zaxpy(mtxd,-prod,hxvec(1,j),1,hxvec(1,i),1)
               endif
             enddo
           endif

           suma = real(zdotc(mtxd,xvec(1,i),1,xvec(1,i),1),REAL64)
           if(xscl(i) > peq) then
             if( suma/xscl(i) < EPS) then
               suma = ZERO
               irow(i) = 0
             else
               suma = UM/sqrt(suma)
             endif
           else
             suma = ZERO
             irow(i) = 0
           endif

           call zdscal(mtxd,suma,xvec(1,i),1)
           if(also) then
             call zdscal(mtxd,suma,hxvec(1,i),1)
           endif

         enddo

       endif

       return
       end subroutine mod_gs_c16
