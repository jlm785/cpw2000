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

!>    given the metric of the lattice vectors g_ij = a_i.a_j
!>    it generates the Niggli lattice

      subroutine metric_niggli(iGin,iG,Mtotal)


!     written by Alvaro Ladeira, 1999
!     modified Jose Luis Martins, Fevereiro 2000
!     corrected bugs, jlm, December 2000.
!     modified for integer logic, April 2004
!     modified, f90, 3 June 2014, f90. JLM
!     modified documentation August 2019. JLM
!     copyright Alvaro Ladeira/Jose Luis Martins/INESC-MN

!     send comments/bug reports to jlmartins@inesc-mn.pt

!     version 4.94


      implicit none

      integer, parameter          :: INT64 = selected_int_kind(15)

!     input

      integer(INT64), intent(in)          ::  iGin(3,3)                  !<  input metric

!     output

      integer(INT64), intent(out)         ::  iG(3,3)                    !<  metric of the Niggli lattice
      integer, intent(out)                ::  Mtotal(3,3)                !<  iG=Mtotal iGin (Mtotal)T, det Mtotal=1

!     local variables

      integer         ::  s(3,3), q(3,3), iq(3,3), t(3,3), it(3,3)
      integer         ::  u(3,3), iu(3,3)

      data q  / 0, 0, 1,  1, 0, 0,  0, 1, 0/
      data iq / 0, 1, 0,  0, 0, 1,  1, 0, 0/
      data t  / 0,-1, 0,  1, 0, 0,  0, 0, 1/
      data it / 0, 1, 0, -1, 0, 0,  0, 0, 1/
      data u  / 1, 0, 0,  1, 1, 0,  0, 0, 1/
      data iu / 1, 0, 0, -1, 1, 0,  0, 0, 1/
      data s  / 1, 0, 0,  0, 1, 0,  0, 0,-1/
      
      integer         ::  Prod1(3,3), Prod2(3,3)

      integer         ::  loop1, loop2, nloop
      integer         ::  icount,icount3

      integer         ::  n,l,lmin,lmax,I1,I2

      integer(INT64)  ::  n30,n60
      integer(INT64)  ::  iG_fact(2,3,3)
      integer(INT64)  ::  minor(4,3),det(6)
      integer(INT64)  ::  ip1,ip2,ip3
      
      logical         ::  lx

!     counters

      integer         ::  i, j, k


!     resets Mtotal

      call metric_Tm_int(0,iG,Mtotal)


      n30 = 1
      do i=1,30
        n30 = n30*2
      enddo
      n60 = n30*n30
      nloop = 1000000



      do i=1,3
      do j=1,3
        iG(i,j) = iGin(i,j)
      enddo
      enddo

!dbg      do i=1,3
!dbg        write(6,701) (iG(j,i),j=1,3)
!dbg      enddo
!dbg 701  format(3i20)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C         TEST: Is G symmetric?                                       C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      If(  iG(1,2) /= iG(2,1) .or. iG(1,3) /= iG(3,1) .or.               &
     &     iG(2,3) /= iG(3,2) ) then

         write(6,*) ' metric_niggli'
         write(6,*) 'Metric non simmetric! '

         stop

      EndIf



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C         TEST:  Is there any non positive diagonal element in  G ?   C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      If( (iG(1,1) <= 0) .or. (iG(2,2) <= 0) .or. (iG(3,3) <= 0) ) then

         write(6,*) ' metric_niggli'
         write(6,*) 'Metric diagonal elements must be positive!',        &
     &                iG(1,1),iG(2,2),iG(3,3)

         stop

      EndIf


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C         TEST:  Is there a danger of overflow                       C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      If((iG(1,1) > n60) .or. (iG(2,2) > n60) .or. (iG(3,3) > n60)) then

         write(6,*) ' metric_niggli'
         write(6,*) 'Metric diagonal elements are larger than 2**60!',   &
     &                iG(1,1),iG(2,2),iG(3,3)

         stop

      EndIf

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C         TEST:  Are minors positive definite                         C
!C         Integer algebra should be safe                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       do i=1,3
       do j=1,3
         iG_fact(1,i,j) = iG(i,j)
         iG_fact(2,i,j) = 0
       enddo
       enddo

       do i=1,3
       do j=1,3
           iG_fact(2,i,j) = iG_fact(2,i,j) + iG_fact(1,i,j) / n30
           iG_fact(1,i,j) = mod(iG_fact(1,i,j),n30)
       enddo
       enddo

       do k=1,3
         I1 = mod(k,3) + 1
         I2 = mod(k+1,3) + 1
         do n=1,3
           minor(n,k) = 0
           lmax = min(n,2)
           lmin = max(n-1,1)
           do l=lmin,lmax
             minor(n,k) = minor(n,k)                                     &
     &                    + iG_fact(l,I1,I1)*iG_fact(n-l+1,I2,I2)        &
     &                    - iG_fact(l,I1,I2)*iG_fact(n-l+1,I2,I1)
           enddo
         enddo
         minor(4,k) = 0
       enddo

       do j=1,3
       do k=1,3
         minor(k+1,j) = minor(k+1,j) + minor(k,j) / n30
         minor(k,j) = mod(minor(k,j),n30)
       enddo
       enddo

       do k=1,3
         ip1 = minor(4,k)*n30 + minor(3,k)
         ip2 = minor(2,k)*n30 + minor(1,k)
         if(ip1 < 0 .or. (ip1 == 0 .and. ip2 < 0)) then
           write(6,*) ' metric_niggli'
           write(6,*) ' Minor ',k, 'is negative'

           stop

         endif
       enddo
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C         TEST:  Is the determinant positive                          C
!C         Integer algebra should be safe                              C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       do k=1,3
         I1 = mod(k,3) + 1
         I2 = mod(k+1,3) + 1
         do n=1,3
           minor(n,k) = 0
           lmax = min(n,2)
           lmin = max(n-1,1)
           do l=lmin,lmax
             minor(n,k) = minor(n,k)                                     &
     &                    + iG_fact(l,I1,2)*iG_fact(n-l+1,I2,3)          &
     &                    - iG_fact(l,I1,3)*iG_fact(n-l+1,I2,2)
           enddo
         enddo
         minor(4,k) = 0
       enddo


       do j=1,3
       do k=1,3
         minor(k+1,j) = minor(k+1,j) + minor(k,j) / n30
         minor(k,j) = mod(minor(k,j),n30)
       enddo
       enddo

       do n=1,5
         det(n) = 0
         lmax = min(n,4)
         lmin = max(n-1,1)
         do k=1,3
         do l=lmin,lmax
           det(n) = det(n) + minor(l,k)*iG_fact(n-l+1,k,1)
         enddo
         enddo
       enddo
       det(6) = 0

       do n=1,5
         det(n+1) = det(n+1) + det(n) / n30
         det(n) = mod(det(n),n30)
       enddo
       ip1 = det(6)*n30+det(5)
       ip2 = det(4)*n30+det(3)
       ip3 = det(2)*n30+det(1)
       if(ip1 < 0 .or. (ip1 == 0 .and. (ip2 < 0 .or.                     &
     &                   (ip2 == 0 .and. ip3 <= 0)))) then
         write(6,*) ' metric_niggli'
         write(6,*) ' determinant is non-positive '

         stop

       endif



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C                                                                     C
!C   STEP 1: IT TAKES THE METRIC TO BUERGER'S FORM                     C
!C                                                                     C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


!dbg         write(6,*) '  entering Buerger cell section'
!dbg         write(6,701) (iG(1,i),i=1,3)
!dbg         write(6,701) (iG(2,i),i=1,3)
!dbg         write(6,701) (iG(3,i),i=1,3)

      do loop1=1,20

!       Put into reverse order

        If( iG(1,1) < iG(2,2) ) then
          If( iG(2,2) < iG(3,3) ) then
!               Gii->(X,Y,Z)                 tm = iq.t
            call metric_Tm_int(1,iG,t)
            call metric_Tm_int(1,iG,iq)
          Else
            If( iG(1,1) <= iG(3,3) ) then
!               Gii->(X,Z,Y),(X,Y,Y),(X,Z,X) tm = q
              call metric_Tm_int(1,iG,q)
            Else
!               Gii->(Y,Z,X)                 tm = t
              call metric_Tm_int(1,iG,t)
            EndIf
          EndIf
        Else
          If( (iG(1,1) <= iG(3,3)) .and. (iG(2,2) /= iG(3,3)) ) then
!               Gii->(Y,X,Z),(Y,X,Y),(X,X,Z) tm = iq
            call metric_Tm_int(1,iG,iq)
          Else
            If( iG(2,2) < iG(3,3) ) then
!               Gii->(Z,X,Y)                 tm = q.t
              call metric_Tm_int(1,iG,t)
              call metric_Tm_int(1,iG,q)
            EndIf
          EndIf
        EndIf

!CC O que sai da estrutura IF-ELSE-ENDIF anterior e' uma metrica numa
!     das formas  Gii-> (Z,Y,X), (Y,Y,X), (Z,X,X) ou (X,X,X)
!CC E que tal comecar a caminhar na direccao Niggli, ate a celula de
! Buerger?...

        icount = 0
        icount3 = 0

        do loop2=1,nloop
          If( 2*iG(1,2) + iG(2,2) < 0) then
!                tm = u
            call metric_Tm_int(1,iG,u)
            icount = 1
          ElseIf( 2*iG(1,2) - iG(2,2) > 0) then
!               tm = iu
            call metric_Tm_int(1,iG,iu)
            icount = 1
          Else

            exit

          EndIf
        enddo


        do loop2=1,nloop
          If( 2*iG(1,3) + iG(3,3) < 0) then
!                 tm = iu.q.t
            call metric_Tm_int(1,iG,t)
            call metric_Tm_int(1,iG,q)
            call metric_Tm_int(1,iG,iu)
            call metric_Tm_int(1,iG,iq)
            call metric_Tm_int(1,iG,it)
            icount = 2
          ElseIf( 2*iG(1,3) - iG(3,3) > 0) then
!                  tm = u.q.t
            call metric_Tm_int(1,iG,t)
            call metric_Tm_int(1,iG,q)
            call metric_Tm_int(1,iG,u)
            call metric_Tm_int(1,iG,iq)
            call metric_Tm_int(1,iG,it)
            icount = 2
          Else

            exit

          EndIf
        enddo



        do loop2=1,nloop
          If( 2*iG(2,3) + iG(3,3) < 0) then
!                    tm = u.q
            call metric_Tm_int(1,iG,q)
            call metric_Tm_int(1,iG,u)
            call metric_Tm_int(1,iG,iq)
            icount = 3
          ElseIf( 2*iG(2,3) - iG(3,3) > 0) then
!                    tm = iu.q
            call metric_Tm_int(1,iG,q)
            call metric_Tm_int(1,iG,iu)
            call metric_Tm_int(1,iG,iq)
            icount = 3
          Else

            exit

          EndIf
        enddo


        If(icount == 0) then

          If(iG(2,2)+iG(3,3)                                             & 
     &          + 2*( iG(1,2)+iG(1,3)+iG(2,3)) < 0) then
            call metric_Tm_int(1,iG,u)
            call metric_Tm_int(1,iG,t)
            call metric_Tm_int(1,iG,q)
            call metric_Tm_int(1,iG,iu)
            icount3 = 1
          ElseIf(iG(2,2)+iG(3,3)                                         &
     &          + 2*( iG(1,2)-iG(1,3)-iG(2,3)) < 0) then
            call metric_Tm_int(1,iG,u)
            call metric_Tm_int(1,iG,t)
            call metric_Tm_int(1,iG,q)
            call metric_Tm_int(1,iG,u)
            icount3 = 2
          ElseIf(iG(2,2)+iG(3,3)                                         &
     &          + 2*(-iG(1,2)+iG(1,3)-iG(2,3)) < 0) then
            call metric_Tm_int(1,iG,iu)
            call metric_Tm_int(1,iG,t)
            call metric_Tm_int(1,iG,q)
            call metric_Tm_int(1,iG,iu)
            icount3 = 3
          ElseIf(iG(2,2)+iG(3,3)                                         &
     &          + 2*(-iG(1,2)-iG(1,3)+iG(2,3)) < 0) then
            call metric_Tm_int(1,iG,iu)
            call metric_Tm_int(1,iG,t)
            call metric_Tm_int(1,iG,q)
            call metric_Tm_int(1,iG,u)
            icount3 = 4
          EndIf

        EndIf

        lx = .TRUE.
        If(icount == 0 .and. icount3 == 0) then

!dbg      write(6,*) 'Chegamos a celula de Buerger!...'

          If( iG(1,1) > iG(2,2) ) then
            If( iG(2,2) > iG(3,3) ) then
!                tm = t.q
              call metric_Tm_int(1,iG,q)
              call metric_Tm_int(1,iG,t)
            Else
!                tm = q
              call metric_Tm_int(1,iG,q)
            EndIf
          ElseIf( iG(2,2) > iG(3,3) ) then
!              tm = iq
            call metric_Tm_int(1,iG,iq)
          EndIf

          lx = .FALSE.
          exit

        endif

      enddo

      if(lx) then
        write(6,*) ' metric_ident'
        write(6,*) ' loop sizes were too small to reach Buerger cell'

        stop
      
      endif


!dbg         write(6,*) '  arrived at Buerger cell'
!dbg         write(6,701) (iG(1,i),i=1,3)
!dbg         write(6,701) (iG(2,i),i=1,3)
!dbg         write(6,701) (iG(3,i),i=1,3)


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C                                                                     C
!C   STEP 2: IT TAKES THE METRIC FROM BUERGER'S TO NIGGLI'S FORM       C
!C                                                                     C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


!CCCCCCCC   It brings the metric to defined + or -  CCCCCCCCCCCC


      If(  ((iG(1,2) == 0) .and. (iG(1,3) > 0) .and. (iG(2,3) < 0))      &
     & .or.((iG(1,2) > 0) .and. (iG(1,3) >= 0) .and. (iG(2,3) <= 0))     &
     & .or.((iG(1,2) < 0) .and. (iG(1,3) < 0) .and. (iG(2,3) > 0)) )then
!      Apply  TM = s.iq.t.t.q.t.t =  -1  0  0
!                                     0  1  0
!                                     0  0  1
        call metric_dotM(Prod1, t,t)
        call metric_dotM(Prod2, q,Prod1)
        call metric_dotM(Prod1, t,Prod2)
        call metric_dotM(Prod2, t,Prod1)
        call metric_dotM(Prod1,iq,Prod2)
        call metric_dotM(Prod2, s,Prod1)
        call metric_Tm_int(1,iG,Prod2)
      ElseIf(((iG(1,2) >= 0) .and. (iG(1,3) < 0) .and. (iG(2,3) > 0))    &
     &   .or.((iG(1,2) > 0) .and. (iG(1,3) == 0) .and. (iG(2,3) > 0))    &
     &   .or.((iG(1,2) > 0) .and. (iG(1,3) < 0) .and. (iG(2,3) == 0))    &
     &   .or.((iG(1,2) < 0) .and. (iG(1,3) > 0) .and. (iG(2,3) < 0))     &
     &                                                           )then
!        Apply  TM = s.iq.t.iq.it =   1  0  0
!                                     0 -1  0
!                                     0  0  1
        call metric_dotM(Prod1,iq,it)
        call metric_dotM(Prod2, t,Prod1)
        call metric_dotM(Prod1,iq,Prod2)
        call metric_dotM(Prod2, s,Prod1)
        call metric_Tm_int(1,iG,Prod2)
      ElseIf(((iG(1,2) <= 0) .and. (iG(1,3) > 0) .and. (iG(2,3) == 0))   &
     &   .or.((iG(1,2) <= 0) .and. (iG(1,3) >= 0) .and. (iG(2,3) > 0))   &
     &   .or.((iG(1,2) > 0) .and. (iG(1,3) < 0) .and. (iG(2,3) < 0))     &
     &                                                           )then
!                   Apply  TM = s =   1  0  0
!                                     0  1  0
!                                     0  0 -1
        call metric_Tm_int(1,iG,s)
      EndIf



!     X   X   X   X   X   X
!      X X     X X     X X
!       X       X       X
!      X X     X X     X X
!     X   X   X   X   X   X


      If( (iG(1,1) == iG(2,2)) .and. (iG(2,2) == iG(3,3)) ) then
! =================== POSITIVE REDUCED =================================
        If( iG(1,2) > 0 .and. iG(1,3) > 0 .and. iG(2,3) > 0 .or.         &
     &      iG(1,2) > 0 .and. iG(1,3) < 0 .and. iG(2,3) < 0 .or.         &
     &      iG(1,2) < 0 .and. iG(1,3) > 0 .and. iG(2,3) < 0 .or.         &
     &      iG(1,2) < 0 .and. iG(1,3) < 0 .and. iG(2,3) > 0 ) then

!              ordena G's ficando S>=T>=U:
          If( (iG(1,2) >= iG(2,3)) .and. (iG(2,3) > iG(1,3)) ) then
!           Apply  TM = s.iq.t.iq =  0 1 0
!                                    1 0 0
!                                    0 0 1
            call metric_dotM(Prod1, t,iq)
            call metric_dotM(Prod2,iq,Prod1)
            call metric_dotM(Prod1, s,Prod2)
            call metric_Tm_int(1,iG,Prod1)
          ElseIf( (iG(2,3) > iG(1,2)) .and. (iG(1,2) >= iG(1,3)) ) then
!                   Apply  TM = q =  0 1 0
!                                    0 0 1
!                                    1 0 0
            call metric_Tm_int(1,iG, q)
          ElseIf( (iG(1,3) > iG(1,2)) .and. (iG(1,2) > iG(1,3)) ) then
!            Apply  TM = s.iq.t.q =  1 0 0
!                                    0 0 1
!                                    0 1 0
            call metric_dotM(Prod1, t, q)
            call metric_dotM(Prod2,iq,Prod1)
            call metric_dotM(Prod1, s,Prod2)
            call metric_Tm_int(1,iG,Prod1)
          ElseIf( iG(2,3) > iG(1,3) ) then
!              Apply  TM = s.iq.t =  0 0 1
!                                    0 1 0
!                                    1 0 0
            call metric_dotM(Prod1,iq, t)
            call metric_dotM(Prod2, s,Prod1)
            call metric_Tm_int(1,iG,Prod2)
          ElseIf( iG(1,2) < iG(1,3) ) then
!                  Apply  TM = iq =  0 0 1
!                                    1 0 0
!                                    0 1 0
            call metric_Tm_int(1,iG,iq)
          EndIf

          If( 2*iG(1,2) == iG(1,1) ) then
            If( 2*iG(1,3) == iG(1,1) ) then
              If( 2*iG(2,3) /= iG(2,2) .and. 2*iG(2,3) < iG(1,3) ) then
!dbg               write(6,*) 'Apply  t15 = s.t.iu.q.it.q'
                call metric_dotM(Prod1,it, q)
                call metric_dotM(Prod2, q,Prod1)
                call metric_dotM(Prod1,iu,Prod2)
                call metric_dotM(Prod2, t,Prod1)
                call metric_dotM(Prod1, s,Prod2)
                call metric_Tm_int(1,iG,Prod1)
              EndIf
            Else
              If( 2*iG(2,3) < iG(1,3) ) then
!dbg               write(6,*) 'Apply  t15 = s.t.iu.q.it.q'
                call metric_dotM(Prod1,it, q)
                call metric_dotM(Prod2, q,Prod1)
                call metric_dotM(Prod1,iu,Prod2)
                call metric_dotM(Prod2, t,Prod1)
                call metric_dotM(Prod1, s,Prod2)
                call metric_Tm_int(1,iG,Prod1)
              EndIf
            EndIf
          EndIf
! =================== NEGATIVE REDUCED =================================
        Else
!            ordena G's ficando  S >= T >= U :
          If( (iG(1,2) <= iG(2,3)) .and. (iG(2,3) < iG(1,3)) ) then
!           Apply  TM = s.iq.t.iq =  0 1 0
!                                    1 0 0
!                                    0 0 1
            call metric_dotM(Prod1, t,iq)
            call metric_dotM(Prod2,iq,Prod1)
            call metric_dotM(Prod1, s,Prod2)
            call metric_Tm_int(1,iG,Prod1)
          ElseIf( (iG(2,3) < iG(1,2)) .and. (iG(1,2) <= iG(1,3)) ) then
!                   Apply  TM = q =  0 1 0
!                                    0 0 1
!                                    1 0 0
            call metric_Tm_int(1,iG, q)
          ElseIf( (iG(1,3) < iG(1,2)) .and. (iG(1,2) < iG(1,3)) ) then
!            Apply  TM = s.iq.t.q =  1 0 0
!                                    0 0 1
!                                    0 1 0
            call metric_dotM(Prod1, t, q)
            call metric_dotM(Prod2,iq,Prod1)
            call metric_dotM(Prod1, s,Prod2)
            call metric_Tm_int(1,iG,Prod1)
          ElseIf( iG(2,3) < iG(1,3) ) then
!              Apply  TM = s.iq.t =  0 0 1
!                                    0 1 0
!                                    1 0 0
            call metric_dotM(Prod1,iq, t)
            call metric_dotM(Prod2, s,Prod1)
            call metric_Tm_int(1,iG,Prod2)
          ElseIf( iG(1,2) > iG(1,3) ) then
!                  Apply  TM = iq =  0 0 1
!                                    1 0 0
!                                    0 1 0
            call metric_Tm_int(1,iG,iq)
          EndIf
! ------------   SUM < (more negative) -X   ----------------------
          If( iG(1,2)+iG(1,3)+iG(2,3) < -iG(1,1) ) then
            If( iG(1,2) == iG(1,3) ) then
              If( iG(1,3) == iG(2,3) ) then
!dbg               write(6,*) 'Apply  t17 = iu.q.t.u.s.q.it.q'
                call metric_dotM(Prod1,it, q)
                call metric_dotM(Prod2, q,Prod1)
                call metric_dotM(Prod1, s,Prod2)
                call metric_dotM(Prod2, u,Prod1)
                call metric_dotM(Prod1, t,Prod2)
                call metric_dotM(Prod2, q,Prod1)
                call metric_dotM(Prod1,iu,Prod2)
                call metric_Tm_int(1,iG,Prod1)
              Else
                If( 2*iG(1,2) == -iG(1,1) ) then
!dbg                 write(6,*) 'Apply  t19 = u.iq.u.t.t.s'
                  call metric_dotM(Prod1, t, s)
                  call metric_dotM(Prod2, t,Prod1)
                  call metric_dotM(Prod1, u,Prod2)
                  call metric_dotM(Prod2,iq,Prod1)
                  call metric_dotM(Prod1, s,Prod2)
                  call metric_Tm_int(1,iG,Prod1)
                Else
!dbg                 write(6,*) 'Apply  t17 = iu.q.t.u.s.q.it.q'
                  call metric_dotM(Prod1,it, q)
                  call metric_dotM(Prod2, q,Prod1)
                  call metric_dotM(Prod1, s,Prod2)
                  call metric_dotM(Prod2, u,Prod1)
                  call metric_dotM(Prod1, t,Prod2)
                  call metric_dotM(Prod2, q,Prod1)
                  call metric_dotM(Prod1,iu,Prod2)
                  call metric_Tm_int(1,iG,Prod1)
                EndIf
              EndIf
            ElseIf( iG(1,3) == iG(2,3) ) then
              If( 2*iG(1,2) == -iG(1,1) ) then
!dbg               write(6,*) 'Apply  t19 = u.iq.u.t.t.s'
                call metric_dotM(Prod1, t, s)
                call metric_dotM(Prod2, t,Prod1)
                call metric_dotM(Prod1, u,Prod2)
                call metric_dotM(Prod2,iq,Prod1)
                call metric_dotM(Prod1, s,Prod2)
                call metric_Tm_int(1,iG,Prod1)
              Else
!dbg               write(6,*) 'Apply  t17 = iu.q.t.u.s.q.it.q'
                call metric_dotM(Prod1,it, q)
                call metric_dotM(Prod2, q,Prod1)
                call metric_dotM(Prod1, s,Prod2)
                call metric_dotM(Prod2, u,Prod1)
                call metric_dotM(Prod1, t,Prod2)
                call metric_dotM(Prod2, q,Prod1)
                call metric_dotM(Prod1,iu,Prod2)
                call metric_Tm_int(1,iG,Prod1)
              EndIf
            Else
              If( 2*iG(1,2) == -iG(1,1) ) then
!dbg               write(6,*) 'Apply  t19 = u.iq.u.t.t.s'
                call metric_dotM(Prod1, t, s)
                call metric_dotM(Prod2, t,Prod1)
                call metric_dotM(Prod1, u,Prod2)
                call metric_dotM(Prod2,iq,Prod1)
                call metric_dotM(Prod1, s,Prod2)
                call metric_Tm_int(1,iG,Prod1)
              Else
!dbg               write(6,*) 'Apply  t20 = iu.it.iq.u'
                call metric_dotM(Prod1,iq, u)
                call metric_dotM(Prod2,it,Prod1)
                call metric_dotM(Prod1,iu,Prod2)
                call metric_Tm_int(1,iG,Prod1)
              EndIf
            EndIf
! -----------  SUM = -X   ---------------------------------------
          ElseIf( iG(1,2)+iG(1,3)+iG(2,3) == -iG(1,1) ) then
            If( 2*iG(1,2) == -iG(1,1) ) then
              If( 2*iG(1,3) == -iG(1,1) ) then
!dbg               write(6,*) 'Apply  t12 = s.t.iu.it'
                call metric_dotM(Prod1,iu,it)
                call metric_dotM(Prod2, t,Prod1)
                call metric_dotM(Prod1, s,Prod2)
                call metric_Tm_int(1,iG,Prod1)
              ElseIf( iG(1,3) == iG(2,3) ) then
!dbg               write(6,*) 'Apply  t18 = s.u'
                call metric_dotM(Prod1,s, u)
                call metric_Tm_int(1,iG,Prod1)
              Else
!dbg               write(6,*) 'Apply  t13 = q.s.t.iu.it'
                call metric_dotM(Prod1,iu,it)
                call metric_dotM(Prod2, t,Prod1)
                call metric_dotM(Prod1, s,Prod2)
                call metric_dotM(Prod2, q,Prod1)
                call metric_Tm_int(1,iG,Prod2)
              EndIf
            EndIf
! ------   SUM > (menos negativo) -X   ------------------------
          Else
            If( 2*iG(1,2) == -iG(1,1) ) then
              If( iG(2,3) == 0 ) then
                If( iG(1,3) /= 0 ) then
!dbg                 write(6,*) 'Apply  t12 = s.t.iu.it'
                  call metric_dotM(Prod1,iu,it)
                  call metric_dotM(Prod2, t,Prod1)
                  call metric_dotM(Prod1, s,Prod2)
                  call metric_Tm_int(1,iG,Prod1)
                EndIf
              Else
                If( iG(1,3) == iG(2,3) ) then
!dbg                 write(6,*) 'Apply  t18 = s.u'
                  call metric_dotM(Prod1,s, u)
                  call metric_Tm_int(1,iG,Prod1)
                Else
!dbg                 write(6,*) 'Apply  t14 = iq.t.iq.t.iu.it'
                  call metric_dotM(Prod1,iu,it)
                  call metric_dotM(Prod2, t,Prod1)
                  call metric_dotM(Prod1,iq,Prod2)
                  call metric_dotM(Prod2, t,Prod1)
                  call metric_dotM(Prod1,iq,Prod2)
                  call metric_Tm_int(1,iG,Prod1)
                EndIf
              EndIf
            EndIf
          EndIf
        EndIf

!    X   X   X   X   XXXXX
!     X X     X X       X
!      X       X       X
!     X X     X X     X
!    X   X   X   X   XXXXX

      Else
        If( (iG(1,1) == iG(2,2)) ) then
! =================== POSITIVE REDUCED =================================
          If( iG(1,2) > 0 .and. iG(1,3) > 0 .and. iG(2,3) > 0 .or.       &
     &        iG(1,2) > 0 .and. iG(1,3) < 0 .and. iG(2,3) < 0 .or.       &
     &        iG(1,2) < 0 .and. iG(1,3) > 0 .and. iG(2,3) < 0 .or.       &
     &        iG(1,2) < 0 .and. iG(1,3) < 0 .and. iG(2,3) > 0 ) then
!                ordena G's ficando T>=U:
            If( iG(2,3) > iG(1,3) ) then
!           Apply  TM = s.iq.t.iq =  0 1 0
!                                    1 0 0
!                                    0 0 1
              call metric_dotM(Prod1, t,iq)
              call metric_dotM(Prod2,iq,Prod1)
              call metric_dotM(Prod1, s,Prod2)
              call metric_Tm_int(1,iG,Prod1)
            EndIf

            If( 2*iG(1,2) == iG(1,1) ) then
              If( 2*iG(1,3) == iG(1,1) ) then
                If(2*iG(2,3) /= iG(2,2) .and. 2*iG(2,3) < iG(1,3)) then
!dbg                 write(6,*) 'Apply  t2 = s.q.iu.iq'
                  call metric_dotM(Prod1,iu,iq)
                  call metric_dotM(Prod2, q,Prod1)
                  call metric_dotM(Prod1, s,Prod2)
                  call metric_Tm_int(1,iG,Prod1)
                EndIf
              Else
                If( 2*iG(2,3) < iG(1,3) ) then
!dbg                 write(6,*) 'Apply  t15 = s.t.iu.q.it.q'
                  call metric_dotM(Prod1,it, q)
                  call metric_dotM(Prod2, q,Prod1)
                  call metric_dotM(Prod1,iu,Prod2)
                  call metric_dotM(Prod2, t,Prod1)
                  call metric_dotM(Prod1, s,Prod2)
                  call metric_Tm_int(1,iG,Prod1)
                EndIf
              EndIf
            Else
              If( 2*iG(1,3) == iG(1,1) ) then
                If( 2*iG(2,3) < iG(1,2) ) then
!dbg                 write(6,*) 'Apply  t2 = s.q.iu.iq'
                  call metric_dotM(Prod1,iu,iq)
                  call metric_dotM(Prod2, q,Prod1)
                  call metric_dotM(Prod1, s,Prod2)
                  call metric_Tm_int(1,iG,Prod1)
                EndIf
              EndIf
            EndIf
! =================== NEGATIVE REDUCED =================================
          Else
!                ordena G's ficando  T >= U :
            If( iG(2,3) < iG(1,3) ) then
!           Apply  TM = s.iq.t.iq =  0 1 0
!                                    1 0 0
!                                    0 0 1
              call metric_dotM(Prod1, t,iq)
              call metric_dotM(Prod2,iq,Prod1)
              call metric_dotM(Prod1, s,Prod2)
              call metric_Tm_int(1,iG,Prod1)
            EndIf
! ------------   SUM < (more negative) -X   -----------------------
            If( iG(1,2)+iG(1,3)+iG(2,3) < -iG(1,1) ) then
              If( 2*iG(1,2) == -iG(1,1) ) then
                If( 2*iG(1,3) == -iG(1,1) ) then
                  If( 2*iG(2,3) == -iG(2,2) ) then
!dbg                   write(6,*) 'Apply  t3 = s.q.u.it.iq.t.t.iu.iq.it'
                    call metric_dotM(Prod1,iq,it)
                    call metric_dotM(Prod2,iu,Prod1)
                    call metric_dotM(Prod1, t,Prod2)
                    call metric_dotM(Prod2, t,Prod1)
                    call metric_dotM(Prod1,iq,Prod2)
                    call metric_dotM(Prod2,it,Prod1)
                    call metric_dotM(Prod1, u,Prod2)
                    call metric_dotM(Prod2, q,Prod1)
                    call metric_dotM(Prod1, s,Prod2)
                    call metric_Tm_int(1,iG,Prod1)
                  Else
!dbg                   write(6,*) 'Apply  t5 = q.iu.q.t.u.iq'
                    call metric_dotM(Prod1, u,iq)
                    call metric_dotM(Prod2, t,Prod1)
                    call metric_dotM(Prod1, q,Prod2)
                    call metric_dotM(Prod2,iu,Prod1)
                    call metric_dotM(Prod1, q,Prod2)
                    call metric_Tm_int(1,iG,Prod1)
                  EndIf
                Else
!dbg                 write(6,*) 'Apply  t10 = s.q.t.t.iu.iq.iu.it'
                  call metric_dotM(Prod1,iu,it)
                  call metric_dotM(Prod2,iq,Prod1)
                  call metric_dotM(Prod1,iu,Prod2)
                  call metric_dotM(Prod2, t,Prod1)
                  call metric_dotM(Prod1, t,Prod2)
                  call metric_dotM(Prod2, q,Prod1)
                  call metric_dotM(Prod1, s,Prod2)
                  call metric_Tm_int(1,iG,Prod1)
                EndIf
              Else
                If( 2*iG(1,3) == -iG(1,1) ) then
                  If( 2*iG(2,3) == -iG(2,2) ) then
!dbg                   write(6,*) 'Apply  t3 = s.q.u.it.iq.t.t.iu.iq.it'
                    call metric_dotM(Prod1,iq,it)
                    call metric_dotM(Prod2,iu,Prod1)
                    call metric_dotM(Prod1, t,Prod2)
                    call metric_dotM(Prod2, t,Prod1)
                    call metric_dotM(Prod1,iq,Prod2)
                    call metric_dotM(Prod2,it,Prod1)
                    call metric_dotM(Prod1, u,Prod2)
                    call metric_dotM(Prod2, q,Prod1)
                    call metric_dotM(Prod1, s,Prod2)
                    call metric_Tm_int(1,iG,Prod1)
                  Else
!dbg                   write(6,*) 'Apply  t5 = q.iu.q.t.u.iq'
                    call metric_dotM(Prod1, u,iq)
                    call metric_dotM(Prod2, t,Prod1)
                    call metric_dotM(Prod1, q,Prod2)
                    call metric_dotM(Prod2,iu,Prod1)
                    call metric_dotM(Prod1, q,Prod2)
                    call metric_Tm_int(1,iG,Prod1)
                  EndIf
                Else
!dbg                 write(6,*) 'Apply  t5 = q.iu.q.t.u.iq'
                  call metric_dotM(Prod1, u,iq)
                  call metric_dotM(Prod2, t,Prod1)
                  call metric_dotM(Prod1, q,Prod2)
                  call metric_dotM(Prod2,iu,Prod1)
                  call metric_dotM(Prod1, q,Prod2)
                  call metric_Tm_int(1,iG,Prod1)
                EndIf
              EndIf
! ------------   SUM = -X   ----------------------------------------
            ElseIf( iG(1,2)+iG(1,3)+iG(2,3) == -iG(1,1) ) then
              If( 2*iG(1,2) == -iG(1,1) ) then
                If( 2*iG(1,3) == -iG(1,1) ) then
!dbg                 write(6,*) 'Apply  t6 = q.u.iq.t.q.it.q.s'
                  call metric_dotM(Prod1, u,iq)
                  call metric_dotM(Prod2, t,Prod1)
                  call metric_dotM(Prod1, q,Prod2)
                  call metric_dotM(Prod2,iu,Prod1)
                  call metric_dotM(Prod1, q,Prod2)
                  call metric_Tm_int(1,iG,Prod1)
                Else
!dbg                 write(6,*) 'Apply  t14 = iq.t.iq.t.iu.it'
                  call metric_dotM(Prod1,iu,it)
                  call metric_dotM(Prod2, t,Prod1)
                  call metric_dotM(Prod1,iq,Prod2)
                  call metric_dotM(Prod2, t,Prod1)
                  call metric_dotM(Prod1,iq,Prod2)
                  call metric_Tm_int(1,iG,Prod1)
                EndIf
              Else
                If( 2*iG(1,3) == -iG(1,1)) then
                  If( 2*iG(2,3) /= -iG(2,2)) then
!dbg                   write(6,*) 'Apply  t6 = q.u.iq.t.q.it.q.s'
                    call metric_dotM(Prod1, u,iq)
                    call metric_dotM(Prod2, t,Prod1)
                    call metric_dotM(Prod1, q,Prod2)
                    call metric_dotM(Prod2,iu,Prod1)
                    call metric_dotM(Prod1, q,Prod2)
                    call metric_Tm_int(1,iG,Prod1)
                  EndIf
                EndIf
              EndIf
! ------------   SUM > (less negative) -X   -------------------------
            Else
              If( 2*iG(1,2) == -iG(1,1) ) then
                If( iG(2,3) == 0 ) then
                  If( iG(1,3) /= 0 ) then
!dbg                   write(6,*) 'Apply  t12 = s.t.iu.it'
                    call metric_dotM(Prod1,iu,it)
                    call metric_dotM(Prod2, t,Prod1)
                    call metric_dotM(Prod1, s,Prod2)
                    call metric_Tm_int(1,iG,Prod1)
                  EndIf
                Else
!dbg                 write(6,*) 'Apply  t14 = iq.t.iq.t.iu.it'
                  call metric_dotM(Prod1,iu,it)
                  call metric_dotM(Prod2, t,Prod1)
                  call metric_dotM(Prod1,iq,Prod2)
                  call metric_dotM(Prod2, t,Prod1)
                  call metric_dotM(Prod1,iq,Prod2)
                  call metric_Tm_int(1,iG,Prod1)
                EndIf
              Else
                If( 2*iG(1,3) == -iG(1,1) ) then
                  If( iG(1,2) /= 0 ) then
!dbg                   write(6,*) 'Apply  t6 = q.u.iq.t.q.it.q.s'
                    call metric_dotM(Prod1, u,iq)
                    call metric_dotM(Prod2, t,Prod1)
                    call metric_dotM(Prod1, q,Prod2)
                    call metric_dotM(Prod2,iu,Prod1)
                    call metric_dotM(Prod1, q,Prod2)
                    call metric_Tm_int(1,iG,Prod1)
                  EndIf
                EndIf
              EndIf
            EndIf
          EndIf

!    X   X   X   X   X   X
!     X X     X X     X X
!      X       X       X
!     X X      X       X
!    X   X     X       X

        Else
          If( (iG(2,2) == iG(3,3)) ) then
! =================== POSITIVE REDUCED =================================
            If( iG(1,2) > 0 .and. iG(1,3) > 0 .and. iG(2,3) > 0 .or.     &
     &          iG(1,2) > 0 .and. iG(1,3) < 0 .and. iG(2,3) < 0 .or.     &
     &          iG(1,2) < 0 .and. iG(1,3) > 0 .and. iG(2,3) < 0 .or.     &
     &          iG(1,2) < 0 .and. iG(1,3) < 0 .and. iG(2,3) > 0 )then
              If( iG(1,3) > iG(1,2) ) then
!               Apply  TM = s.iq.t.q =  1 0 0
!                                       0 0 1
!                                       0 1 0
                call metric_dotM(Prod1, t, q)
                call metric_dotM(Prod2,iq,Prod1)
                call metric_dotM(Prod1, s,Prod2)
                call metric_Tm_int(1,iG,Prod1)
              EndIf
              If( 2*iG(2,3) == iG(2,2) ) then
                If(2*iG(1,3) /= iG(1,1) .and. 2*iG(2,3) < iG(1,3)) then
!dbg                 write(6,*) 'Apply  t15 = s.t.iu.q.it.q'
                  call metric_dotM(Prod1,it, q)
                  call metric_dotM(Prod2, q,Prod1)
                  call metric_dotM(Prod1,iu,Prod2)
                  call metric_dotM(Prod2, t,Prod1)
                  call metric_dotM(Prod1, s,Prod2)
                  call metric_Tm_int(1,iG,Prod1)
                Else
                  If( 2*iG(1,2) == iG(1,1) ) then
                    If( 2*iG(2,3) < iG(1,2) ) then
!dbg                     write(6,*) 'Apply  t15 = s.t.iu.q.it.q'
                      call metric_dotM(Prod1,it, q)
                      call metric_dotM(Prod2, q,Prod1)
                      call metric_dotM(Prod1,iu,Prod2)
                      call metric_dotM(Prod2, t,Prod1)
                      call metric_dotM(Prod1, s,Prod2)
                      call metric_Tm_int(1,iG,Prod1)
                    EndIf
                  EndIf
                EndIf
              Else
                If( 2*iG(1,2) == iG(1,1) ) then
                  If( 2*iG(1,3) /= iG(1,1) ) then
                    If( 2*iG(1,3) < iG(1,2) ) then
!dbg                     write(6,*) 'Apply  t1 = s.q.it.iq.t.t.u.iq.it'
                      call metric_dotM(Prod1,iq,it)
                      call metric_dotM(Prod2, u,Prod1)
                      call metric_dotM(Prod1, t,Prod2)
                      call metric_dotM(Prod2, t,Prod1)
                      call metric_dotM(Prod1,iq,Prod2)
                      call metric_dotM(Prod2,it,Prod1)
                      call metric_dotM(Prod1, q,Prod2)
                      call metric_dotM(Prod2, s,Prod1)
                      call metric_Tm_int(1,iG,Prod2)
                    EndIf
                  EndIf
                Else
                  If( 2*iG(1,3)  <  iG(1,2) ) then
!dbg                   write(6,*) 'Apply  t1 = s.q.it.iq.t.t.u.iq.it'
                    call metric_dotM(Prod1,iq,it)
                    call metric_dotM(Prod2, u,Prod1)
                    call metric_dotM(Prod1, t,Prod2)
                    call metric_dotM(Prod2, t,Prod1)
                    call metric_dotM(Prod1,iq,Prod2)
                    call metric_dotM(Prod2,it,Prod1)
                    call metric_dotM(Prod1, q,Prod2)
                    call metric_dotM(Prod2, s,Prod1)
                    call metric_Tm_int(1,iG,Prod2)
                  EndIf
                EndIf
              EndIf
! =================== NEGATIVE REDUCED =================================
            Else
              If( iG(1,3) < iG(1,2) ) then
!               Apply  TM = s.iq.t.q =  1 0 0
!                                       0 0 1
!                                       0 1 0
                call metric_dotM(Prod1, t, q)
                call metric_dotM(Prod2,iq,Prod1)
                call metric_dotM(Prod1, s,Prod2)
                call metric_Tm_int(1,iG,Prod1)
              EndIf
! ------------   SUM < (more negative) -(X+Y)/2   --------------------
              If( 2*(iG(1,2)+iG(1,3)+iG(2,3)) <  -(iG(1,1)+iG(2,2)) )    &
     &                                                           then
                If( 2*iG(2,3) == -iG(2,2) ) then
                  If( 2*iG(1,2) == -iG(1,1) ) then
                    If( 2*iG(1,3) == -iG(1,1) ) then
!dbg                     write(6,*) 'Apply  t4 = iq.iu.q.t.u.iq'
                      call metric_dotM(Prod1, u,iq)
                      call metric_dotM(Prod2, u,Prod1)
                      call metric_dotM(Prod1, t,Prod2)
                      call metric_dotM(Prod2, q,Prod1)
                      call metric_dotM(Prod1,iu,Prod2)
                      call metric_dotM(Prod2,iq,Prod1)
                      call metric_Tm_int(1,iG,Prod2)
                    Else
!dbg                     write(6,*) 'Apply  t16 = iq.iu.q.t.u.s.q.it.q'
                      call metric_dotM(Prod1,it, q)
                      call metric_dotM(Prod2, q,Prod1)
                      call metric_dotM(Prod1, s,Prod2)
                      call metric_dotM(Prod2, u,Prod1)
                      call metric_dotM(Prod1, t,Prod2)
                      call metric_dotM(Prod2, q,Prod1)
                      call metric_dotM(Prod1,iu,Prod2)
                      call metric_dotM(Prod2,iq,Prod1)
                      call metric_Tm_int(1,iG,Prod2)
                    EndIf
                  Else
!dbg                   write(6,*) 'Apply  t4 = iq.iu.q.t.u.iq'
                    call metric_dotM(Prod1, u,iq)
                    call metric_dotM(Prod2, u,Prod1)
                    call metric_dotM(Prod1, t,Prod2)
                    call metric_dotM(Prod2, q,Prod1)
                    call metric_dotM(Prod1,iu,Prod2)
                    call metric_dotM(Prod2,iq,Prod1)
                    call metric_Tm_int(1,iG,Prod2)
                  EndIf
                Else
                  If( 2*iG(1,2) == -iG(1,1) ) then
                    If( 2*iG(1,3) == -iG(1,1) ) then
!dbg                     write(6,*) 'Apply  t4 = iq.iu.q.t.u.iq'
                      call metric_dotM(Prod1, u,iq)
                      call metric_dotM(Prod2, u,Prod1)
                      call metric_dotM(Prod1, t,Prod2)
                      call metric_dotM(Prod2, q,Prod1)
                      call metric_dotM(Prod1,iu,Prod2)
                      call metric_dotM(Prod2,iq,Prod1)
                      call metric_Tm_int(1,iG,Prod2)
                    Else
!dbg                     write(6,*) 'Apply  t11 = s.iq.iu.iq.iu.it'
                      call metric_dotM(Prod1,iu,it)
                      call metric_dotM(Prod2,iq,Prod1)
                      call metric_dotM(Prod1,iu,Prod2)
                      call metric_dotM(Prod2,iq,Prod1)
                      call metric_dotM(Prod1, s,Prod2)
                      call metric_Tm_int(1,iG,Prod1)
                    EndIf
                  Else
!dbg                   write(6,*) 'Apply  t4 = iq.iu.q.t.u.iq'
                    call metric_dotM(Prod1, u,iq)
                    call metric_dotM(Prod2, u,Prod1)
                    call metric_dotM(Prod1, t,Prod2)
                    call metric_dotM(Prod2, q,Prod1)
                    call metric_dotM(Prod1,iu,Prod2)
                    call metric_dotM(Prod2,iq,Prod1)
                    call metric_Tm_int(1,iG,Prod2)
                  EndIf
                EndIf
! ------------   SUM = -(X+Y)/2   ------------------------------------
              ElseIf(2*(iG(1,2)+iG(1,3)+iG(2,3)) == -(iG(1,1)+iG(2,2)))  &
     &                                                             then
                If( 2*iG(2,3) == -iG(2,2) ) then
                  If( 2*iG(1,2) == -iG(1,1) ) then
!dbg                   write(6,*) 'Apply  t9 = q.it.iq.t.t.iu.iq.iu.it'
                    call metric_dotM(Prod1,iu,it)
                    call metric_dotM(Prod2,iq,Prod1)
                    call metric_dotM(Prod1,iu,Prod2)
                    call metric_dotM(Prod2, t,Prod1)
                    call metric_dotM(Prod1, t,Prod2)
                    call metric_dotM(Prod2,iq,Prod1)
                    call metric_dotM(Prod1,it,Prod2)
                    call metric_dotM(Prod2, q,Prod1)
                    call metric_Tm_int(1,iG,Prod2)
                  Else
!dbg                   write(6,*) 'Apply  t8 = iq.u.iq.it.s'
                    call metric_dotM(Prod1,it, s)
                    call metric_dotM(Prod2,iq,Prod1)
                    call metric_dotM(Prod1, u,Prod2)
                    call metric_dotM(Prod2,iq,Prod1)
                    call metric_Tm_int(1,iG,Prod2)
                  EndIf
                Else
                  If( 2*iG(1,2) == -iG(1,1)) then
                    If( 2*iG(1,3) == -iG(1,1)) then
!dbg                     write(6,*) 'Apply  t6 = q.u.iq.t.q.it.q.s'
                      call metric_dotM(Prod1, u,iq)
                      call metric_dotM(Prod2, t,Prod1)
                      call metric_dotM(Prod1, q,Prod2)
                      call metric_dotM(Prod2,iu,Prod1)
                      call metric_dotM(Prod1, q,Prod2)
                      call metric_Tm_int(1,iG,Prod1)
                    Else
!dbg                     write(6,*) 'Apply  t12 = s.t.iu.it'
                      call metric_dotM(Prod1,iu,it)
                      call metric_dotM(Prod2, t,Prod1)
                      call metric_dotM(Prod1, s,Prod2)
                      call metric_Tm_int(1,iG,Prod1)
                    EndIf
                  Else
                    If( iG(1,2) == iG(1,3)) then
                      If( 2*iG(1,3)+iG(1,2) > -iG(1,1) ) then
!dbg                       write(6,*) 'Apply  t4 = iq.iu.q.t.u.iq'
                        call metric_dotM(Prod1, u,iq)
                        call metric_dotM(Prod2, u,Prod1)
                        call metric_dotM(Prod1, t,Prod2)
                        call metric_dotM(Prod2, q,Prod1)
                        call metric_dotM(Prod1,iu,Prod2)
                        call metric_dotM(Prod2,iq,Prod1)
                        call metric_Tm_int(1,iG,Prod2)
                      EndIf
                    Else
                      If( 2*iG(1,3)+iG(1,2) > -iG(1,1) ) then
!dbg                       write(6,*) 'Apply  t3 = s.q.u.it.iq.t.t.iu
!c     &                                                    .iq.it'
                        call metric_dotM(Prod1,iq,it)
                        call metric_dotM(Prod2,iu,Prod1)
                        call metric_dotM(Prod1, t,Prod2)
                        call metric_dotM(Prod2, t,Prod1)
                        call metric_dotM(Prod1,iq,Prod2)
                        call metric_dotM(Prod2,it,Prod1)
                        call metric_dotM(Prod1, u,Prod2)
                        call metric_dotM(Prod2, q,Prod1)
                        call metric_dotM(Prod1, s,Prod2)
                        call metric_Tm_int(1,iG,Prod1)
                      EndIf
                    EndIf
                  EndIf
                EndIf
! ------------   SUM > (less negative) -(X+Y)/2   --------------------
              Else
                If( 2*iG(1,2) == -iG(1,1) ) then
                  If( iG(1,3) /= 0 ) then
!dbg                   write(6,*) 'Apply  t12 = s.t.iu.it'
                    call metric_dotM(Prod1,iu,it)
                    call metric_dotM(Prod2, t,Prod1)
                    call metric_dotM(Prod1, s,Prod2)
                    call metric_Tm_int(1,iG,Prod1)
                  EndIf
                Else
                  If( 2*iG(2,3) == -iG(2,2) ) then
                    If( iG(1,3) == 0 ) then
                      If( iG(1,2) /= 0 ) then
!dbg                       write(6,*) 'Apply  t7 = s.iq.t.u.iq.it.s'
                        call metric_dotM(Prod1,it, s)
                        call metric_dotM(Prod2,iq,Prod1)
                        call metric_dotM(Prod1, u,Prod2)
                        call metric_dotM(Prod2, t,Prod1)
                        call metric_dotM(Prod1,iq,Prod2)
                        call metric_dotM(Prod2, s,Prod1)
                        call metric_Tm_int(1,iG,Prod2)
                      EndIf
                    Else
!dbg                     write(6,*) 'Apply  t8 = iq.u.iq.it.s'
                      call metric_dotM(Prod1,it, s)
                      call metric_dotM(Prod2,iq,Prod1)
                      call metric_dotM(Prod1, u,Prod2)
                      call metric_dotM(Prod2,iq,Prod1)
                      call metric_Tm_int(1,iG,Prod2)
                    EndIf
                  EndIf
                EndIf
              EndIf
            EndIf

!     X   X   X   X   XXXXX
!      X X     X X       X
!       X       X       X
!      X X      X      X
!     X   X     X     XXXXX

          Else
! =================== POSITIVE REDUCED =================================
            If( iG(1,2) > 0 .and. iG(1,3) > 0 .and. iG(2,3) > 0 .or.     &
     &          iG(1,2) > 0 .and. iG(1,3) < 0 .and. iG(2,3) < 0 .or.     &
     &          iG(1,2) < 0 .and. iG(1,3) > 0 .and. iG(2,3) < 0 .or.     &
     &          iG(1,2) < 0 .and. iG(1,3) < 0 .and. iG(2,3) > 0 )then
               If( 2*iG(1,2) == iG(1,1) ) then
                If( 2*iG(1,3) /= iG(1,1) ) then
                  If( 2*iG(2,3) /= iG(2,2) ) then
                    If( 2*iG(2,3) < iG(1,3) ) then
!dbg                     write(6,*) 'Apply  t15 = s.t.iu.q.it.q'
                      call metric_dotM(Prod1,it, q)
                      call metric_dotM(Prod2, q,Prod1)
                      call metric_dotM(Prod1,iu,Prod2)
                      call metric_dotM(Prod2, t,Prod1)
                      call metric_dotM(Prod1, s,Prod2)
                      call metric_Tm_int(1,iG,Prod1)
                    EndIf
                  Else
                    If( 2*iG(1,3) < iG(1,2) ) then
!dbg                     write(6,*) 'Apply  t1 = s.q.it.iq.t.t.u.iq.it'
                      call metric_dotM(Prod1,iq,it)
                      call metric_dotM(Prod2, u,Prod1)
                      call metric_dotM(Prod1, t,Prod2)
                      call metric_dotM(Prod2, t,Prod1)
                      call metric_dotM(Prod1,iq,Prod2)
                      call metric_dotM(Prod2,it,Prod1)
                      call metric_dotM(Prod1, q,Prod2)
                      call metric_dotM(Prod2, s,Prod1)
                      call metric_Tm_int(1,iG,Prod2)
                    EndIf
                  EndIf
                Else
                  If( 2*iG(2,3) /= iG(2,2) .and.                         &
     &                2*iG(2,3) < iG(1,3) ) then
!dbg                   write(6,*) 'Apply  t2 = s.q.iu.iq'
                    call metric_dotM(Prod1,iu,iq)
                    call metric_dotM(Prod2, q,Prod1)
                    call metric_dotM(Prod1, s,Prod2)
                    call metric_Tm_int(1,iG,Prod1)
                  EndIf
                EndIf
              Else
                If( 2*iG(1,3) == iG(1,1) ) then
                  If( 2*iG(2,3) /= iG(2,2) .and.                         &
     &                2*iG(2,3) < iG(1,2) ) then
!dbg                   write(6,*) 'Apply  t2 = s.q.iu.iq'
                    call metric_dotM(Prod1,iu,iq)
                    call metric_dotM(Prod2, q,Prod1)
                    call metric_dotM(Prod1, s,Prod2)
                    call metric_Tm_int(1,iG,Prod1)
                  EndIf
                Else
                  If( 2*iG(2,3) == iG(2,2) ) then
                    If( 2*iG(1,3) < iG(1,2) ) then
!dbg                     write(6,*) 'Apply  t1 = s.q.it.iq.t.t.u.iq.it'
                      call metric_dotM(Prod1,iq,it)
                      call metric_dotM(Prod2, u,Prod1)
                      call metric_dotM(Prod1, t,Prod2)
                      call metric_dotM(Prod2, t,Prod1)
                      call metric_dotM(Prod1,iq,Prod2)
                      call metric_dotM(Prod2,it,Prod1)
                      call metric_dotM(Prod1, q,Prod2)
                      call metric_dotM(Prod2, s,Prod1)
                      call metric_Tm_int(1,iG,Prod2)
                    EndIf
                  EndIf
                EndIf
              EndIf
! =================== NEGATIVE REDUCED =================================
            Else
! ------------   SUM < (more negative) -(X+Y)/2   ---------------------
              If( 2*(iG(1,2)+iG(1,3)+iG(2,3)) < -(iG(1,1)+iG(2,2)) )     &
     &                                                             then
                If( 2*iG(2,3) == -iG(2,2) ) then
                  If( 2*iG(1,2) == -iG(1,1) ) then
                    If( 2*iG(1,3) == -iG(1,1) ) then
!dbg                     write(6,*) 'Apply  t4 = iq.iu.q.t.u.iq'
                      call metric_dotM(Prod1, u,iq)
                      call metric_dotM(Prod2, u,Prod1)
                      call metric_dotM(Prod1, t,Prod2)
                      call metric_dotM(Prod2, q,Prod1)
                      call metric_dotM(Prod1,iu,Prod2)
                      call metric_dotM(Prod2,iq,Prod1)
                      call metric_Tm_int(1,iG,Prod2)
                    Else
!dbg                     write(6,*) 'Apply t9 = q.it.iq.t.t.iu.iq.iu.it'
                      call metric_dotM(Prod1,iu,it)
                      call metric_dotM(Prod2,iq,Prod1)
                      call metric_dotM(Prod1,iu,Prod2)
                      call metric_dotM(Prod2, t,Prod1)
                      call metric_dotM(Prod1, t,Prod2)
                      call metric_dotM(Prod2,iq,Prod1)
                      call metric_dotM(Prod1,it,Prod2)
                      call metric_dotM(Prod2, q,Prod1)
                      call metric_Tm_int(1,iG,Prod2)
                    EndIf
                  Else
!dbg                   write(6,*) 'Apply  t3 = s.q.u.it.iq.t.t.iu
!dbg     &                                                      .iq.it'
                    call metric_dotM(Prod1,iq,it)
                    call metric_dotM(Prod2,iu,Prod1)
                    call metric_dotM(Prod1, t,Prod2)
                    call metric_dotM(Prod2, t,Prod1)
                    call metric_dotM(Prod1,iq,Prod2)
                    call metric_dotM(Prod2,it,Prod1)
                    call metric_dotM(Prod1, u,Prod2)
                    call metric_dotM(Prod2, q,Prod1)
                    call metric_dotM(Prod1, s,Prod2)
                    call metric_Tm_int(1,iG,Prod1)
                  EndIf
                Else
                  If( 2*iG(1,2) == -iG(1,1) ) then
                    If( 2*iG(1,3) == -iG(1,1) ) then
!dbg                     write(6,*) 'Apply t3 = s.q.u.it.iq.t.t.iu.iq.it'
                      call metric_dotM(Prod1,iq,it)
                      call metric_dotM(Prod2,iu,Prod1)
                      call metric_dotM(Prod1, t,Prod2)
                      call metric_dotM(Prod2, t,Prod1)
                      call metric_dotM(Prod1,iq,Prod2)
                      call metric_dotM(Prod2,it,Prod1)
                      call metric_dotM(Prod1, u,Prod2)
                      call metric_dotM(Prod2, q,Prod1)
                      call metric_dotM(Prod1, s,Prod2)
                      call metric_Tm_int(1,iG,Prod1)
                    Else
!dbg                     write(6,*) 'Apply  t9 = q.it.iq.t.t.iu.iq.iu.it'
                      call metric_dotM(Prod1,iu,it)
                      call metric_dotM(Prod2,iq,Prod1)
                      call metric_dotM(Prod1,iu,Prod2)
                      call metric_dotM(Prod2, t,Prod1)
                      call metric_dotM(Prod1, t,Prod2)
                      call metric_dotM(Prod2,iq,Prod1)
                      call metric_dotM(Prod1,it,Prod2)
                      call metric_dotM(Prod2, q,Prod1)
                      call metric_Tm_int(1,iG,Prod2)
                    EndIf
                  Else
!dbg                   write(6,*) 'Apply t3 = s.q.u.it.iq.t.t.iu.iq.it'
                    call metric_dotM(Prod1,iq,it)
                    call metric_dotM(Prod2,iu,Prod1)
                    call metric_dotM(Prod1, t,Prod2)
                    call metric_dotM(Prod2, t,Prod1)
                    call metric_dotM(Prod1,iq,Prod2)
                    call metric_dotM(Prod2,it,Prod1)
                    call metric_dotM(Prod1, u,Prod2)
                    call metric_dotM(Prod2, q,Prod1)
                    call metric_dotM(Prod1, s,Prod2)
                    call metric_Tm_int(1,iG,Prod1)
                  EndIf
                EndIf
! ------------   SUM = -(X+Y)/2   --------------------------------------
              ElseIf(2*(iG(1,2)+iG(1,3)+iG(2,3)) == -(iG(1,1)+iG(2,2)))  &
     &                                                             then
                If( 2*iG(2,3) == -iG(2,2) ) then
                  If( 2*iG(1,2) == -iG(1,1) ) then
!dbg                   write(6,*) 'Apply  t9 = q.it.iq.t.t.iu.iq.iu.it'
                    call metric_dotM(Prod1,iu,it)
                    call metric_dotM(Prod2,iq,Prod1)
                    call metric_dotM(Prod1,iu,Prod2)
                    call metric_dotM(Prod2, t,Prod1)
                    call metric_dotM(Prod1, t,Prod2)
                    call metric_dotM(Prod2,iq,Prod1)
                    call metric_dotM(Prod1,it,Prod2)
                    call metric_dotM(Prod2, q,Prod1)
                    call metric_Tm_int(1,iG,Prod2)
                  ElseIf ( 2*iG(1,3) /= -iG(1,1) ) then
!dbg                   write(6,*) 'Apply  t7 = s.iq.t.u.iq.it.s'
                    call metric_dotM(Prod1,it, s)
                    call metric_dotM(Prod2,iq,Prod1)
                    call metric_dotM(Prod1, u,Prod2)
                    call metric_dotM(Prod2, t,Prod1)
                    call metric_dotM(Prod1,iq,Prod2)
                    call metric_dotM(Prod2, s,Prod1)
                    call metric_Tm_int(1,iG,Prod2)
                  EndIf
                Else
                  If( 2*iG(1,2) == -iG(1,1)) then
!dbg                   write(6,*) 'Apply  t12 = s.t.iu.it'
                    call metric_dotM(Prod1,iu,it)
                    call metric_dotM(Prod2, t,Prod1)
                    call metric_dotM(Prod1, s,Prod2)
                    call metric_Tm_int(1,iG,Prod1)
                  Else
                    If( 2*iG(1,3) == -iG(1,1)) then
!dbg                     write(6,*) 'Apply  t6 = q.u.iq.t.q.it.q.s'
                      call metric_dotM(Prod1, u,iq)
                      call metric_dotM(Prod2, t,Prod1)
                      call metric_dotM(Prod1, q,Prod2)
                      call metric_dotM(Prod2,iu,Prod1)
                      call metric_dotM(Prod1, q,Prod2)
                      call metric_Tm_int(1,iG,Prod1)
                    Else
                      If( 2*iG(1,3)+iG(1,2) > -iG(1,1) ) then
!dbg                       write(6,*) 'Apply  t3 = s.q.u.it.iq.t.t.iu
!dbg     &                                                   .iq.it'
                        call metric_dotM(Prod1,iq,it)
                        call metric_dotM(Prod2,iu,Prod1)
                        call metric_dotM(Prod1, t,Prod2)
                        call metric_dotM(Prod2, t,Prod1)
                        call metric_dotM(Prod1,iq,Prod2)
                        call metric_dotM(Prod2,it,Prod1)
                        call metric_dotM(Prod1, u,Prod2)
                        call metric_dotM(Prod2, q,Prod1)
                        call metric_dotM(Prod1, s,Prod2)
                        call metric_Tm_int(1,iG,Prod1)
                      EndIf
                    EndIf
                  EndIf
                EndIf
! ------------   SUM > (less negative) -(X+Y)/2   ------------------
              Else
                If( 2*iG(1,2) == -iG(1,1) ) then
                  If( iG(1,3) /= 0 ) then
!dbg                   write(6,*) 'Apply  t12 = s.t.iu.it'
                    call metric_dotM(Prod1,iu,it)
                    call metric_dotM(Prod2, t,Prod1)
                    call metric_dotM(Prod1, s,Prod2)
                    call metric_Tm_int(1,iG,Prod1)
                  EndIf
                Else
                  If( 2*iG(2,3) == -iG(2,2) ) then
                    If( iG(1,2) /= 0 ) then
!dbg                     write(6,*) 'Apply  t7 = s.iq.t.u.iq.it.s'
                      call metric_dotM(Prod1,it, s)
                      call metric_dotM(Prod2,iq,Prod1)
                      call metric_dotM(Prod1, u,Prod2)
                      call metric_dotM(Prod2, t,Prod1)
                      call metric_dotM(Prod1,iq,Prod2)
                      call metric_dotM(Prod2, s,Prod1)
                      call metric_Tm_int(1,iG,Prod2)
                    EndIf
                  Else
                    If( 2*iG(1,3) == -iG(1,1) ) then
                      If( iG(1,2) /= 0 ) then
!dbg                       write(6,*) 'Apply  t6 = q.u.iq.t.q.it.q.s'
                        call metric_dotM(Prod1, u,iq)
                        call metric_dotM(Prod2, t,Prod1)
                        call metric_dotM(Prod1, q,Prod2)
                        call metric_dotM(Prod2,iu,Prod1)
                        call metric_dotM(Prod1, q,Prod2)
                        call metric_Tm_int(1,iG,Prod1)
                      EndIf
                    EndIf
                  EndIf
                EndIf
              EndIf
            EndIf
          EndIf
        EndIf
      EndIf



!dbg      write(6,*) '  arrived at Niggli cell'
!dbg      write(6,701) (iG(1,i),i=1,3)
!dbg      write(6,701) (iG(2,i),i=1,3)
!dbg      write(6,701) (iG(3,i),i=1,3)



      call metric_Tm_int(2,iG,Mtotal)


      return
      end

