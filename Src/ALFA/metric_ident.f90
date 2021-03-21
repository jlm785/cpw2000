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
!>    it identifies the lattice, generates the Niggli lattice
!>    and gives the canonical primitive lattice vectors that
!>    generate g_ij

      subroutine metric_ident(Gin,G,Gout,bravais,Mtotal,avec,aconv,      &
     &                        avecnig,tol,iprint)



!     output:
!     G          metric of the Niggli lattice
!     Gout       same as Gin but without "noise"
!     bravais    bravais lattice (1=sc,2=fcc,3=bcc,)
!     Mtotal     integer matrix  G=Mtotal Gin (Mtotal)T, det Mtotal=1
!     avec(n,j)  j-th canonical primitive lattice vectors that generate Gin
!                (n-th cartesian coordinate)

!     written by Alvaro Ladeira, 1999
!     modified Jose Luis Martins, Fevereiro 2000
!     corrected bugs, jlm, December 2000.
!     modified by Jose Luis Martins to use integer logic
!     in the niggli recognition. April 2004
!     modified documentation August 2019. JLM
!     copyright Alvaro Ladeira/Jose Luis Martins/INESC-MN

!     send comments/bug reports to jlmartins@inesc-mn.pt

!     version 4.94


      implicit none

      integer, parameter          :: REAL64 = selected_real_kind(12)
      integer, parameter          :: INT64 = selected_int_kind(15)

!     input

      real(REAL64), intent(in)           ::  Gin(3,3)                    !<  input metric
      real(REAL64), intent(in)           ::  tol                         !<  tolerance for considering a difference is negligible
      integer, intent(in)                ::  iprint                      !<  prints information if iprint /= 0

!     output


      real(REAL64), intent(out)          ::  G(3,3)                      !<  metric of the Niggli lattice
      real(REAL64), intent(out)          ::  Gout(3,3)                   !<  same as Gin but without "noise"
      integer, intent(out)               ::  bravais                     !<  bravais lattice (1=sc,2=fcc,3=bcc,4=st,5=ct,6=so,7=bco,8=fco,9=bco,10=hex,11=romb,12=sm,13=cm,14=tri)
      integer, intent(out)               ::  Mtotal(3,3)                 !<  G=Mtotal Gin (Mtotal)T, det Mtotal=1
      real(REAL64), intent(out)          ::  avec(3,3)                   !<  primitive lattice vectors that generate Gin in canonical orientation
      real(REAL64), intent(out)          ::  aconv(3,3)                  !<  conventional lattice vectors
      real(REAL64), intent(out)          ::  avecnig(3,3)                !<  primitive lattice vectors of the Niggli cell in canonical orientation

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
      real(REAL64)    ::  error
      real(REAL64)    ::  scal, detG, cubrdetG, scltol, skew
      integer         ::  rede
      integer         ::  icount

      integer         ::  invM(3,3)
      real(REAL64)    ::  Gavec(3,3)

      integer(INT64)  ::  iGin(3,3)
      integer(INT64)  ::  iG(3,3)
      
      integer         ::  Mtot1(3,3),Mtot2(3,3)
      real(REAL64)    ::  scalint
      
      logical         ::  lmetric                                        !  True if g is a metric
      

!     constants

      real(REAL64), parameter :: ZERO = 0.0_REAL64
      real(REAL64), parameter :: EPS = 1.0d-12

!     counters

      integer         ::  i, j, k



      do i=1,3
      do j=1,3
        G(i,j) = Gin(i,j)
      enddo
      enddo

!     resets Mtotal

      call metric_Tm(0,G,Mtotal)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C        TEST: Is G a metric?                                         C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      call  metric_test_g(G,lmetric)


      if(.NOT. lmetric) stop


      scal = ZERO
      do i=1,3
      do j=1,3
        scal = max(scal,abs(G(i,j)))
      enddo
      enddo


      detG = G(1,1)*G(2,2)*G(3,3) - G(1,1)*G(2,3)*G(3,2)                 &
     &     + G(1,2)*G(2,3)*G(3,1) - G(1,2)*G(2,1)*G(3,3)                 &
     &     + G(1,3)*G(2,1)*G(3,2) - G(1,3)*G(2,2)*G(3,1)

      cubrdetG = detG**(1.0/3.0)

      skew = scal/cubrdetG

!     checks range of tol

      scltol = tol
      if(tol > 0.01) then
        scltol = 0.01
        if(iprint >= 1) write(6,*) 'metric_ident   tol=0.01'
      endif

      if(tol < 5.0*eps) then
        scltol = 5.0*eps
        if(iprint >= 1) write(6,*) 'metric_ident   tol= ',5.0*eps
      endif



      if( cubrdetG < 10.0*scltol*scal) then

         write(6,*) ' metric_ident'
         write(6,*) ' det(G),scale,tol ',detG,scal,scltol
         write(6,*) 'The cell is too skewed for the choice of tol'

         stop

      endif

      scltol = cubrdetG*scltol

!     CHANGES INTO INTEGER ARITHMETIC

      scalint = ZERO
      do i=1,3
      do j=1,3
        scalint = max(scalint,abs(Gin(i,j)))
      enddo
      enddo


      scalint = scalint / 2.0d0**50

      do i=1,3
      do j=1,3
        iGin(i,j) = int(G(i,j)/scalint,INT64)
      enddo
      enddo

      iGin(1,2) = (iGin(1,2)+iGin(2,1))/2
      iGin(2,1) = iGin(1,2)
      iGin(2,3) = (iGin(2,3)+iGin(3,2))/2
      iGin(3,2) = iGin(2,3)
      iGin(3,1) = (iGin(3,1)+iGin(1,3))/2
      iGin(1,3) = iGin(3,1)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C                                                                     C
!C  STEP 1+2 OBTAINS THE NIGGLI CELL                                   C
!C                                                                     C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      call metric_niggli(iGin,iG,Mtot1)

      do i=1,3
      do j=1,3
        G(i,j) = dble(iG(i,j))*scalint
      enddo
      enddo


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C                                                                     C
!C  STEP 3: TRAPS SPECIAL CASES OF NOISE EFFECTS                       C
!C                                                                     C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      rede=0
      bravais=0

!     FCC

      If( abs(G(1,1)-G(2,2)) < scltol .and.                              &
     &    abs(G(2,2)-G(3,3)) < scltol .and.                              &
     &    abs(G(1,2) - G(1,1) / 2) < scltol .and.                        &
     &      abs(G(1,3)) < scltol .and.                                   &
     &      abs(G(2,3) - G(2,2) / 2) < scltol ) then
          G(1,1) = (G(1,1)+G(2,2)+G(3,3)) / 3
          G(2,2) = G(1,1)
          G(3,3) = G(1,1)
          G(1,2) = G(1,1) / 2
          G(2,1) = G(1,2)
          G(1,3) = ZERO
          G(3,1) = ZERO
          G(2,3) = G(1,1) / 2
          G(3,2) = G(2,3)
        call metric_dotM(Prod1,iq,it)
        call metric_dotM(Prod2, t,Prod1)
        call metric_dotM(Prod1,iq,Prod2)
        call metric_dotM(Prod2, s,Prod1)
        call metric_Tm(1,G,Prod2)
!dbg             write(6,*) 'face centred cubic  special-4'
      EndIf

      If( abs(G(1,1)-G(2,2)) < scltol .and.                              &
     &    abs(G(2,2)-G(3,3)) < scltol .and.                              &
     &    abs(G(1,2) - G(1,1) / 2) < scltol .and.                        &
     &      abs(G(1,3) - G(1,1) / 2) < scltol .and.                      &
     &      abs(G(2,3)) < scltol ) then
          G(1,1) = (G(1,1)+G(2,2)+G(3,3)) / 3
          G(2,2) = G(1,1)
          G(3,3) = G(1,1)
          G(1,2) = G(1,1) / 2
          G(2,1) = G(1,2)
          G(1,3) = G(1,1) / 2
          G(3,1) = G(1,3)
          G(2,3) = ZERO
          G(3,2) = ZERO
        call metric_dotM(Prod1, t,t)
        call metric_dotM(Prod2, q,Prod1)
        call metric_dotM(Prod1, t,Prod2)
        call metric_dotM(Prod2, t,Prod1)
        call metric_dotM(Prod1,iq,Prod2)
        call metric_dotM(Prod2, s,Prod1)
        call metric_Tm(1,G,Prod2)
!dbg             write(6,*) 'face centred cubic  special-5'
      EndIf

      If( abs(G(1,1)-G(2,2)) < scltol .and.                              &
     &    abs(G(2,2)-G(3,3)) < scltol .and.                              &
     &    abs(G(1,2)) < scltol .and.                                     &
     &      abs(G(1,3) - G(1,1) / 2) < scltol .and.                      &
     &      abs(G(2,3) - G(2,2) / 2) < scltol ) then
          G(1,1) = (G(1,1)+G(2,2)+G(3,3)) / 3
          G(2,2) = G(1,1)
          G(3,3) = G(1,1)
          G(1,2) = ZERO
          G(2,1) = ZERO
          G(1,3) = G(1,1) / 2
          G(3,1) = G(1,3)
          G(2,3) = G(1,1) / 2
          G(3,2) = G(2,3)
             call metric_Tm(1,G, s)
!dbg             write(6,*) 'face centred cubic  special-6'
      EndIf

      If( abs(G(1,1)-G(2,2)) < scltol .and.                              &
     &    abs(G(2,2)-G(3,3)) < scltol .and.                              &
     &    abs(G(1,2)) < scltol .and.                                     &
     &      abs(G(1,3) + G(1,1) / 2) < scltol .and.                      &
     &      abs(G(2,3) + G(2,2) / 2) < scltol ) then
          G(1,1) = (G(1,1)+G(2,2)+G(3,3)) / 3
          G(2,2) = G(1,1)
          G(3,3) = G(1,1)
          G(1,2) = ZERO
          G(2,1) = ZERO
          G(1,3) = -G(1,1) / 2
          G(3,1) = G(1,3)
          G(2,3) = -G(1,1) / 2
          G(3,2) = G(2,3)
             call metric_Tm(1,G, iq)
!dbg             write(6,*) 'face centred cubic  special-2'
      EndIf

      If( abs(G(1,1)-G(2,2)) < scltol .and.                              &
     &    abs(G(2,2)-G(3,3)) < scltol .and.                              &
     &    abs(G(1,2) + G(1,1) / 2) < scltol .and.                        &
     &      abs(G(1,3)) < scltol .and.                                   &
     &      abs(G(2,3) + G(2,2) / 2) < scltol ) then
          G(1,1) = (G(1,1)+G(2,2)+G(3,3)) / 3
          G(2,2) = G(1,1)
          G(3,3) = G(1,1)
          G(1,2) = -G(1,1) / 2
          G(2,1) = G(1,2)
          G(1,3) = ZERO
          G(3,1) = ZERO
          G(2,3) = -G(1,1) / 2
          G(3,2) = G(2,3)
            call metric_dotM(Prod1, t,iq)
            call metric_dotM(Prod2,iq,Prod1)
            call metric_dotM(Prod1, s,Prod2)
            call metric_Tm(1,G,Prod1)
!dbg             write(6,*) 'face centred cubic  special-3'
      EndIf

      If( abs(G(1,1)-G(2,2)) < scltol .and.                              &
     &    abs(G(2,2)-G(3,3)) < scltol .and.                              &
     &    abs(G(1,2) + G(1,1) / 2) < scltol .and.                        &
     &      abs(G(1,3) + G(1,1) / 2) < scltol .and.                      &
     &      abs(G(2,3)) < scltol ) then
          G(1,1) = (G(1,1)+G(2,2)+G(3,3)) / 3
          G(2,2) = G(1,1)
          G(3,3) = G(1,1)
          G(1,2) = -G(1,1) / 2
          G(2,1) = G(1,2)
          G(1,3) = -G(1,1) / 2
          G(3,1) = G(1,3)
          G(2,3) = ZERO
          G(3,2) = ZERO
                call metric_dotM(Prod1,iu,it)
                call metric_dotM(Prod2, t,Prod1)
                call metric_dotM(Prod1, s,Prod2)
                call metric_Tm(1,G,Prod1)
!dbg             write(6,*) 'face centred cubic  special-1'
      EndIf


!     HEXAGONAL

        If( abs(G(1,1) - G(2,2)) < scltol .and.                          &
     &      abs(G(1,2) - G(1,1) / 2) < scltol .and.                      &
     &      abs(G(2,3)) < scltol .and.                                   &
     &      abs(G(1,3)) < scltol ) then
          G(1,1) = (G(1,1) + G(2,2)) / 2
          G(2,2) = G(1,1)
          G(1,2) = G(1,1) / 2
          G(2,1) = G(1,2)
          G(2,3) = ZERO
          G(3,2) = ZERO
          G(1,3) = ZERO
          G(3,1) = ZERO
!dbg           write(6,*) 'simple hexagonal special-1'
          call metric_dotM(Prod1, t,t)
          call metric_dotM(Prod2, q,Prod1)
          call metric_dotM(Prod1, t,Prod2)
          call metric_dotM(Prod2, t,Prod1)
          call metric_dotM(Prod1,iq,Prod2)
          call metric_dotM(Prod2, s,Prod1)
          call metric_Tm(1,G,Prod2)
        EndIf

       If( abs(G(3,3) - G(2,2)) < scltol .and.                           &
     &      abs(G(2,3) - G(2,2) / 2) < scltol .and.                      &
     &      abs(G(1,3)) < scltol .and.                                   &
     &      abs(G(1,2)) < scltol ) then
          G(3,3) = (G(3,3) + G(2,2)) / 2
          G(2,2) = G(3,3)
          G(2,3) = G(2,2) / 2
          G(3,2) = G(2,3)
          G(1,3) = ZERO
          G(3,1) = ZERO
          G(1,2) = ZERO
          G(2,1) = ZERO
!dbg           write(6,*) 'simple hexagonal special-3'
          call metric_Tm(1,G,s)
       EndIf

       If( abs(G(3,3) - G(1,1)) < scltol .and.                           &
     &      abs(G(1,3) - G(1,1) / 2) < scltol .and.                      &
     &      abs(G(2,3)) < scltol .and.                                   &
     &      abs(G(1,2)) < scltol ) then
          G(3,3) = (G(3,3) + G(1,1)) / 2
          G(1,1) = G(3,3)
          G(1,2) = ZERO
          G(2,1) = ZERO
          G(2,3) = ZERO
          G(3,2) = ZERO
          G(1,3) = G(3,3) / 2
          G(3,1) = G(1,3)
!dbg           write(6,*) 'simple hexagonal special-2'
          call metric_Tm(1,G,s)
        Endif

!     CENTERED TETRAGONAL

        If( bravais == 0 .and.                                           &
     &     abs(G(1,1)-G(2,2)) < scltol .and.                             &
     &     abs(G(2,2)-G(3,3)) < scltol .and.                             &
     &     abs(G(1,2) - G(2,3)) < scltol .and.                           &
     &     abs(G(1,1) + G(2,3) + G(1,2) + G(1,3)) < scltol ) then
          G(1,1) = (G(1,1)+G(2,2)+G(3,3)) / 3
          G(2,2) = G(1,1)
          G(3,3) = G(1,1)
          G(1,2) = (G(1,2)+G(2,3)) / 2
          G(2,1) = G(1,2)
          G(2,3) = G(1,2)
          G(3,2) = G(1,2)
          G(1,3) = -G(1,1) - G(1,2) - G(2,3)
          G(3,1) = G(1,3)
!dbg            write(6,*) 'centred tetragonal special-3'
            call metric_Tm(1,G,q)
        Endif

        If( bravais == 0 .and.                                           &
     &     abs(G(1,1) - G(2,2)) < scltol .and.                           &
     &     abs(G(1,2)) < scltol .and.                                    &
     &     abs(G(1,3) - G(1,1) / 2) < scltol .and.                       &
     &     abs(G(2,3) - G(2,2) / 2) < scltol) then 
          G(1,1) = (G(1,1) + G(2,2)) / 2
          G(2,2) = G(1,1)
          G(1,2) = ZERO
          G(2,1) = ZERO
          G(1,3) = G(1,1) / 2
          G(3,1) = G(1,3)
          G(2,3) = G(2,2) / 2
          G(3,2) = G(2,3)
!dbg           write(6,*) 'centred tetragonal special-1'
          call metric_Tm(1,G,s)
        Endif

!     RHOMBOHEDRAL

      If( abs(G(1,1) - G(2,2)) < scltol .and.                            &
     &    abs(G(1,2) - G(1,1) / 2) < scltol .and.                        &
     &    abs(G(1,3)) < scltol .and.                                     &
     &    abs(G(2,3) - G(2,2) / 2) < scltol ) then
          G(1,1) = (G(1,1) + G(2,2)) / 2
          G(2,2) = G(1,1)
          G(1,2) = G(1,1) / 2
          G(2,1) = G(1,2)
          G(1,3) = ZERO
          G(3,1) = ZERO
          G(2,3) = G(2,2) / 2
          G(3,2) = G(2,3)
!dbg           write(6,*) 'rombohedrica special-1'
            call metric_Tm(1,G,iu)
            call metric_Tm(1,G,q)
            call metric_Tm(1,G,s)
            call metric_Tm(1,G,iq)
      EndIf

      If( abs(G(1,1) - G(2,2)) < scltol .and.                            &
     &    abs(G(1,2) + G(1,1) / 2) < scltol .and.                        &
     &    abs(G(1,3)) < scltol .and.                                     &
     &    abs(G(2,3) + G(2,2) / 2) < scltol ) then
          G(1,1) = (G(1,1) + G(2,2)) / 2
          G(2,2) = G(1,1)
          G(1,2) = -G(1,1) / 2
          G(2,1) = G(1,2)
          G(1,3) = ZERO
          G(3,1) = ZERO
          G(2,3) = -G(2,2) / 2
          G(3,2) = G(2,3)
!dbg           write(6,*) 'rombohedrica special-2'
            call metric_Tm(1,G,u)
            call metric_Tm(1,G,s)
      EndIf

      If( abs(G(1,1) - G(2,2)) < scltol .and.                            &
     &    abs(G(1,2) + G(1,1) / 2) < scltol .and.                        &
     &    abs(G(2,3)) < scltol .and.                                     &
     &    abs(G(1,3) + G(2,2) / 2) < scltol ) then
          G(1,1) = (G(1,1) + G(2,2)) / 2
          G(2,2) = G(1,1)
          G(1,2) = -G(1,1) / 2
          G(2,1) = G(1,2)
          G(2,3) = ZERO
          G(3,2) = ZERO
          G(1,3) = -G(1,1) / 2
          G(3,1) = G(1,3)
!dbg           write(6,*) 'rombohedrica special-3'
            call metric_Tm(1,G,t)
            call metric_Tm(1,G,iu)
            call metric_Tm(1,G,q)
            call metric_Tm(1,G,s)
            call metric_Tm(1,G,iq)
      EndIf

      If( abs(G(1,1) - G(2,2)) < scltol .and.                            &
     &    abs(G(1,2) - G(1,1) / 2) < scltol .and.                        &
     &    abs(G(2,3)) < scltol .and.                                     &
     &    abs(G(1,3) - G(2,2) / 2) < scltol ) then
          G(1,1) = (G(1,1) + G(2,2)) / 2
          G(2,2) = G(1,1)
          G(1,2) = G(1,1) / 2
          G(2,1) = G(1,2)
          G(2,3) = ZERO
          G(3,2) = ZERO
          G(1,3) = G(2,2) / 2
          G(3,1) = G(1,3)
!dbg           write(6,*) 'rombohedrica special-4'
            call metric_Tm(1,G,t)
            call metric_Tm(1,G,u)
            call metric_Tm(1,G,s)
      EndIf

!     BASE CENTERED ORTHORHOMBIC

       If( abs(G(1,3) - G(1,1) / 2) < scltol .and.                       &
     &     abs(G(2,3)) < scltol .and.                                    &
     &     abs(G(1,2)) < scltol) then 
          G(1,3) = G(1,1) / 2
          G(3,1) = G(1,3)
          G(2,3) = ZERO
          G(3,2) = ZERO
          G(1,2) = ZERO
          G(2,1) = ZERO
!dbg           write(6,*) 'base centred orthorhombic special-1'
            call metric_Tm(1,G,s)
      EndIf

       If( abs(G(1,2) - G(1,1) / 2) < scltol .and.                       &
     &     abs(G(1,3)) < scltol .and.                                    &
     &     abs(G(2,3)) < scltol) then 
          G(1,2) = G(1,1) / 2
          G(2,1) = G(1,2)
          G(1,3) = ZERO
          G(3,1) = ZERO
          G(2,3) = ZERO
          G(3,2) = ZERO
!dbg           write(6,*) 'base centred orthorhombic special-2'
            call metric_Tm(1,G,q)
            call metric_Tm(1,G,s)
            call metric_Tm(1,G,iq)
       EndIf

       If( abs(G(2,3) - G(2,2) / 2) < scltol .and.                       &
     &     abs(G(1,3)) < scltol .and.                                    &
     &     abs(G(1,2)) < scltol) then 
          G(2,3) = G(2,2) / 2
          G(3,2) = G(2,3)
          G(1,3) = ZERO
          G(3,1) = ZERO
          G(1,2) = ZERO
          G(2,1) = ZERO
!dbg           write(6,*) 'base centred orthorhombic special-3'
            call metric_Tm(1,G,s)
       EndIf

!      BODY CENTRED ORTHORHOMBIC

       If( abs(G(3,3) - G(2,2)) < scltol .and.                           &
     &     abs(G(1,2) + G(1,1) / 2) < scltol .and.                       &
     &     abs(G(1,3) + G(1,1) / 2) < scltol ) then
          G(3,3) = (G(3,3) + G(2,2)) / 2
          G(2,2) = G(3,3)
          G(1,2) = -G(1,1) / 2
          G(2,1) = G(1,2)
          G(1,3) = -G(1,1) / 2
          G(3,1) = G(1,3)
!dbg           write(6,*) 'body centred orthorombic special-1'
            call metric_Tm(1,G,q)
            call metric_Tm(1,G,s)
            call metric_Tm(1,G,iq)
       EndIf

       If( abs(G(1,3) - G(1,1) / 2) < scltol .and.                       &
     &     abs(G(2,3) - G(2,2) / 2) < scltol .and.                       &
     &     abs(G(1,2)) < scltol ) then
          G(1,3) = G(1,1) / 2
          G(3,1) = G(1,3)
          G(2,3) = G(2,2) / 2
          G(3,2) = G(2,3)
          G(1,2) = ZERO
          G(2,1) = ZERO
!dbg           write(6,*) 'body centred orthorhombic-6'
            call metric_Tm(1,G,s)
       EndIf

!dbg      write(6,*) '  arrived at modified Niggli cell'
!dbg      write(6,777) (G(1,i),i=1,3)
!dbg      write(6,777) (G(2,i),i=1,3)
!dbg      write(6,777) (G(3,i),i=1,3)



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C                                                                     C
!C  STEP 4: IDENTIFIES LATTICE FROM NIGGLI                             C
!C                                                                     C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!     SEARCHES FOR CUBIC LATTICES

         If( abs(G(1,1)-G(2,2)) < scltol .and.                           &
     &      abs(G(2,2)-G(3,3)) < scltol .and.                            &
     &      abs(G(1,2) - G(1,1) / 2) < scltol .and.                      &
     &      abs(G(1,3) - G(1,1) / 2) < scltol .and.                      &
     &      abs(G(2,3) - G(2,2) / 2) < scltol ) then
          G(1,1) = (G(1,1)+G(2,2)+G(3,3)) / 3
          G(2,2) = G(1,1)
          G(3,3) = G(1,1)
          G(1,2) = G(1,1) / 2
          G(2,1) = G(1,2)
          G(1,3) = G(1,1) / 2
          G(3,1) = G(1,3)
          G(2,3) = G(1,1) / 2
          G(3,2) = G(2,3)
!dbg              write(6,*) 'face centred cubic'
          rede=1
          bravais=2
        EndIf

!       sc

      If( bravais == 0 .and.                                             &
     &      abs(G(1,1)-G(2,2)) < scltol .and.                            &
     &      abs(G(2,2)-G(3,3)) < scltol .and.                            &
     &      abs(G(1,2)) < scltol .and.                                   &
     &      abs(G(1,3)) < scltol .and.                                   &
     &      abs(G(2,3)) < scltol ) then
          G(1,1) = (G(1,1)+G(2,2)+G(3,3)) / 3
          G(2,2) = G(1,1)
          G(3,3) = G(1,1)
          G(1,2) = ZERO
          G(2,1) = ZERO
          G(1,3) = ZERO
          G(3,1) = ZERO
          G(2,3) = ZERO
          G(3,2) = ZERO
!dbg           write(6,*) 'simple cubic'
           rede=16
           bravais=1
         EndIf

!       bcc

      If( bravais == 0 .and.                                             &
     &      abs(G(1,1)-G(2,2)) < scltol .and.                            &
     &      abs(G(2,2)-G(3,3)) < scltol .and.                            &
     &      abs(G(1,2) + G(1,1)/3 ) < scltol .and.                       &
     &      abs(G(1,3) + G(1,1)/3 ) < scltol .and.                       &
     &      abs(G(2,3) + G(2,2)/3 ) < scltol ) then
          G(1,1) = (G(1,1)+G(2,2)+G(3,3)) / 3
          G(2,2) = G(1,1)
          G(3,3) = G(1,1)
          G(1,2) = -G(1,1) / 3
          G(2,1) = G(1,2)
          G(1,3) = -G(1,1) / 3
          G(3,1) = G(1,3)
          G(2,3) = -G(1,1) / 3
          G(3,2) = G(2,3)
!dbg           write(6,*) 'body centred cubic'
          rede=11
          bravais=3
      EndIf

!     SEARCHES FOR RHOMBOHEDRAL LATTICES

      If( bravais == 0 .and.                                             &
     &    (abs(G(1,1)-G(2,2)) < scltol) .and.                            &
     &    (abs(G(2,2)-G(3,3)) < scltol) .and.                            &
     &    (abs(G(1,2)-G(1,3)) < scltol) .and.                            &
     &    (abs(G(1,3)-G(2,3)) < scltol) ) then
          G(1,1) = (G(1,1)+G(2,2)+G(3,3)) / 3
          G(2,2) = G(1,1)
          G(3,3) = G(1,1)
          G(2,3) = (G(2,3) + G(1,3) + G(1,2)) / 3
          G(1,3) = G(2,3)
          G(1,2) = G(2,3)
          G(3,2) = G(2,3)
          G(3,1) = G(2,3)
          G(2,1) = G(2,3)
          If( G(1,2) > 0) then
            rede=7
!dbg             write(6,*) 'rombohedric-1'
          Else
            rede=17
!dbg             write(6,*) 'rombohedric-2'
          EndIf
          bravais=11
        EndIf

      If( bravais == 0 .and.                                             &
     &    abs(G(1,1) - G(2,2)) < scltol .and.                            &
     &    abs(G(1,2) - G(1,1) / 2) < scltol .and.                        &
     &    abs(G(1,3) - G(1,1) / 2) < scltol .and.                        &
     &    abs(G(2,3) - G(2,2) / 2) < scltol ) then
          G(1,1) = (G(1,1) + G(2,2)) / 2
          G(2,2) = G(1,1)
          G(1,2) = G(1,1) / 2
          G(2,1) = G(1,2)
          G(1,3) = G(1,1) / 2
          G(3,1) = G(1,3)
          G(2,3) = G(2,2) / 2
          G(3,2) = G(2,3)
!dbg           write(6,*) 'rombohedrica-3'
          rede=24
          bravais=11
        EndIf

        If( bravais == 0 .and.                                           &
     &      abs(G(2,2) - G(3,3)) < scltol .and.                          &
     &      abs(G(1,2) + G(1,1)/3) < scltol .and.                        &
     &      abs(G(1,2) - G(1,3)) < scltol .and.                          &
     &      abs(G(2,3) + (G(2,2)+G(1,2)) / 2) < scltol) then
          G(2,2) = (G(2,2) + G(3,3)) / 2
          G(3,3) = G(2,2)
          G(1,2) =-G(1,1)/3
          G(2,1) = G(1,2)
          G(3,1) = G(1,2)
          G(1,3) = G(1,2)
          G(2,3) =-(G(2,2)+G(1,2)) / 2
          G(3,2) = G(2,3)
!dbg           write(6,*) 'rombohedric-4'
          rede=71
          bravais=11
        EndIf

!     SEARCHES FOR HEXAGONAL LATTICES

        If( bravais == 0 .and.                                           &
     &      abs(G(1,1) - G(2,2)) < scltol .and.                          &
     &      abs(G(1,2) + G(1,1) / 2) < scltol .and.                      &
     &      abs(G(2,3)) < scltol .and.                                   &
     &      abs(G(1,3)) < scltol ) then
          G(1,1) = (G(1,1) + G(2,2)) / 2
          G(2,2) = G(1,1)
          G(1,2) = -G(1,1) / 2
          G(2,1) = G(1,2)
          G(2,3) = ZERO
          G(3,2) = ZERO
          G(1,3) = ZERO
          G(3,1) = ZERO
!dbg           write(6,*) 'simple hexagonal-1'
          rede=38
          bravais=10
        EndIf

       If( bravais == 0 .and.                                            &
     &      abs(G(3,3) - G(2,2)) < scltol .and.                          &
     &      abs(G(2,3) + G(2,2) / 2) < scltol .and.                      &
     &      abs(G(1,3)) < scltol .and.                                   &
     &      abs(G(1,2)) < scltol ) then
          G(3,3) = (G(3,3) + G(2,2)) / 2
          G(2,2) = G(3,3)
          G(2,3) = -G(2,2) / 2
          G(3,2) = G(2,3)
          G(1,3) = ZERO
          G(3,1) = ZERO
          G(1,2) = ZERO
          G(2,1) = ZERO
!dbg           write(6,*) 'simple hexagonal-3'
           rede=64
           bravais=10
         EndIf

       If( bravais == 0 .and.                                            &
     &      abs(G(3,3) - G(1,1)) < scltol .and.                          &
     &      abs(G(1,3) + G(1,1) / 2) < scltol .and.                      &
     &      abs(G(2,3)) < scltol .and.                                   &
     &      abs(G(1,2)) < scltol ) then
          G(3,3) = (G(3,3) + G(1,1)) / 2
          G(1,1) = G(3,3)
          G(1,2) = ZERO
          G(2,1) = ZERO
          G(2,3) = ZERO
          G(3,2) = ZERO
          G(1,3) = G(1,1) / 2
          G(3,1) = G(1,3)
!dbg           write(6,*) 'simple hexagonal-2'
          rede=15
          bravais=10
        Endif

!     SEARCHES FOR TETRAGONAL LATTICES

        If( bravais == 0 .and.                                           &
     &      abs(G(1,1) - G(2,2)) < scltol .and.                          &
     &      abs(G(1,2)) < scltol .and.                                   &
     &      abs(G(2,3)) < scltol .and.                                   &
     &      abs(G(1,3)) < scltol ) then
          G(1,1) = (G(1,1) + G(2,2)) / 2
          G(2,2) = G(1,1)
          G(1,2) = ZERO
          G(2,1) = ZERO
          G(2,3) = ZERO
          G(3,2) = ZERO
          G(1,3) = ZERO
          G(3,1) = ZERO
!dbg           write(6,*) 'simple tetragonal-1'
          rede=41
          bravais=4
        EndIf

        If( bravais == 0 .and.                                           &
     &      abs(G(3,3) - G(2,2)) < scltol .and.                          &
     &      abs(G(1,2)) < scltol .and.                                   &
     &      abs(G(2,3)) < scltol .and.                                   &
     &      abs(G(1,3)) < scltol ) then
          G(3,3) = (G(3,3) + G(2,2)) / 2
          G(2,2) = G(3,3)
          G(1,2) = ZERO
          G(2,1) = ZERO
          G(2,3) = ZERO
          G(3,2) = ZERO
          G(1,3) = ZERO
          G(3,1) = ZERO
!dbg           write(6,*) 'simple tetragonal-2'
          rede=65
          bravais=4
        EndIf

        If( bravais == 0 .and.                                           &
     &      abs(G(3,3) - G(2,2)) < scltol .and.                          &
     &      abs(G(1,2) - G(1,1) / 2) < scltol .and.                      &
     &      abs(G(1,3) - G(1,1) / 2) < scltol .and.                      &
     &      abs(G(2,3) - G(1,3) / 2) < scltol ) then
          G(3,3) = (G(3,3) + G(2,2)) / 2
          G(2,2) = G(3,3)
          G(1,2) = G(1,1) / 2
          G(2,1) = G(1,2)
          G(1,3) = G(1,1) / 2
          G(3,1) = G(1,3)
          G(2,3) = G(1,3) / 2
          G(3,2) = G(2,3)
!dbg                  write(6,*) 'centred tetragonal-4'
          rede=2
          bravais=5
        EndIf

        If( bravais == 0 .and.                                           &
     &     abs(G(1,1)-G(2,2)) < scltol .and.                             &
     &     abs(G(2,2)-G(3,3)) < scltol .and.                             &
     &     abs(G(1,2) - G(1,3)) < scltol .and.                           &
     &     abs(G(1,1) + G(2,3) + G(1,2) + G(1,3)) < scltol ) then
          G(1,1) = (G(1,1)+G(2,2)+G(3,3)) / 3
          G(2,2) = G(1,1)
          G(3,3) = G(1,1)
          G(1,2) = (G(1,3)+G(1,2)) / 2
          G(2,1) = G(1,2)
          G(1,3) = G(1,2)
          G(3,1) = G(1,2)
          G(2,3) = -G(1,1) - G(1,2) - G(1,3)
          G(3,2) = G(2,3)
!dbg           write(6,*) 'centred tetragonal-2'
          rede=12
          bravais=5
        Endif

        If( bravais == 0 .and.                                           &
     &     abs(G(1,1)-G(2,2)) < scltol .and.                             &
     &     abs(G(2,2)-G(3,3)) < scltol .and.                             &
     &     abs(G(1,3) - G(2,3)) < scltol .and.                           &
     &     abs(G(1,1) + G(2,3) + G(1,2) + G(1,3)) < scltol ) then
          G(1,1) = (G(1,1)+G(2,2)+G(3,3)) / 3
          G(2,2) = G(1,1)
          G(3,3) = G(1,1)
          G(1,3) = (G(1,3)+G(2,3)) / 2
          G(3,1) = G(1,3)
          G(2,3) = G(1,3)
          G(3,2) = G(1,3)
          G(1,2) = -G(1,1) - G(2,3) - G(1,3)
          G(2,1) = G(1,2)
!dbg           write(6,*) 'centred tetragonal-3'
          rede=13
          bravais=5
        Endif

        If( bravais == 0 .and.                                           &
     &     abs(G(1,1) - G(2,2)) < scltol .and.                           &
     &     abs(G(1,2)) < scltol .and.                                    &
     &     abs(G(1,3) + G(1,1) / 2) < scltol .and.                       &
     &      abs(G(2,3) + G(2,2) / 2) < scltol) then 
          G(1,1) = (G(1,1) + G(2,2)) / 2
          G(2,2) = G(1,1)
          G(1,2) = ZERO
          G(2,1) = ZERO
          G(1,3) = -G(1,1) / 2
          G(3,1) = G(1,3)
          G(2,3) = -G(2,2) / 2
          G(3,2) = G(2,3)
!dbg           write(6,*) 'centred tetragonal-1'
          rede=35
          bravais=5
        Endif

        If( bravais == 0 .and.                                           &
     &     abs(G(2,2) - G(3,3)) < scltol .and.                           &
     &     abs(G(1,2) - G(1,3)) < scltol .and.                           &
     &     abs(G(1,3) - G(1,1) / 2) < scltol .and.                       &
     &     abs(G(2,3) - G(1,3) / 2) < scltol) then 
          G(3,3) = (G(3,3) + G(2,2)) / 2
          G(2,2) = G(3,3)
          G(1,3) = G(1,1) / 2
          G(3,1) = G(1,3)
          G(1,2) = G(1,3)
          G(2,1) = G(1,3)
          G(2,3) = G(1,3) / 2
          G(3,2) = G(2,3)
!dbg           write(6,*) 'centred tetragonal-5'
          rede=49
          bravais=5
        EndIf

!     SEARCHES FOR ORTHORHOMBIC LATTICES

        If( bravais == 0 .and.                                           &
     &     abs(G(1,2)) < scltol .and.                                    &
     &     abs(G(1,3)) < scltol .and.                                    &
     &     abs(G(2,3)) < scltol) then 
          G(1,2) = ZERO
          G(2,1) = ZERO
          G(1,3) = ZERO
          G(3,1) = ZERO
          G(2,3) = ZERO
          G(3,2) = ZERO
!dbg           write(6,*) 'simple orthorhombic'
          rede=95
          bravais=6
        EndIf

        If( bravais == 0 .and.                                           &
     &     abs(G(1,1) - G(2,2)) < scltol .and.                           &
     &     abs(G(1,3)) < scltol .and.                                    &
     &     abs(G(2,3)) < scltol) then 
          G(1,1) = (G(1,1) + G(2,2)) / 2
          G(2,2) = G(1,1)
          G(1,3) = ZERO
          G(3,1) = ZERO
          G(2,3) = ZERO
          G(3,2) = ZERO
!dbg           write(6,*) 'base centred orthorhombic-2'
          rede=20
          bravais=7
        EndIf

        If( bravais == 0 .and.                                           &
     &     abs(G(1,2) + G(1,1) / 2) < scltol .and.                       &
     &     abs(G(1,3)) < scltol .and.                                    &
     &     abs(G(2,3)) < scltol) then 
          G(1,2) = -G(1,1) / 2
          G(2,1) = G(1,2)
          G(1,3) = ZERO
          G(3,1) = ZERO
          G(2,3) = ZERO
          G(3,2) = ZERO
!dbg           write(6,*) 'base centred orthorhombic-4'
          rede=62
          bravais=7
        EndIf

        If( bravais == 0 .and.                                           &
     &     abs(G(3,3) - G(2,2)) < scltol .and.                           &
     &     abs(G(1,3)) < scltol .and.                                    &
     &     abs(G(2,1)) < scltol) then 
          G(3,3) = (G(3,3) + G(2,2)) / 2
          G(2,2) = G(3,3)
          G(1,3) = ZERO
          G(3,1) = ZERO
          G(2,1) = ZERO
          G(1,2) = ZERO
!dbg           write(6,*) 'base centred orthorhombic-3'
          rede=69
          bravais=7
        EndIf

        If( bravais == 0 .and.                                           &
     &     abs(G(2,3) + G(2,2) / 2) < scltol .and.                       &
     &     abs(G(1,3)) < scltol .and.                                    &
     &     abs(G(1,2)) < scltol) then 
          G(2,3) = -G(2,2) / 2
          G(3,2) = G(2,3)
          G(1,3) = ZERO
          G(3,1) = ZERO
          G(1,2) = ZERO
          G(2,1) = ZERO
!dbg           write(6,*) 'base centred orthorhombic-8'
          rede=91
          bravais=7
        EndIf

        If( bravais == 0 .and.                                           &
     &     abs(G(1,3) + G(1,1) / 2) < scltol .and.                       &
     &     abs(G(2,3)) < scltol .and.                                    &
     &     abs(G(1,2)) < scltol) then 
          G(1,3) = -G(1,1) / 2
          G(3,1) = G(1,3)
          G(2,3) = ZERO
          G(3,2) = ZERO
          G(1,2) = ZERO
          G(2,1) = ZERO
!dbg           write(6,*) 'base centred orthorhombic-7'
          rede=93
          bravais=7
        EndIf

        If( bravais == 0 .and.                                           &
     &     abs(G(1,2) - G(1,1) / 2) < scltol .and.                       &
     &     abs(G(1,3) - G(1,1) / 2) < scltol .and.                       &
     &     abs(G(2,3) - G(1,3) / 2) < scltol ) then
          G(1,2) = G(1,1) / 2
          G(2,1) = G(1,2)
          G(1,3) = G(1,1) / 2
          G(3,1) = G(1,3)
          G(2,3) = G(1,3) / 2
          G(3,2) = G(2,3)
!dbg           write(6,*) 'face centred orthorhombic-1'
          rede=25
          bravais=8
        EndIf

        If( bravais == 0 .and.                                           &
     &     abs(G(1,1) - G(2,2)) < scltol .and.                           &
     &     abs(G(1,3) - G(2,3)) < scltol .and.                           &
     &     abs(G(1,2)+G(1,3)+G(2,3) + G(1,1)) < scltol ) then
          G(1,1) = (G(1,1) + G(2,2)) / 2
          G(2,2) = G(1,1)
          G(1,3) = (G(1,3) + G(2,3)) / 2
          G(3,1) = G(1,3)
          G(2,3) = G(1,3)
          G(3,2) = G(1,3)
          G(1,2) = -(G(1,3) + G(2,3) + G(1,1))
          G(2,1) = G(1,2)
!dbg           write(6,*) 'face centred orthorhombic-2'
          rede=36
          bravais=8
        EndIf

        If( bravais == 0 .and.                                           &
     &     abs(G(3,3) - G(2,2)) < scltol .and.                           &
     &     abs(G(1,2) - G(1,1) / 2) < scltol .and.                       &
     &     abs(G(1,3) - G(1,1) / 2) < scltol ) then
          G(3,3) = (G(3,3) + G(2,2)) / 2
          G(2,2) = G(3,3)
          G(1,2) = G(1,1) / 2
          G(2,1) = G(1,2)
          G(1,3) = G(1,1) / 2
          G(3,1) = G(1,3)
!dbg           write(6,*) 'body centred orthorombic-4'
          rede=3
          bravais=9
        EndIf

        If( bravais == 0 .and.                                           &
     &     abs(G(1,1) - G(2,2)) < scltol .and.                           &
     &     abs(G(2,2) - G(3,3)) < scltol .and.                           &
     &     abs(G(1,2)+G(1,3)+G(2,3) + G(1,1)) < scltol ) then
          G(1,1) = (G(1,1)+G(2,2)+G(3,3)) / 3
          G(2,2) = G(1,1)
          G(3,3) = G(1,1)
          G(1,2) = -(G(1,3) + G(2,3) + G(1,1))
          G(2,1) = G(1,2)
!dbg           write(6,*) 'body centred orthorhombic-5'
          rede=14
          bravais=9
        EndIf

        If( bravais == 0 .and.                                           &
     &     abs(G(1,3) + G(1,1) / 2) < scltol .and.                       &
     &     abs(G(2,3) + G(2,2) / 2) < scltol .and.                       &
     &     abs(G(1,2)) < scltol ) then
          G(1,3) =-G(1,1) / 2
          G(3,1) = G(1,3)
          G(2,3) =-G(2,2) / 2
          G(3,2) = G(2,3)
          G(1,2) = ZERO
          G(2,1) = ZERO
!dbg           write(6,*) 'body centred orthorhombic-6'
          rede=87
          bravais=9
        EndIf

!     SEARCHES FOR SIMPLE MONOCLINIC LATTICES

        If( bravais == 0 .and.                                           &
     &     abs(G(1,2)) < scltol .and.                                    &
     &     abs(G(2,3)) < scltol) then 
          G(1,2) = ZERO
          G(2,1) = ZERO
          G(2,3) = ZERO
          G(3,2) = ZERO
!dbg           write(6,*) 'simple monoclinic-1'
          rede=42
          bravais=12
        EndIf

        If( bravais == 0 .and.                                           &
     &     abs(G(2,3)) < scltol .and.                                    &
     &     abs(G(1,3)) < scltol) then 
          G(1,3) = ZERO
          G(3,1) = ZERO
          G(2,3) = ZERO
          G(3,2) = ZERO
!dbg           write(6,*) 'simple monoclinic-2'
          rede=66
          bravais=12
        EndIf

        If( bravais == 0 .and.                                           &
     &     abs(G(1,2)) < scltol .and.                                    &
     &     abs(G(1,3)) < scltol) then 
          G(1,2) = ZERO
          G(2,1) = ZERO
          G(1,3) = ZERO
          G(3,1) = ZERO
!dbg           write(6,*) 'simple monoclinic-3'
          rede=96
          bravais=12
        EndIf

!     SEARCHES FOR CENTERED MONOCLINIC LATTICES

        If( bravais == 0 .and.                                           &
     &     abs(G(1,1)-G(2,2)) < scltol .and.                             &
     &     abs(G(1,3)-G(2,3)) < scltol) then 
          G(1,1) = (G(1,1) + G(2,2)) / 2
          G(2,2) = G(1,1)
          G(1,3) = (G(1,3) + G(2,3)) / 2
          G(3,1) = G(1,3)
          G(2,3) = G(1,3)
          G(3,2) = G(1,3)
!dbg           write(6,*) 'centred monoclinic-1'
          rede=4
          bravais=13
        EndIf

        If( bravais == 0 .and.                                           &
     &     abs(G(3,3)-G(2,2)) < scltol .and.                             &
     &     abs(G(1,3)-G(1,2)) < scltol) then 
          G(3,3) = (G(3,3) + G(2,2)) / 2
          G(2,2) = G(3,3)
          G(1,3) = (G(1,3) + G(1,2)) / 2
          G(3,1) = G(1,3)
          G(1,2) = G(1,3)
          G(2,1) = G(1,3)
!dbg           write(6,*) 'centred monoclinic-2'
          rede=5
          bravais=13
        EndIf

        If( bravais == 0 .and.                                           &
     &     abs(G(1,2) + G(1,1) / 2) < scltol .and.                       &
     &     abs(G(1,3) + G(1,1) / 2) < scltol) then 
          G(1,2) =-G(1,1) / 2
          G(2,1) = G(1,2)
          G(1,3) =-G(1,1) / 2
          G(3,1) = G(1,3)
!dbg           write(6,*) 'centred monoclinic-3'
          rede=8
          bravais=13
        EndIf

        If( bravais == 0 .and.                                           &
     &     abs(G(1,2) - G(1,1) / 2) < scltol .and.                       &
     &     abs(G(1,3) - G(1,1) / 2) < scltol) then 
          G(1,2) = G(1,1) / 2
          G(2,1) = G(1,2)
          G(1,3) = G(1,1) / 2
          G(3,1) = G(1,3)
!dbg           write(6,*) 'centred monoclinic-4'
          rede=9
          bravais=13
        EndIf

        If( bravais == 0 .and.                                           &
     &     abs(G(1,1) - G(2,2)) < scltol .and.                           &
     &     abs(G(1,2)+G(1,3)+G(2,3) + G(1,1)) < scltol ) then
          G(1,1) = (G(1,1) + G(2,2)) / 2
          G(2,2) = G(1,1)
          G(1,2) = -(G(1,3) + G(2,3) + G(1,1))
          G(2,1) = G(1,2)
!dbg           write(6,*) 'centred monoclinic-5'
          rede=18
          bravais=13
        EndIf

        If( bravais == 0 .and.                                           &
     &     abs(G(1,2) - G(1,1) / 2) < scltol .and.                       &
     &     abs(G(1,3)) < scltol) then 
          G(1,2) = G(1,1) / 2
          G(2,1) = G(1,2)
          G(1,3) = ZERO
          G(3,1) = ZERO
!dbg           write(6,*) 'centred monoclinic-6'
          rede=19
          bravais=13
        EndIf

        If( bravais == 0 .and.                                           &
     &     abs(G(1,3) - G(1,1) / 2) < scltol .and.                       &
     &     abs(G(1,2)) < scltol) then 
          G(1,3) = G(1,1) / 2
          G(3,1) = G(1,3)
          G(1,2) = ZERO
          G(2,1) = ZERO
!dbg           write(6,*) 'centred monoclinic-7'
          rede=21
          bravais=13
        EndIf

        If( bravais == 0 .and.                                           &
     &     abs(G(1,2) + G(1,1) / 2) < scltol .and.                       &
     &     abs(G(1,3)) < scltol) then 
          G(1,2) =-G(1,1) / 2
          G(2,1) = G(1,2)
          G(1,3) = ZERO
          G(3,1) = ZERO
!dbg           write(6,*) 'centred monoclinic-8'
          rede=26
          bravais=13
        EndIf

        If( bravais == 0 .and.                                           &
     &     abs(G(1,2) - G(1,1) / 2) < scltol .and.                       &
     &     abs(G(2,3) - G(1,3) / 2) < scltol) then 
          G(1,2) = G(1,1) / 2
          G(2,1) = G(1,2)
          G(2,3) = G(1,3) / 2
          G(3,2) = G(2,3)
!dbg           write(6,*) 'centred monoclinic-9'
          rede=27
          bravais=13
        EndIf

        If( bravais == 0 .and.                                           &
     &     abs(G(2,3) + G(2,2) / 2) < scltol .and.                       &
     &     abs(G(1,2)) < scltol) then 
          G(2,3) =-G(2,2) / 2
          G(3,2) = G(2,3)
          G(1,2) = ZERO
          G(2,1) = ZERO
!dbg           write(6,*) 'centred monoclinic-10'
          rede=31
          bravais=13
        EndIf

        If( bravais == 0 .and.                                           &
     &     abs(G(2,3) - G(2,2) / 2) < scltol .and.                       &
     &     abs(G(1,3) - G(1,2) / 2) < scltol) then 
          G(2,3) = G(2,2) / 2
          G(3,2) = G(2,3)
          G(1,3) = G(1,2) / 2
          G(3,1) = G(1,3)
!dbg           write(6,*) 'centred monoclinic-11'
          rede=33
          bravais=13
        EndIf

        If( bravais == 0 .and.                                           &
     &     abs(G(1,3) - G(1,1) / 2) < scltol .and.                       &
     &     abs(G(2,3) - G(1,2) / 2) < scltol) then 
          G(1,3) = G(1,1) / 2
          G(3,1) = G(1,3)
          G(2,3) = G(1,2) / 2
          G(3,2) = G(2,3)
!dbg           write(6,*) 'centred monoclinic-12'
          rede=37
          bravais=13
        EndIf

        If( bravais == 0 .and.                                           &
     &     abs(G(2,3) - G(2,2) / 2) < scltol .and.                       &
     &     abs(G(1,2)) < scltol) then 
          G(2,3) = G(2,2) / 2
          G(3,2) = G(2,3)
          G(1,2) = ZERO
          G(2,1) = ZERO
!dbg           write(6,*) 'centred monoclinic-13'
          rede=40
          bravais=13
        EndIf

        If( bravais == 0 .and.                                           &
     &     abs(G(1,3) + G(1,1) / 2) < scltol .and.                       &
     &     abs(G(1,2)) < scltol) then 
          G(1,3) =-G(1,1) / 2
          G(3,1) = G(1,3)
          G(1,2) = ZERO
          G(2,1) = ZERO
!dbg           write(6,*) 'centred monoclinic-14'
          rede=43
          bravais=13
        EndIf

        If( bravais == 0 .and.                                           &
     &     abs(G(1,3) + G(1,2) / 2 + G(1,1) / 2) < scltol .and.          &
     &     abs(G(2,3) + G(1,2) / 2 + G(2,2) / 2) < scltol) then 
          G(1,3) =-( G(1,2) + G(1,1) ) / 2
          G(3,1) = G(1,3)
          G(2,3) =-( G(1,2) + G(2,2) ) / 2
          G(3,2) = G(2,3)
!dbg           write(6,*) 'centred monoclinic-15'
          rede=47
          bravais=13
        EndIf

!     WHATEVER IS LEFT IS SHOULD BE TRICLINIC

      If( bravais == 0) then
!dbg            write(6,*) 'triclinic-1'
        rede=6
        bravais=14
      ENDIF

!     recovers the total transformation


      call metric_Tm(2,G,Mtot2)

      call metric_dotM(Mtotal,Mtot2,Mtot1)



!     makes sure that the determinant is 1
      
       invm(1,1) = Mtotal(2,2)*Mtotal(3,3) - Mtotal(3,2)*Mtotal(2,3)
       invm(1,2) = Mtotal(3,2)*Mtotal(1,3) - Mtotal(1,2)*Mtotal(3,3)
       invm(1,3) = Mtotal(1,2)*Mtotal(2,3) - Mtotal(2,2)*Mtotal(1,3)
       invm(2,1) = Mtotal(2,3)*Mtotal(3,1) - Mtotal(3,3)*Mtotal(2,1)
       invm(2,2) = Mtotal(3,3)*Mtotal(1,1) - Mtotal(1,3)*Mtotal(3,1)
       invm(2,3) = Mtotal(1,3)*Mtotal(2,1) - Mtotal(2,3)*Mtotal(1,1)
       invm(3,1) = Mtotal(2,1)*Mtotal(3,2) - Mtotal(3,1)*Mtotal(2,2)
       invm(3,2) = Mtotal(3,1)*Mtotal(1,2) - Mtotal(1,1)*Mtotal(3,2)
       invm(3,3) = Mtotal(1,1)*Mtotal(2,2) - Mtotal(2,1)*Mtotal(1,2)
      
       icount = invm(1,1)*Mtotal(1,1) + invm(1,2)*Mtotal(2,1) +          &
     &          invm(1,3)*Mtotal(3,1)

      if( icount == -1) then

        do i=1,3
        do j=1,3
          Mtotal(i,j) = -Mtotal(i,j)
        enddo
        enddo

      elseif (icount /= 1) then
    
        write(6,*) ' metric_ident'
        write(6,*) ' determinant of Mtotal is ', icount

        stop
      endif

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C                                                                     C
!C  STEP 5: GENERATES CONVENTIONAL LATTICE VECTORS AND PRINTS          C
!C                                                                     C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      call Metric_Print(G,rede,avecnig,aconv,iprint)
      
      do i=1,3
      do j=1,3
        avec(i,j) = ZERO
        do k=1,3
          avec(i,j) = avec(i,j) + invm(j,k)*avecnig(i,k)
        enddo
      enddo
      enddo
      
      do i=1,3
      do j=1,3
        Gout(i,j) = G(i,j)
      enddo
      enddo

      call metric_Tm(1,Gout,invm)

!     checks if all is correct

      error = ZERO
      do i=1,3
      do j=1,3
        error = max(error,abs(Gin(i,j)-Gout(i,j)))
      enddo
      enddo

      if( error > 5.0*skew*scltol) then

        write(6,*) ' metric_ident'
        write(6,*) ' Gin ' , ((Gin(i,j),i=1,3),j=1,3)
        write(6,*) ' Gout ' , ((Gout(i,j),i=1,3),j=1,3)

        write(6,*) error,scltol,skew,'error,scltol'

        stop

      endif 

      do i=1,3
      do j=1,3
        Gavec(i,j) = ZERO
        do k=1,3
          Gavec(i,j) = Gavec(i,j) + avec(k,j)*avec(k,i)
        enddo
      enddo
      enddo

      error = ZERO
      do i=1,3
      do j=1,3
        error = max(error,abs(Gavec(i,j)-Gout(i,j)))
      enddo
      enddo

      if( error > 5.0*skew*scltol) then

        write(6,*) ' metric_ident'
        write(6,*) ' Gavec ' , ((Gavec(i,j),i=1,3),j=1,3)
        write(6,*) ' Gout ' , ((Gout(i,j),i=1,3),j=1,3)

        write(6,*) error,scltol,skew,'error,scltol'

        stop

      endif 

      return
      end subroutine metric_ident

