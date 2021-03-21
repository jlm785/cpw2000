       SUBROUTINE VPASS2(A,B,C,D,TRIGS,TRIGI,IFAC,LA,N,
     1  ISKIP,ISKIP1,ISKIP3,NUMBER,ISIGN)
C
C      WRITTEN JULY 24, 1992 BY NORMAN J. TROULLIER JR.
C      COPYRIGHT NORMAN J. TROULLIER JR.
C      DISTRIBUTION ONLY BY WRITTEN PERMISSION
C      NOTE: *vdir nodep ARE NEC-SX3 DIRECTIVES
C
       IMPLICIT REAL*8 (A-H,O-Z)
C
       DIMENSION A(*),C(*),TRIGS(*),TRIGI(*),B(*),D(*)
C
       PARAMETER (SIN60=0.8660254037844386D0)
       PARAMETER (SIN72=0.9510565162951536D0)
       PARAMETER (SIN36=0.5877852522924732D0)
       PARAMETER (SQ54=0.5590169943749475D0)
       PARAMETER (PO25=0.25D0)
       PARAMETER (SQ2=0.7071067811865476D0)
       PARAMETER (PO5=0.5D0)
       PARAMETER (SIN20=0.3420201433256687D0)
       PARAMETER (SIN40=0.6427876096865393D0)
       PARAMETER (SIN80=0.9848077530122080D0)
       PARAMETER (COS20=0.9396926207859083D0)
       PARAMETER (COS40=0.7660444431189781D0)
       PARAMETER (COS80=0.1736481776669304D0)
       PARAMETER (COS36=0.8090169943749475D0)
       PARAMETER (COS72=0.3090169943749475D0)
C
       M=N/IFAC
       MSKIP=M*ISKIP
       IB=MSKIP
       LASKIP=LA*ISKIP1
       JB=LASKIP
       I=1
       J=1
       JUMP=(IFAC-1)*LASKIP
C
       IF(IFAC .EQ. 10) THEN
C
C        FACTOR 10
C
         IC=IB+MSKIP
         ID=IC+MSKIP
         IE=ID+MSKIP
         IF=IE+MSKIP
         IG=IF+MSKIP
         IH=IG+MSKIP
         II=IH+MSKIP
         IJ=II+MSKIP
         JC=JB+LASKIP
         JD=JC+LASKIP
         JE=JD+LASKIP
         JF=JE+LASKIP
         JG=JF+LASKIP
         JH=JG+LASKIP
         JI=JH+LASKIP
         JJ=JI+LASKIP
         IF(ISIGN.LT.0) THEN
           DO L=1,LA
             IADD2=-ISKIP3+J
*vdir nodep             
CDIR$ IVDEP
CDEC$ IVDEP
             DO NU=0,NUMBER-1
               IADD1=NU+I
               IADD2=IADD2+ISKIP3
               T1REAL=A(IADD1)+A(IADD1+IF)
               T1IMAG=B(IADD1)+B(IADD1+IF)
               T2REAL=A(IADD1)-A(IADD1+IF)
               T2IMAG=B(IADD1)-B(IADD1+IF)
               T3REAL=A(IADD1+IB)+A(IADD1+IG)
               T3IMAG=B(IADD1+IB)+B(IADD1+IG)
               T4REAL=A(IADD1+IB)-A(IADD1+IG)
               T4IMAG=B(IADD1+IB)-B(IADD1+IG)
               T5REAL=A(IADD1+IC)+A(IADD1+IH)
               T5IMAG=B(IADD1+IC)+B(IADD1+IH)
               T6REAL=A(IADD1+IC)-A(IADD1+IH)
               T6IMAG=B(IADD1+IC)-B(IADD1+IH)
               T7REAL=A(IADD1+ID)+A(IADD1+II)
               T7IMAG=B(IADD1+ID)+B(IADD1+II)
               T8REAL=A(IADD1+ID)-A(IADD1+II)
               T8IMAG=B(IADD1+ID)-B(IADD1+II)
               T9REAL=A(IADD1+IE)+A(IADD1+IJ)
               T9IMAG=B(IADD1+IE)+B(IADD1+IJ)
               T10REAL=A(IADD1+IE)-A(IADD1+IJ)
               T10IMAG=B(IADD1+IE)-B(IADD1+IJ)
               B1REAL=T3REAL+T9REAL
               B1IMAG=T3IMAG+T9IMAG
               B2REAL=T3REAL-T9REAL
               B2IMAG=T3IMAG-T9IMAG
               B3REAL=T5REAL+T7REAL
               B3IMAG=T5IMAG+T7IMAG
               B4REAL=T5REAL-T7REAL
               B4IMAG=T5IMAG-T7IMAG
               B5REAL=T4REAL+T10REAL
               B5IMAG=T4IMAG+T10IMAG
               B6REAL=T4REAL-T10REAL
               B6IMAG=T4IMAG-T10IMAG
               B7REAL=T6REAL+T8REAL
               B7IMAG=T6IMAG+T8IMAG
               B8REAL=T6REAL-T8REAL
               B8IMAG=T6IMAG-T8IMAG
               C1REAL=T2REAL+B6REAL*COS36+B8REAL*COS72
               C1IMAG=T2IMAG+B6IMAG*COS36+B8IMAG*COS72
               C2REAL=B5REAL*SIN36+B7REAL*SIN72
               C2IMAG=B5IMAG*SIN36+B7IMAG*SIN72
               C3REAL=T1REAL+B1REAL*COS72-B3REAL*COS36
               C3IMAG=T1IMAG+B1IMAG*COS72-B3IMAG*COS36
               C4REAL=B2REAL*SIN72+B4REAL*SIN36
               C4IMAG=B2IMAG*SIN72+B4IMAG*SIN36
               C5REAL=T2REAL-B6REAL*COS72-B8REAL*COS36
               C5IMAG=T2IMAG-B6IMAG*COS72-B8IMAG*COS36
               C6REAL=B5REAL*SIN72-B7REAL*SIN36
               C6IMAG=B5IMAG*SIN72-B7IMAG*SIN36
               C7REAL=T1REAL-B1REAL*COS36+B3REAL*COS72
               C7IMAG=T1IMAG-B1IMAG*COS36+B3IMAG*COS72
               C8REAL=B2REAL*SIN36-B4REAL*SIN72
               C8IMAG=B2IMAG*SIN36-B4IMAG*SIN72
               C(IADD2)   =T1REAL+B1REAL+B3REAL
               D(IADD2)   =T1IMAG+B1IMAG+B3IMAG
               C(IADD2+JB)=C1REAL-C2IMAG
               D(IADD2+JB)=C1IMAG+C2REAL
               C(IADD2+JC)=C3REAL-C4IMAG
               D(IADD2+JC)=C3IMAG+C4REAL
               C(IADD2+JD)=C5REAL-C6IMAG
               D(IADD2+JD)=C5IMAG+C6REAL
               C(IADD2+JE)=C7REAL-C8IMAG
               D(IADD2+JE)=C7IMAG+C8REAL
               C(IADD2+JF)=T2REAL-B6REAL+B8REAL
               D(IADD2+JF)=T2IMAG-B6IMAG+B8IMAG
               C(IADD2+JG)=C7REAL+C8IMAG
               D(IADD2+JG)=C7IMAG-C8REAL
               C(IADD2+JH)=C5REAL+C6IMAG
               D(IADD2+JH)=C5IMAG-C6REAL
               C(IADD2+JI)=C3REAL+C4IMAG
               D(IADD2+JI)=C3IMAG-C4REAL
               C(IADD2+JJ)=C1REAL+C2IMAG
               D(IADD2+JJ)=C1IMAG-C2REAL
             ENDDO
             I=I+ISKIP
             J=J+ISKIP1
           ENDDO
           J=J+JUMP
           DO K=LA,M-LA,LA
             K2=K+K
             KK1=K2+1
             KK2=KK1+K2
             KK3=KK2+K2
             KK4=KK3+K2
             KK5=KK4+K2
             KK6=KK5+K2
             KK7=KK6+K2
             KK8=KK7+K2
             KK9=KK8+K2
             DO L=1,LA
               IADD2=-ISKIP3+J
*vdir nodep             
CDIR$ IVDEP
CDEC$ IVDEP
               DO NU=0,NUMBER-1
                 IADD1=NU+I
                 IADD2=IADD2+ISKIP3
                 T1REAL=A(IADD1)+A(IADD1+IF)
                 T1IMAG=B(IADD1)+B(IADD1+IF)
                 T2REAL=A(IADD1)-A(IADD1+IF)
                 T2IMAG=B(IADD1)-B(IADD1+IF)
                 T3REAL=A(IADD1+IB)+A(IADD1+IG)
                 T3IMAG=B(IADD1+IB)+B(IADD1+IG)
                 T4REAL=A(IADD1+IB)-A(IADD1+IG)
                 T4IMAG=B(IADD1+IB)-B(IADD1+IG)
                 T5REAL=A(IADD1+IC)+A(IADD1+IH)
                 T5IMAG=B(IADD1+IC)+B(IADD1+IH)
                 T6REAL=A(IADD1+IC)-A(IADD1+IH)
                 T6IMAG=B(IADD1+IC)-B(IADD1+IH)
                 T7REAL=A(IADD1+ID)+A(IADD1+II)
                 T7IMAG=B(IADD1+ID)+B(IADD1+II)
                 T8REAL=A(IADD1+ID)-A(IADD1+II)
                 T8IMAG=B(IADD1+ID)-B(IADD1+II)
                 T9REAL=A(IADD1+IE)+A(IADD1+IJ)
                 T9IMAG=B(IADD1+IE)+B(IADD1+IJ)
                 T10REAL=A(IADD1+IE)-A(IADD1+IJ)
                 T10IMAG=B(IADD1+IE)-B(IADD1+IJ)
                 B1REAL=T3REAL+T9REAL
                 B1IMAG=T3IMAG+T9IMAG
                 B2REAL=T3REAL-T9REAL
                 B2IMAG=T3IMAG-T9IMAG
                 B3REAL=T5REAL+T7REAL
                 B3IMAG=T5IMAG+T7IMAG
                 B4REAL=T5REAL-T7REAL
                 B4IMAG=T5IMAG-T7IMAG
                 B5REAL=T4REAL+T10REAL
                 B5IMAG=T4IMAG+T10IMAG
                 B6REAL=T4REAL-T10REAL
                 B6IMAG=T4IMAG-T10IMAG
                 B7REAL=T6REAL+T8REAL
                 B7IMAG=T6IMAG+T8IMAG
                 B8REAL=T6REAL-T8REAL
                 B8IMAG=T6IMAG-T8IMAG
                 C1REAL=T2REAL+B6REAL*COS36+B8REAL*COS72
                 C1IMAG=T2IMAG+B6IMAG*COS36+B8IMAG*COS72
                 C2REAL=B5REAL*SIN36+B7REAL*SIN72
                 C2IMAG=B5IMAG*SIN36+B7IMAG*SIN72
                 C3REAL=T1REAL+B1REAL*COS72-B3REAL*COS36
                 C3IMAG=T1IMAG+B1IMAG*COS72-B3IMAG*COS36
                 C4REAL=B2REAL*SIN72+B4REAL*SIN36
                 C4IMAG=B2IMAG*SIN72+B4IMAG*SIN36
                 C5REAL=T2REAL-B6REAL*COS72-B8REAL*COS36
                 C5IMAG=T2IMAG-B6IMAG*COS72-B8IMAG*COS36
                 C6REAL=B5REAL*SIN72-B7REAL*SIN36
                 C6IMAG=B5IMAG*SIN72-B7IMAG*SIN36
                 C7REAL=T1REAL-B1REAL*COS36+B3REAL*COS72
                 C7IMAG=T1IMAG-B1IMAG*COS36+B3IMAG*COS72
                 C8REAL=B2REAL*SIN36-B4REAL*SIN72
                 C8IMAG=B2IMAG*SIN36-B4IMAG*SIN72
                 X1REAL=C1REAL-C2IMAG
                 X1IMAG=C1IMAG+C2REAL
                 X2REAL=C3REAL-C4IMAG
                 X2IMAG=C3IMAG+C4REAL
                 X3REAL=C5REAL-C6IMAG
                 X3IMAG=C5IMAG+C6REAL
                 X4REAL=C7REAL-C8IMAG
                 X4IMAG=C7IMAG+C8REAL
                 X5REAL=T2REAL-B6REAL+B8REAL
                 X5IMAG=T2IMAG-B6IMAG+B8IMAG
                 X6REAL=C7REAL+C8IMAG
                 X6IMAG=C7IMAG-C8REAL
                 X7REAL=C5REAL+C6IMAG
                 X7IMAG=C5IMAG-C6REAL
                 X8REAL=C3REAL+C4IMAG
                 X8IMAG=C3IMAG-C4REAL
                 X9REAL=C1REAL+C2IMAG
                 X9IMAG=C1IMAG-C2REAL
                 C(IADD2)   =T1REAL+B1REAL+B3REAL
                 D(IADD2)   =T1IMAG+B1IMAG+B3IMAG
                 C(IADD2+JB)=TRIGS(KK1)*X1REAL-TRIGI(KK1)*X1IMAG
                 D(IADD2+JB)=TRIGI(KK1)*X1REAL+TRIGS(KK1)*X1IMAG
                 C(IADD2+JC)=TRIGS(KK2)*X2REAL-TRIGI(KK2)*X2IMAG
                 D(IADD2+JC)=TRIGI(KK2)*X2REAL+TRIGS(KK2)*X2IMAG
                 C(IADD2+JD)=TRIGS(KK3)*X3REAL-TRIGI(KK3)*X3IMAG
                 D(IADD2+JD)=TRIGI(KK3)*X3REAL+TRIGS(KK3)*X3IMAG
                 C(IADD2+JE)=TRIGS(KK4)*X4REAL-TRIGI(KK4)*X4IMAG
                 D(IADD2+JE)=TRIGI(KK4)*X4REAL+TRIGS(KK4)*X4IMAG
                 C(IADD2+JF)=TRIGS(KK5)*X5REAL-TRIGI(KK5)*X5IMAG
                 D(IADD2+JF)=TRIGI(KK5)*X5REAL+TRIGS(KK5)*X5IMAG
                 C(IADD2+JG)=TRIGS(KK6)*X6REAL-TRIGI(KK6)*X6IMAG
                 D(IADD2+JG)=TRIGI(KK6)*X6REAL+TRIGS(KK6)*X6IMAG
                 C(IADD2+JH)=TRIGS(KK7)*X7REAL-TRIGI(KK7)*X7IMAG
                 D(IADD2+JH)=TRIGI(KK7)*X7REAL+TRIGS(KK7)*X7IMAG
                 C(IADD2+JI)=TRIGS(KK8)*X8REAL-TRIGI(KK8)*X8IMAG
                 D(IADD2+JI)=TRIGI(KK8)*X8REAL+TRIGS(KK8)*X8IMAG
                 C(IADD2+JJ)=TRIGS(KK9)*X9REAL-TRIGI(KK9)*X9IMAG
                 D(IADD2+JJ)=TRIGI(KK9)*X9REAL+TRIGS(KK9)*X9IMAG
               ENDDO
               I=I+ISKIP
               J=J+ISKIP1
             ENDDO
             J=J+JUMP
           ENDDO
         ELSE
           DO L=1,LA
             IADD2=-ISKIP3+J
*vdir nodep             
CDIR$ IVDEP
CDEC$ IVDEP
             DO NU=0,NUMBER-1
               IADD1=NU+I
               IADD2=IADD2+ISKIP3
               T1REAL=A(IADD1)+A(IADD1+IF)
               T1IMAG=B(IADD1)+B(IADD1+IF)
               T2REAL=A(IADD1)-A(IADD1+IF)
               T2IMAG=B(IADD1)-B(IADD1+IF)
               T3REAL=A(IADD1+IB)+A(IADD1+IG)
               T3IMAG=B(IADD1+IB)+B(IADD1+IG)
               T4REAL=A(IADD1+IB)-A(IADD1+IG)
               T4IMAG=B(IADD1+IB)-B(IADD1+IG)
               T5REAL=A(IADD1+IC)+A(IADD1+IH)
               T5IMAG=B(IADD1+IC)+B(IADD1+IH)
               T6REAL=A(IADD1+IC)-A(IADD1+IH)
               T6IMAG=B(IADD1+IC)-B(IADD1+IH)
               T7REAL=A(IADD1+ID)+A(IADD1+II)
               T7IMAG=B(IADD1+ID)+B(IADD1+II)
               T8REAL=A(IADD1+ID)-A(IADD1+II)
               T8IMAG=B(IADD1+ID)-B(IADD1+II)
               T9REAL=A(IADD1+IE)+A(IADD1+IJ)
               T9IMAG=B(IADD1+IE)+B(IADD1+IJ)
               T10REAL=A(IADD1+IE)-A(IADD1+IJ)
               T10IMAG=B(IADD1+IE)-B(IADD1+IJ)
               B1REAL=T3REAL+T9REAL
               B1IMAG=T3IMAG+T9IMAG
               B2REAL=T3REAL-T9REAL
               B2IMAG=T3IMAG-T9IMAG
               B3REAL=T5REAL+T7REAL
               B3IMAG=T5IMAG+T7IMAG
               B4REAL=T5REAL-T7REAL
               B4IMAG=T5IMAG-T7IMAG
               B5REAL=T4REAL+T10REAL
               B5IMAG=T4IMAG+T10IMAG
               B6REAL=T4REAL-T10REAL
               B6IMAG=T4IMAG-T10IMAG
               B7REAL=T6REAL+T8REAL
               B7IMAG=T6IMAG+T8IMAG
               B8REAL=T6REAL-T8REAL
               B8IMAG=T6IMAG-T8IMAG
               C1REAL=T2REAL+B6REAL*COS36+B8REAL*COS72
               C1IMAG=T2IMAG+B6IMAG*COS36+B8IMAG*COS72
               C2REAL=B5REAL*SIN36+B7REAL*SIN72
               C2IMAG=B5IMAG*SIN36+B7IMAG*SIN72
               C3REAL=T1REAL+B1REAL*COS72-B3REAL*COS36
               C3IMAG=T1IMAG+B1IMAG*COS72-B3IMAG*COS36
               C4REAL=B2REAL*SIN72+B4REAL*SIN36
               C4IMAG=B2IMAG*SIN72+B4IMAG*SIN36
               C5REAL=T2REAL-B6REAL*COS72-B8REAL*COS36
               C5IMAG=T2IMAG-B6IMAG*COS72-B8IMAG*COS36
               C6REAL=B5REAL*SIN72-B7REAL*SIN36
               C6IMAG=B5IMAG*SIN72-B7IMAG*SIN36
               C7REAL=T1REAL-B1REAL*COS36+B3REAL*COS72
               C7IMAG=T1IMAG-B1IMAG*COS36+B3IMAG*COS72
               C8REAL=B2REAL*SIN36-B4REAL*SIN72
               C8IMAG=B2IMAG*SIN36-B4IMAG*SIN72
               C(IADD2)   =T1REAL+B1REAL+B3REAL
               D(IADD2)   =T1IMAG+B1IMAG+B3IMAG
               C(IADD2+JB)=C1REAL+C2IMAG
               D(IADD2+JB)=C1IMAG-C2REAL
               C(IADD2+JC)=C3REAL+C4IMAG
               D(IADD2+JC)=C3IMAG-C4REAL
               C(IADD2+JD)=C5REAL+C6IMAG
               D(IADD2+JD)=C5IMAG-C6REAL
               C(IADD2+JE)=C7REAL+C8IMAG
               D(IADD2+JE)=C7IMAG-C8REAL
               C(IADD2+JF)=T2REAL-B6REAL+B8REAL
               D(IADD2+JF)=T2IMAG-B6IMAG+B8IMAG
               C(IADD2+JG)=C7REAL-C8IMAG
               D(IADD2+JG)=C7IMAG+C8REAL
               C(IADD2+JH)=C5REAL-C6IMAG
               D(IADD2+JH)=C5IMAG+C6REAL
               C(IADD2+JI)=C3REAL-C4IMAG
               D(IADD2+JI)=C3IMAG+C4REAL
               C(IADD2+JJ)=C1REAL-C2IMAG
               D(IADD2+JJ)=C1IMAG+C2REAL
             ENDDO
             I=I+ISKIP
             J=J+ISKIP1
           ENDDO
           J=J+JUMP
           DO K=LA,M-LA,LA
             K2=K+K
             KK1=K2+1
             KK2=KK1+K2
             KK3=KK2+K2
             KK4=KK3+K2
             KK5=KK4+K2
             KK6=KK5+K2
             KK7=KK6+K2
             KK8=KK7+K2
             KK9=KK8+K2
             DO L=1,LA
               IADD2=-ISKIP3+J
*vdir nodep             
CDIR$ IVDEP
CDEC$ IVDEP
               DO NU=0,NUMBER-1
                 IADD1=NU+I
                 IADD2=IADD2+ISKIP3
                 T1REAL=A(IADD1)+A(IADD1+IF)
                 T1IMAG=B(IADD1)+B(IADD1+IF)
                 T2REAL=A(IADD1)-A(IADD1+IF)
                 T2IMAG=B(IADD1)-B(IADD1+IF)
                 T3REAL=A(IADD1+IB)+A(IADD1+IG)
                 T3IMAG=B(IADD1+IB)+B(IADD1+IG)
                 T4REAL=A(IADD1+IB)-A(IADD1+IG)
                 T4IMAG=B(IADD1+IB)-B(IADD1+IG)
                 T5REAL=A(IADD1+IC)+A(IADD1+IH)
                 T5IMAG=B(IADD1+IC)+B(IADD1+IH)
                 T6REAL=A(IADD1+IC)-A(IADD1+IH)
                 T6IMAG=B(IADD1+IC)-B(IADD1+IH)
                 T7REAL=A(IADD1+ID)+A(IADD1+II)
                 T7IMAG=B(IADD1+ID)+B(IADD1+II)
                 T8REAL=A(IADD1+ID)-A(IADD1+II)
                 T8IMAG=B(IADD1+ID)-B(IADD1+II)
                 T9REAL=A(IADD1+IE)+A(IADD1+IJ)
                 T9IMAG=B(IADD1+IE)+B(IADD1+IJ)
                 T10REAL=A(IADD1+IE)-A(IADD1+IJ)
                 T10IMAG=B(IADD1+IE)-B(IADD1+IJ)
                 B1REAL=T3REAL+T9REAL
                 B1IMAG=T3IMAG+T9IMAG
                 B2REAL=T3REAL-T9REAL
                 B2IMAG=T3IMAG-T9IMAG
                 B3REAL=T5REAL+T7REAL
                 B3IMAG=T5IMAG+T7IMAG
                 B4REAL=T5REAL-T7REAL
                 B4IMAG=T5IMAG-T7IMAG
                 B5REAL=T4REAL+T10REAL
                 B5IMAG=T4IMAG+T10IMAG
                 B6REAL=T4REAL-T10REAL
                 B6IMAG=T4IMAG-T10IMAG
                 B7REAL=T6REAL+T8REAL
                 B7IMAG=T6IMAG+T8IMAG
                 B8REAL=T6REAL-T8REAL
                 B8IMAG=T6IMAG-T8IMAG
                 C1REAL=T2REAL+B6REAL*COS36+B8REAL*COS72
                 C1IMAG=T2IMAG+B6IMAG*COS36+B8IMAG*COS72
                 C2REAL=B5REAL*SIN36+B7REAL*SIN72
                 C2IMAG=B5IMAG*SIN36+B7IMAG*SIN72
                 C3REAL=T1REAL+B1REAL*COS72-B3REAL*COS36
                 C3IMAG=T1IMAG+B1IMAG*COS72-B3IMAG*COS36
                 C4REAL=B2REAL*SIN72+B4REAL*SIN36
                 C4IMAG=B2IMAG*SIN72+B4IMAG*SIN36
                 C5REAL=T2REAL-B6REAL*COS72-B8REAL*COS36
                 C5IMAG=T2IMAG-B6IMAG*COS72-B8IMAG*COS36
                 C6REAL=B5REAL*SIN72-B7REAL*SIN36
                 C6IMAG=B5IMAG*SIN72-B7IMAG*SIN36
                 C7REAL=T1REAL-B1REAL*COS36+B3REAL*COS72
                 C7IMAG=T1IMAG-B1IMAG*COS36+B3IMAG*COS72
                 C8REAL=B2REAL*SIN36-B4REAL*SIN72
                 C8IMAG=B2IMAG*SIN36-B4IMAG*SIN72
                 X1REAL=C1REAL+C2IMAG
                 X1IMAG=C1IMAG-C2REAL
                 X2REAL=C3REAL+C4IMAG
                 X2IMAG=C3IMAG-C4REAL
                 X3REAL=C5REAL+C6IMAG
                 X3IMAG=C5IMAG-C6REAL
                 X4REAL=C7REAL+C8IMAG
                 X4IMAG=C7IMAG-C8REAL
                 X5REAL=T2REAL-B6REAL+B8REAL
                 X5IMAG=T2IMAG-B6IMAG+B8IMAG
                 X6REAL=C7REAL-C8IMAG
                 X6IMAG=C7IMAG+C8REAL
                 X7REAL=C5REAL-C6IMAG
                 X7IMAG=C5IMAG+C6REAL
                 X8REAL=C3REAL-C4IMAG
                 X8IMAG=C3IMAG+C4REAL
                 X9REAL=C1REAL-C2IMAG
                 X9IMAG=C1IMAG+C2REAL
                 C(IADD2)   = T1REAL+B1REAL+B3REAL
                 D(IADD2)   = T1IMAG+B1IMAG+B3IMAG
                 C(IADD2+JB)= TRIGS(KK1)*X1REAL+TRIGI(KK1)*X1IMAG
                 D(IADD2+JB)=-TRIGI(KK1)*X1REAL+TRIGS(KK1)*X1IMAG
                 C(IADD2+JC)= TRIGS(KK2)*X2REAL+TRIGI(KK2)*X2IMAG
                 D(IADD2+JC)=-TRIGI(KK2)*X2REAL+TRIGS(KK2)*X2IMAG
                 C(IADD2+JD)= TRIGS(KK3)*X3REAL+TRIGI(KK3)*X3IMAG
                 D(IADD2+JD)=-TRIGI(KK3)*X3REAL+TRIGS(KK3)*X3IMAG
                 C(IADD2+JE)= TRIGS(KK4)*X4REAL+TRIGI(KK4)*X4IMAG
                 D(IADD2+JE)=-TRIGI(KK4)*X4REAL+TRIGS(KK4)*X4IMAG
                 C(IADD2+JF)= TRIGS(KK5)*X5REAL+TRIGI(KK5)*X5IMAG
                 D(IADD2+JF)=-TRIGI(KK5)*X5REAL+TRIGS(KK5)*X5IMAG
                 C(IADD2+JG)= TRIGS(KK6)*X6REAL+TRIGI(KK6)*X6IMAG
                 D(IADD2+JG)=-TRIGI(KK6)*X6REAL+TRIGS(KK6)*X6IMAG
                 C(IADD2+JH)= TRIGS(KK7)*X7REAL+TRIGI(KK7)*X7IMAG
                 D(IADD2+JH)=-TRIGI(KK7)*X7REAL+TRIGS(KK7)*X7IMAG
                 C(IADD2+JI)= TRIGS(KK8)*X8REAL+TRIGI(KK8)*X8IMAG
                 D(IADD2+JI)=-TRIGI(KK8)*X8REAL+TRIGS(KK8)*X8IMAG
                 C(IADD2+JJ)= TRIGS(KK9)*X9REAL+TRIGI(KK9)*X9IMAG
                 D(IADD2+JJ)=-TRIGI(KK9)*X9REAL+TRIGS(KK9)*X9IMAG
               ENDDO
               I=I+ISKIP
               J=J+ISKIP1
             ENDDO
             J=J+JUMP
           ENDDO
         ENDIF
       ELSEIF(IFAC .EQ.  9) THEN
C
C        FACTOR 9
C
         IC=IB+MSKIP
         ID=IC+MSKIP
         IE=ID+MSKIP
         IF=IE+MSKIP
         IG=IF+MSKIP
         IH=IG+MSKIP
         II=IH+MSKIP
         JC=JB+LASKIP
         JD=JC+LASKIP
         JE=JD+LASKIP
         JF=JE+LASKIP
         JG=JF+LASKIP
         JH=JG+LASKIP
         JI=JH+LASKIP
         IF(ISIGN.LT.0) THEN
           DO L=1,LA
             IADD2=-ISKIP3+J
*vdir nodep             
CDIR$ IVDEP
CDEC$ IVDEP
             DO NU=0,NUMBER-1
               IADD1=NU+I
               IADD2=IADD2+ISKIP3
               T1REAL=A(IADD1+IB)+A(IADD1+II)
               T1IMAG=B(IADD1+IB)+B(IADD1+II)
               T2REAL=A(IADD1+IB)-A(IADD1+II)
               T2IMAG=B(IADD1+IB)-B(IADD1+II)
               T3REAL=A(IADD1+IC)+A(IADD1+IH)
               T3IMAG=B(IADD1+IC)+B(IADD1+IH)
               T4REAL=A(IADD1+IC)-A(IADD1+IH)
               T4IMAG=B(IADD1+IC)-B(IADD1+IH)
               T5REAL=A(IADD1+ID)+A(IADD1+IG)
               T5IMAG=B(IADD1+ID)+B(IADD1+IG)
               T6REAL=A(IADD1+ID)-A(IADD1+IG)
               T6IMAG=B(IADD1+ID)-B(IADD1+IG)
               T7REAL=A(IADD1+IE)+A(IADD1+IF)
               T7IMAG=B(IADD1+IE)+B(IADD1+IF)
               T8REAL=A(IADD1+IE)-A(IADD1+IF)
               T8IMAG=B(IADD1+IE)-B(IADD1+IF)
               B1REAL=A(IADD1)+T5REAL
               B1IMAG=B(IADD1)+T5IMAG
               B2REAL=T1REAL+T3REAL+T7REAL
               B2IMAG=T1IMAG+T3IMAG+T7IMAG
               B3REAL=B1REAL-PO5*B2REAL
               B3IMAG=B1IMAG-PO5*B2IMAG
               B4REAL=(T2REAL-T4REAL+T8REAL)*SIN60
               B4IMAG=(T2IMAG-T4IMAG+T8IMAG)*SIN60
               B5REAL=-PO5*T5REAL+A(IADD1)
               B5IMAG=-PO5*T5IMAG+B(IADD1)
               B6REAL=T6REAL*SIN60
               B6IMAG=T6IMAG*SIN60
               B7REAL=T1REAL*COS40+T3REAL*COS80+B5REAL-T7REAL*COS20
               B7IMAG=T1IMAG*COS40+T3IMAG*COS80+B5IMAG-T7IMAG*COS20
               B8REAL=T2REAL*SIN40+T4REAL*SIN80+B6REAL+T8REAL*SIN20
               B8IMAG=T2IMAG*SIN40+T4IMAG*SIN80+B6IMAG+T8IMAG*SIN20
               B9REAL=T1REAL*COS80-T3REAL*COS20+B5REAL+T7REAL*COS40
               B9IMAG=T1IMAG*COS80-T3IMAG*COS20+B5IMAG+T7IMAG*COS40
               B10REAL=T2REAL*SIN80+T4REAL*SIN20-B6REAL-T8REAL*SIN40
               B10IMAG=T2IMAG*SIN80+T4IMAG*SIN20-B6IMAG-T8IMAG*SIN40
               B11REAL=-T1REAL*COS20+T3REAL*COS40+B5REAL+T7REAL*COS80
               B11IMAG=-T1IMAG*COS20+T3IMAG*COS40+B5IMAG+T7IMAG*COS80
               B12REAL=T2REAL*SIN20-T4REAL*SIN40+B6REAL-T8REAL*SIN80
               B12IMAG=T2IMAG*SIN20-T4IMAG*SIN40+B6IMAG-T8IMAG*SIN80
               C(IADD2)   =B1REAL+B2REAL
               D(IADD2)   =B1IMAG+B2IMAG
               C(IADD2+JB)=B7REAL-B8IMAG
               D(IADD2+JB)=B7IMAG+B8REAL
               C(IADD2+JC)=B9REAL-B10IMAG
               D(IADD2+JC)=B9IMAG+B10REAL
               C(IADD2+JD)=B3REAL-B4IMAG
               D(IADD2+JD)=B3IMAG+B4REAL
               C(IADD2+JE)=B11REAL-B12IMAG
               D(IADD2+JE)=B11IMAG+B12REAL
               C(IADD2+JF)=B11REAL+B12IMAG
               D(IADD2+JF)=B11IMAG-B12REAL
               C(IADD2+JG)=B3REAL+B4IMAG
               D(IADD2+JG)=B3IMAG-B4REAL
               C(IADD2+JH)=B9REAL+B10IMAG
               D(IADD2+JH)=B9IMAG-B10REAL
               C(IADD2+JI)=B7REAL+B8IMAG
               D(IADD2+JI)=B7IMAG-B8REAL 
             ENDDO
             I=I+ISKIP
             J=J+ISKIP1
           ENDDO
           J=J+JUMP
           DO K=LA,M-LA,LA
             K2=K+K
             KK1=K2+1
             KK2=KK1+K2
             KK3=KK2+K2
             KK4=KK3+K2
             KK5=KK4+K2
             KK6=KK5+K2
             KK7=KK6+K2
             KK8=KK7+K2
             DO L=1,LA
               IADD2=-ISKIP3+J
*vdir nodep             
CDIR$ IVDEP
CDEC$ IVDEP
               DO NU=0,NUMBER-1
                 IADD1=NU+I
                 IADD2=IADD2+ISKIP3
                 T1REAL=A(IADD1+IB)+A(IADD1+II)
                 T1IMAG=B(IADD1+IB)+B(IADD1+II)
                 T2REAL=A(IADD1+IB)-A(IADD1+II)
                 T2IMAG=B(IADD1+IB)-B(IADD1+II)
                 T3REAL=A(IADD1+IC)+A(IADD1+IH)
                 T3IMAG=B(IADD1+IC)+B(IADD1+IH)
                 T4REAL=A(IADD1+IC)-A(IADD1+IH)
                 T4IMAG=B(IADD1+IC)-B(IADD1+IH)
                 T5REAL=A(IADD1+ID)+A(IADD1+IG)
                 T5IMAG=B(IADD1+ID)+B(IADD1+IG)
                 T6REAL=A(IADD1+ID)-A(IADD1+IG)
                 T6IMAG=B(IADD1+ID)-B(IADD1+IG)
                 T7REAL=A(IADD1+IE)+A(IADD1+IF)
                 T7IMAG=B(IADD1+IE)+B(IADD1+IF)
                 T8REAL=A(IADD1+IE)-A(IADD1+IF)
                 T8IMAG=B(IADD1+IE)-B(IADD1+IF)
                 B1REAL=A(IADD1)+T5REAL
                 B1IMAG=B(IADD1)+T5IMAG
                 B2REAL=T1REAL+T3REAL+T7REAL
                 B2IMAG=T1IMAG+T3IMAG+T7IMAG
                 B3REAL=B1REAL-PO5*B2REAL
                 B3IMAG=B1IMAG-PO5*B2IMAG
                 B4REAL=(T2REAL-T4REAL+T8REAL)*SIN60
                 B4IMAG=(T2IMAG-T4IMAG+T8IMAG)*SIN60
                 B5REAL=-PO5*T5REAL+A(IADD1)
                 B5IMAG=-PO5*T5IMAG+B(IADD1)
                 B6REAL=T6REAL*SIN60
                 B6IMAG=T6IMAG*SIN60
                 B7REAL=T1REAL*COS40+T3REAL*COS80+B5REAL-T7REAL*COS20
                 B7IMAG=T1IMAG*COS40+T3IMAG*COS80+B5IMAG-T7IMAG*COS20
                 B8REAL=T2REAL*SIN40+T4REAL*SIN80+B6REAL+T8REAL*SIN20
                 B8IMAG=T2IMAG*SIN40+T4IMAG*SIN80+B6IMAG+T8IMAG*SIN20
                 B9REAL=T1REAL*COS80-T3REAL*COS20+B5REAL+T7REAL*COS40
                 B9IMAG=T1IMAG*COS80-T3IMAG*COS20+B5IMAG+T7IMAG*COS40
                 B10REAL=T2REAL*SIN80+T4REAL*SIN20-B6REAL-T8REAL*SIN40
                 B10IMAG=T2IMAG*SIN80+T4IMAG*SIN20-B6IMAG-T8IMAG*SIN40
                 B11REAL=-T1REAL*COS20+T3REAL*COS40+B5REAL+T7REAL*COS80
                 B11IMAG=-T1IMAG*COS20+T3IMAG*COS40+B5IMAG+T7IMAG*COS80
                 B12REAL=T2REAL*SIN20-T4REAL*SIN40+B6REAL-T8REAL*SIN80
                 B12IMAG=T2IMAG*SIN20-T4IMAG*SIN40+B6IMAG-T8IMAG*SIN80
                 X1REAL=B7REAL-B8IMAG
                 X1IMAG=B7IMAG+B8REAL
                 X2REAL=B9REAL-B10IMAG
                 X2IMAG=B9IMAG+B10REAL
                 X3REAL=B3REAL-B4IMAG
                 X3IMAG=B3IMAG+B4REAL
                 X4REAL=B11REAL-B12IMAG
                 X4IMAG=B11IMAG+B12REAL
                 X5REAL=B11REAL+B12IMAG
                 X5IMAG=B11IMAG-B12REAL
                 X6REAL=B3REAL+B4IMAG
                 X6IMAG=B3IMAG-B4REAL
                 X7REAL=B9REAL+B10IMAG
                 X7IMAG=B9IMAG-B10REAL
                 X8REAL=B7REAL+B8IMAG
                 X8IMAG=B7IMAG-B8REAL 
                 C(IADD2)   =B1REAL+B2REAL
                 D(IADD2)   =B1IMAG+B2IMAG
                 C(IADD2+JB)=TRIGS(KK1)*X1REAL-TRIGI(KK1)*X1IMAG
                 D(IADD2+JB)=TRIGI(KK1)*X1REAL+TRIGS(KK1)*X1IMAG
                 C(IADD2+JC)=TRIGS(KK2)*X2REAL-TRIGI(KK2)*X2IMAG
                 D(IADD2+JC)=TRIGI(KK2)*X2REAL+TRIGS(KK2)*X2IMAG
                 C(IADD2+JD)=TRIGS(KK3)*X3REAL-TRIGI(KK3)*X3IMAG
                 D(IADD2+JD)=TRIGI(KK3)*X3REAL+TRIGS(KK3)*X3IMAG
                 C(IADD2+JE)=TRIGS(KK4)*X4REAL-TRIGI(KK4)*X4IMAG
                 D(IADD2+JE)=TRIGI(KK4)*X4REAL+TRIGS(KK4)*X4IMAG
                 C(IADD2+JF)=TRIGS(KK5)*X5REAL-TRIGI(KK5)*X5IMAG
                 D(IADD2+JF)=TRIGI(KK5)*X5REAL+TRIGS(KK5)*X5IMAG
                 C(IADD2+JG)=TRIGS(KK6)*X6REAL-TRIGI(KK6)*X6IMAG
                 D(IADD2+JG)=TRIGI(KK6)*X6REAL+TRIGS(KK6)*X6IMAG
                 C(IADD2+JH)=TRIGS(KK7)*X7REAL-TRIGI(KK7)*X7IMAG
                 D(IADD2+JH)=TRIGI(KK7)*X7REAL+TRIGS(KK7)*X7IMAG
                 C(IADD2+JI)=TRIGS(KK8)*X8REAL-TRIGI(KK8)*X8IMAG
                 D(IADD2+JI)=TRIGI(KK8)*X8REAL+TRIGS(KK8)*X8IMAG
               ENDDO
               I=I+ISKIP
               J=J+ISKIP1
             ENDDO
             J=J+JUMP
           ENDDO
         ELSE
           DO L=1,LA
             IADD2=-ISKIP3+J
*vdir nodep             
CDIR$ IVDEP
CDEC$ IVDEP
             DO NU=0,NUMBER-1
               IADD1=NU+I
               IADD2=IADD2+ISKIP3
               T1REAL=A(IADD1+IB)+A(IADD1+II)
               T1IMAG=B(IADD1+IB)+B(IADD1+II)
               T2REAL=A(IADD1+IB)-A(IADD1+II)
               T2IMAG=B(IADD1+IB)-B(IADD1+II)
               T3REAL=A(IADD1+IC)+A(IADD1+IH)
               T3IMAG=B(IADD1+IC)+B(IADD1+IH)
               T4REAL=A(IADD1+IC)-A(IADD1+IH)
               T4IMAG=B(IADD1+IC)-B(IADD1+IH)
               T5REAL=A(IADD1+ID)+A(IADD1+IG)
               T5IMAG=B(IADD1+ID)+B(IADD1+IG)
               T6REAL=A(IADD1+ID)-A(IADD1+IG)
               T6IMAG=B(IADD1+ID)-B(IADD1+IG)
               T7REAL=A(IADD1+IE)+A(IADD1+IF)
               T7IMAG=B(IADD1+IE)+B(IADD1+IF)
               T8REAL=A(IADD1+IE)-A(IADD1+IF)
               T8IMAG=B(IADD1+IE)-B(IADD1+IF)
               B1REAL=A(IADD1)+T5REAL
               B1IMAG=B(IADD1)+T5IMAG
               B2REAL=T1REAL+T3REAL+T7REAL
               B2IMAG=T1IMAG+T3IMAG+T7IMAG
               B3REAL=B1REAL-PO5*B2REAL
               B3IMAG=B1IMAG-PO5*B2IMAG
               B4REAL=(T2REAL-T4REAL+T8REAL)*SIN60
               B4IMAG=(T2IMAG-T4IMAG+T8IMAG)*SIN60
               B5REAL=-PO5*T5REAL+A(IADD1)
               B5IMAG=-PO5*T5IMAG+B(IADD1)
               B6REAL=T6REAL*SIN60
               B6IMAG=T6IMAG*SIN60
               B7REAL=T1REAL*COS40+T3REAL*COS80+B5REAL-T7REAL*COS20
               B7IMAG=T1IMAG*COS40+T3IMAG*COS80+B5IMAG-T7IMAG*COS20
               B8REAL=T2REAL*SIN40+T4REAL*SIN80+B6REAL+T8REAL*SIN20
               B8IMAG=T2IMAG*SIN40+T4IMAG*SIN80+B6IMAG+T8IMAG*SIN20
               B9REAL=T1REAL*COS80-T3REAL*COS20+B5REAL+T7REAL*COS40
               B9IMAG=T1IMAG*COS80-T3IMAG*COS20+B5IMAG+T7IMAG*COS40
               B10REAL=T2REAL*SIN80+T4REAL*SIN20-B6REAL-T8REAL*SIN40
               B10IMAG=T2IMAG*SIN80+T4IMAG*SIN20-B6IMAG-T8IMAG*SIN40
               B11REAL=-T1REAL*COS20+T3REAL*COS40+B5REAL+T7REAL*COS80
               B11IMAG=-T1IMAG*COS20+T3IMAG*COS40+B5IMAG+T7IMAG*COS80
               B12REAL=T2REAL*SIN20-T4REAL*SIN40+B6REAL-T8REAL*SIN80
               B12IMAG=T2IMAG*SIN20-T4IMAG*SIN40+B6IMAG-T8IMAG*SIN80
               C(IADD2)   =B1REAL+B2REAL
               D(IADD2)   =B1IMAG+B2IMAG
               C(IADD2+JB)=B7REAL+B8IMAG
               D(IADD2+JB)=B7IMAG-B8REAL
               C(IADD2+JC)=B9REAL+B10IMAG
               D(IADD2+JC)=B9IMAG-B10REAL
               C(IADD2+JD)=B3REAL+B4IMAG
               D(IADD2+JD)=B3IMAG-B4REAL
               C(IADD2+JE)=B11REAL+B12IMAG
               D(IADD2+JE)=B11IMAG-B12REAL
               C(IADD2+JF)=B11REAL-B12IMAG
               D(IADD2+JF)=B11IMAG+B12REAL
               C(IADD2+JG)=B3REAL-B4IMAG
               D(IADD2+JG)=B3IMAG+B4REAL
               C(IADD2+JH)=B9REAL-B10IMAG
               D(IADD2+JH)=B9IMAG+B10REAL
               C(IADD2+JI)=B7REAL-B8IMAG
               D(IADD2+JI)=B7IMAG+B8REAL 
             ENDDO
             I=I+ISKIP
             J=J+ISKIP1
           ENDDO
           J=J+JUMP
           DO K=LA,M-LA,LA
             K2=K+K
             KK1=K2+1
             KK2=KK1+K2
             KK3=KK2+K2
             KK4=KK3+K2
             KK5=KK4+K2
             KK6=KK5+K2
             KK7=KK6+K2
             KK8=KK7+K2
             DO L=1,LA
               IADD2=-ISKIP3+J
*vdir nodep             
CDIR$ IVDEP
CDEC$ IVDEP
               DO NU=0,NUMBER-1
                 IADD1=NU+I
                 IADD2=IADD2+ISKIP3
                 T1REAL=A(IADD1+IB)+A(IADD1+II)
                 T1IMAG=B(IADD1+IB)+B(IADD1+II)
                 T2REAL=A(IADD1+IB)-A(IADD1+II)
                 T2IMAG=B(IADD1+IB)-B(IADD1+II)
                 T3REAL=A(IADD1+IC)+A(IADD1+IH)
                 T3IMAG=B(IADD1+IC)+B(IADD1+IH)
                 T4REAL=A(IADD1+IC)-A(IADD1+IH)
                 T4IMAG=B(IADD1+IC)-B(IADD1+IH)
                 T5REAL=A(IADD1+ID)+A(IADD1+IG)
                 T5IMAG=B(IADD1+ID)+B(IADD1+IG)
                 T6REAL=A(IADD1+ID)-A(IADD1+IG)
                 T6IMAG=B(IADD1+ID)-B(IADD1+IG)
                 T7REAL=A(IADD1+IE)+A(IADD1+IF)
                 T7IMAG=B(IADD1+IE)+B(IADD1+IF)
                 T8REAL=A(IADD1+IE)-A(IADD1+IF)
                 T8IMAG=B(IADD1+IE)-B(IADD1+IF)
                 B1REAL=A(IADD1)+T5REAL
                 B1IMAG=B(IADD1)+T5IMAG
                 B2REAL=T1REAL+T3REAL+T7REAL
                 B2IMAG=T1IMAG+T3IMAG+T7IMAG
                 B3REAL=B1REAL-PO5*B2REAL
                 B3IMAG=B1IMAG-PO5*B2IMAG
                 B4REAL=(T2REAL-T4REAL+T8REAL)*SIN60
                 B4IMAG=(T2IMAG-T4IMAG+T8IMAG)*SIN60
                 B5REAL=-PO5*T5REAL+A(IADD1)
                 B5IMAG=-PO5*T5IMAG+B(IADD1)
                 B6REAL=T6REAL*SIN60
                 B6IMAG=T6IMAG*SIN60
                 B7REAL=T1REAL*COS40+T3REAL*COS80+B5REAL-T7REAL*COS20
                 B7IMAG=T1IMAG*COS40+T3IMAG*COS80+B5IMAG-T7IMAG*COS20
                 B8REAL=T2REAL*SIN40+T4REAL*SIN80+B6REAL+T8REAL*SIN20
                 B8IMAG=T2IMAG*SIN40+T4IMAG*SIN80+B6IMAG+T8IMAG*SIN20
                 B9REAL=T1REAL*COS80-T3REAL*COS20+B5REAL+T7REAL*COS40
                 B9IMAG=T1IMAG*COS80-T3IMAG*COS20+B5IMAG+T7IMAG*COS40
                 B10REAL=T2REAL*SIN80+T4REAL*SIN20-B6REAL-T8REAL*SIN40
                 B10IMAG=T2IMAG*SIN80+T4IMAG*SIN20-B6IMAG-T8IMAG*SIN40
                 B11REAL=-T1REAL*COS20+T3REAL*COS40+B5REAL+T7REAL*COS80
                 B11IMAG=-T1IMAG*COS20+T3IMAG*COS40+B5IMAG+T7IMAG*COS80
                 B12REAL=T2REAL*SIN20-T4REAL*SIN40+B6REAL-T8REAL*SIN80
                 B12IMAG=T2IMAG*SIN20-T4IMAG*SIN40+B6IMAG-T8IMAG*SIN80
                 X1REAL=B7REAL+B8IMAG
                 X1IMAG=B7IMAG-B8REAL
                 X2REAL=B9REAL+B10IMAG
                 X2IMAG=B9IMAG-B10REAL
                 X3REAL=B3REAL+B4IMAG
                 X3IMAG=B3IMAG-B4REAL
                 X4REAL=B11REAL+B12IMAG
                 X4IMAG=B11IMAG-B12REAL
                 X5REAL=B11REAL-B12IMAG
                 X5IMAG=B11IMAG+B12REAL
                 X6REAL=B3REAL-B4IMAG
                 X6IMAG=B3IMAG+B4REAL
                 X7REAL=B9REAL-B10IMAG
                 X7IMAG=B9IMAG+B10REAL
                 X8REAL=B7REAL-B8IMAG
                 X8IMAG=B7IMAG+B8REAL 
                 C(IADD2)   = B1REAL+B2REAL
                 D(IADD2)   = B1IMAG+B2IMAG
                 C(IADD2+JB)= TRIGS(KK1)*X1REAL+TRIGI(KK1)*X1IMAG
                 D(IADD2+JB)=-TRIGI(KK1)*X1REAL+TRIGS(KK1)*X1IMAG
                 C(IADD2+JC)= TRIGS(KK2)*X2REAL+TRIGI(KK2)*X2IMAG
                 D(IADD2+JC)=-TRIGI(KK2)*X2REAL+TRIGS(KK2)*X2IMAG
                 C(IADD2+JD)= TRIGS(KK3)*X3REAL+TRIGI(KK3)*X3IMAG
                 D(IADD2+JD)=-TRIGI(KK3)*X3REAL+TRIGS(KK3)*X3IMAG
                 C(IADD2+JE)= TRIGS(KK4)*X4REAL+TRIGI(KK4)*X4IMAG
                 D(IADD2+JE)=-TRIGI(KK4)*X4REAL+TRIGS(KK4)*X4IMAG
                 C(IADD2+JF)= TRIGS(KK5)*X5REAL+TRIGI(KK5)*X5IMAG
                 D(IADD2+JF)=-TRIGI(KK5)*X5REAL+TRIGS(KK5)*X5IMAG
                 C(IADD2+JG)= TRIGS(KK6)*X6REAL+TRIGI(KK6)*X6IMAG
                 D(IADD2+JG)=-TRIGI(KK6)*X6REAL+TRIGS(KK6)*X6IMAG
                 C(IADD2+JH)= TRIGS(KK7)*X7REAL+TRIGI(KK7)*X7IMAG
                 D(IADD2+JH)=-TRIGI(KK7)*X7REAL+TRIGS(KK7)*X7IMAG
                 C(IADD2+JI)= TRIGS(KK8)*X8REAL+TRIGI(KK8)*X8IMAG
                 D(IADD2+JI)=-TRIGI(KK8)*X8REAL+TRIGS(KK8)*X8IMAG
               ENDDO
               I=I+ISKIP
               J=J+ISKIP1
             ENDDO
             J=J+JUMP
           ENDDO
         ENDIF
       ELSEIF (IFAC .EQ. 8) THEN
C
C        FACTOR 8
C
         IC=IB+MSKIP
         ID=IC+MSKIP
         IE=ID+MSKIP
         IF=IE+MSKIP
         IG=IF+MSKIP
         IH=IG+MSKIP
         JC=JB+LASKIP
         JD=JC+LASKIP
         JE=JD+LASKIP
         JF=JE+LASKIP
         JG=JF+LASKIP
         JH=JG+LASKIP
         IF(ISIGN.LT.0) THEN
           DO L=1,LA
             IADD2=-ISKIP3+J
*vdir nodep             
CDIR$ IVDEP
CDEC$ IVDEP
             DO NU=0,NUMBER-1
               IADD2=IADD2+ISKIP3
               IADD1=NU+I
               T1REAL=A(IADD1)+A(IADD1+IE)
               T1IMAG=B(IADD1)+B(IADD1+IE)
               T2REAL=A(IADD1)-A(IADD1+IE)
               T2IMAG=B(IADD1)-B(IADD1+IE)
               T3REAL=A(IADD1+IB)+A(IADD1+IF)
               T3IMAG=B(IADD1+IB)+B(IADD1+IF)
               T4REAL=A(IADD1+IB)-A(IADD1+IF)
               T4IMAG=B(IADD1+IB)-B(IADD1+IF)
               T5REAL=A(IADD1+IC)+A(IADD1+IG)
               T5IMAG=B(IADD1+IC)+B(IADD1+IG)
               T6REAL=A(IADD1+IC)-A(IADD1+IG)
               T6IMAG=B(IADD1+IC)-B(IADD1+IG)
               T7REAL=A(IADD1+ID)+A(IADD1+IH)
               T7IMAG=B(IADD1+ID)+B(IADD1+IH)
               T8REAL=A(IADD1+ID)-A(IADD1+IH)
               T8IMAG=B(IADD1+ID)-B(IADD1+IH)
               B1REAL=T1REAL+T5REAL
               B1IMAG=T1IMAG+T5IMAG
               B2REAL=T1REAL-T5REAL
               B2IMAG=T1IMAG-T5IMAG
               B3REAL=T3REAL+T7REAL
               B3IMAG=T3IMAG+T7IMAG
               B4REAL=T7IMAG-T3IMAG
               B4IMAG=T3REAL-T7REAL
               B5REAL=T2REAL-T6IMAG
               B5IMAG=T2IMAG+T6REAL
               B6REAL=T2REAL+T6IMAG
               B6IMAG=T2IMAG-T6REAL
               T1    =SQ2*(T4REAL+T4IMAG)
               T2    =SQ2*(T4REAL-T4IMAG)
               T3    =SQ2*(T8REAL+T8IMAG)
               T4    =SQ2*(T8REAL-T8IMAG)
               B7REAL=T2-T3
               B7IMAG=T1+T4
               B8REAL=T1-T4
               B8IMAG=-(T2+T3)
               C(IADD2)   =B1REAL+B3REAL
               D(IADD2)   =B1IMAG+B3IMAG
               C(IADD2+JB)=B5REAL+B7REAL
               D(IADD2+JB)=B5IMAG+B7IMAG
               C(IADD2+JC)=B2REAL+B4REAL
               D(IADD2+JC)=B2IMAG+B4IMAG
               C(IADD2+JD)=B6REAL-B8REAL
               D(IADD2+JD)=B6IMAG-B8IMAG
               C(IADD2+JE)=B1REAL-B3REAL
               D(IADD2+JE)=B1IMAG-B3IMAG
               C(IADD2+JF)=B5REAL-B7REAL
               D(IADD2+JF)=B5IMAG-B7IMAG
               C(IADD2+JG)=B2REAL-B4REAL
               D(IADD2+JG)=B2IMAG-B4IMAG
               C(IADD2+JH)=B6REAL+B8REAL
               D(IADD2+JH)=B6IMAG+B8IMAG
             ENDDO
             I=I+ISKIP
             J=J+ISKIP1
           ENDDO
           J=J+JUMP
           DO K=LA,M-LA,LA
             K2=K+K
             KK1=K2+1
             KK2=KK1+K2
             KK3=KK2+K2
             KK4=KK3+K2
             KK5=KK4+K2
             KK6=KK5+K2
             KK7=KK6+K2
             DO L=1,LA
               IADD2=-ISKIP3+J
*vdir nodep             
CDIR$ IVDEP
CDEC$ IVDEP
               DO NU=0,NUMBER-1
                 IADD2=IADD2+ISKIP3
                 IADD1=NU+I
                 T1REAL=A(IADD1)+A(IADD1+IE)
                 T1IMAG=B(IADD1)+B(IADD1+IE)
                 T2REAL=A(IADD1)-A(IADD1+IE)
                 T2IMAG=B(IADD1)-B(IADD1+IE)
                 T3REAL=A(IADD1+IB)+A(IADD1+IF)
                 T3IMAG=B(IADD1+IB)+B(IADD1+IF)
                 T4REAL=A(IADD1+IB)-A(IADD1+IF)
                 T4IMAG=B(IADD1+IB)-B(IADD1+IF)
                 T5REAL=A(IADD1+IC)+A(IADD1+IG)
                 T5IMAG=B(IADD1+IC)+B(IADD1+IG)
                 T6REAL=A(IADD1+IC)-A(IADD1+IG)
                 T6IMAG=B(IADD1+IC)-B(IADD1+IG)
                 T7REAL=A(IADD1+ID)+A(IADD1+IH)
                 T7IMAG=B(IADD1+ID)+B(IADD1+IH)
                 T8REAL=A(IADD1+ID)-A(IADD1+IH)
                 T8IMAG=B(IADD1+ID)-B(IADD1+IH)
                 B1REAL=T1REAL+T5REAL
                 B1IMAG=T1IMAG+T5IMAG
                 B2REAL=T1REAL-T5REAL
                 B2IMAG=T1IMAG-T5IMAG
                 B3REAL=T3REAL+T7REAL
                 B3IMAG=T3IMAG+T7IMAG
                 B4REAL=T7IMAG-T3IMAG
                 B4IMAG=T3REAL-T7REAL
                 B5REAL=T2REAL-T6IMAG
                 B5IMAG=T2IMAG+T6REAL
                 B6REAL=T2REAL+T6IMAG
                 B6IMAG=T2IMAG-T6REAL
                 T1    =SQ2*(T4REAL+T4IMAG)
                 T2    =SQ2*(T4REAL-T4IMAG)
                 T3    =SQ2*(T8REAL+T8IMAG)
                 T4    =SQ2*(T8REAL-T8IMAG)
                 B7REAL=T2-T3
                 B7IMAG=T1+T4
                 B8REAL=T1-T4
                 B8IMAG=-(T2+T3)
                 X1REAL=B5REAL+B7REAL
                 X1IMAG=B5IMAG+B7IMAG
                 X2REAL=B2REAL+B4REAL
                 X2IMAG=B2IMAG+B4IMAG
                 X3REAL=B6REAL-B8REAL
                 X3IMAG=B6IMAG-B8IMAG
                 X4REAL=B1REAL-B3REAL
                 X4IMAG=B1IMAG-B3IMAG
                 X5REAL=B5REAL-B7REAL
                 X5IMAG=B5IMAG-B7IMAG
                 X6REAL=B2REAL-B4REAL
                 X6IMAG=B2IMAG-B4IMAG
                 X7REAL=B6REAL+B8REAL
                 X7IMAG=B6IMAG+B8IMAG
                 C(IADD2)   =B1REAL+B3REAL
                 D(IADD2)   =B1IMAG+B3IMAG
                 C(IADD2+JB)=TRIGS(KK1)*X1REAL-TRIGI(KK1)*X1IMAG
                 D(IADD2+JB)=TRIGI(KK1)*X1REAL+TRIGS(KK1)*X1IMAG
                 C(IADD2+JC)=TRIGS(KK2)*X2REAL-TRIGI(KK2)*X2IMAG
                 D(IADD2+JC)=TRIGI(KK2)*X2REAL+TRIGS(KK2)*X2IMAG
                 C(IADD2+JD)=TRIGS(KK3)*X3REAL-TRIGI(KK3)*X3IMAG
                 D(IADD2+JD)=TRIGI(KK3)*X3REAL+TRIGS(KK3)*X3IMAG
                 C(IADD2+JE)=TRIGS(KK4)*X4REAL-TRIGI(KK4)*X4IMAG
                 D(IADD2+JE)=TRIGI(KK4)*X4REAL+TRIGS(KK4)*X4IMAG
                 C(IADD2+JF)=TRIGS(KK5)*X5REAL-TRIGI(KK5)*X5IMAG
                 D(IADD2+JF)=TRIGI(KK5)*X5REAL+TRIGS(KK5)*X5IMAG
                 C(IADD2+JG)=TRIGS(KK6)*X6REAL-TRIGI(KK6)*X6IMAG
                 D(IADD2+JG)=TRIGI(KK6)*X6REAL+TRIGS(KK6)*X6IMAG
                 C(IADD2+JH)=TRIGS(KK7)*X7REAL-TRIGI(KK7)*X7IMAG
                 D(IADD2+JH)=TRIGI(KK7)*X7REAL+TRIGS(KK7)*X7IMAG
               ENDDO
               I=I+ISKIP
               J=J+ISKIP1
             ENDDO
             J=J+JUMP
           ENDDO
         ELSE
           DO L=1,LA
             IADD2=-ISKIP3+J
*vdir nodep             
CDIR$ IVDEP
CDEC$ IVDEP
             DO NU=0,NUMBER-1
               IADD2=IADD2+ISKIP3
               IADD1=NU+I
               T1REAL=A(IADD1)+A(IADD1+IE)
               T1IMAG=B(IADD1)+B(IADD1+IE)
               T2REAL=A(IADD1)-A(IADD1+IE)
               T2IMAG=B(IADD1)-B(IADD1+IE)
               T3REAL=A(IADD1+IB)+A(IADD1+IF)
               T3IMAG=B(IADD1+IB)+B(IADD1+IF)
               T4REAL=A(IADD1+IB)-A(IADD1+IF)
               T4IMAG=B(IADD1+IB)-B(IADD1+IF)
               T5REAL=A(IADD1+IC)+A(IADD1+IG)
               T5IMAG=B(IADD1+IC)+B(IADD1+IG)
               T6REAL=A(IADD1+IC)-A(IADD1+IG)
               T6IMAG=B(IADD1+IC)-B(IADD1+IG)
               T7REAL=A(IADD1+ID)+A(IADD1+IH)
               T7IMAG=B(IADD1+ID)+B(IADD1+IH)
               T8REAL=A(IADD1+ID)-A(IADD1+IH)
               T8IMAG=B(IADD1+ID)-B(IADD1+IH)
               B1REAL=T1REAL+T5REAL
               B1IMAG=T1IMAG+T5IMAG
               B2REAL=T1REAL-T5REAL
               B2IMAG=T1IMAG-T5IMAG
               B3REAL=T3REAL+T7REAL
               B3IMAG=T3IMAG+T7IMAG
               B4REAL=T7IMAG-T3IMAG
               B4IMAG=T3REAL-T7REAL
               B5REAL=T2REAL-T6IMAG
               B5IMAG=T2IMAG+T6REAL
               B6REAL=T2REAL+T6IMAG
               B6IMAG=T2IMAG-T6REAL
               T1    =SQ2*(T4REAL+T4IMAG)
               T2    =SQ2*(T4REAL-T4IMAG)
               T3    =SQ2*(T8REAL+T8IMAG)
               T4    =SQ2*(T8REAL-T8IMAG)
               B7REAL=T2-T3
               B7IMAG=T1+T4
               B8REAL=T1-T4
               B8IMAG=-(T2+T3)
               C(IADD2)   =B1REAL+B3REAL
               D(IADD2)   =B1IMAG+B3IMAG
               C(IADD2+JB)=B6REAL+B8REAL
               D(IADD2+JB)=B6IMAG+B8IMAG
               C(IADD2+JC)=B2REAL-B4REAL
               D(IADD2+JC)=B2IMAG-B4IMAG
               C(IADD2+JD)=B5REAL-B7REAL
               D(IADD2+JD)=B5IMAG-B7IMAG
               C(IADD2+JE)=B1REAL-B3REAL
               D(IADD2+JE)=B1IMAG-B3IMAG
               C(IADD2+JF)=B6REAL-B8REAL
               D(IADD2+JF)=B6IMAG-B8IMAG
               C(IADD2+JG)=B2REAL+B4REAL
               D(IADD2+JG)=B2IMAG+B4IMAG
               C(IADD2+JH)=B5REAL+B7REAL
               D(IADD2+JH)=B5IMAG+B7IMAG
             ENDDO
             I=I+ISKIP
             J=J+ISKIP1
           ENDDO
           J=J+JUMP
           DO K=LA,M-LA,LA
             K2=K+K
             KK1=K2+1
             KK2=KK1+K2
             KK3=KK2+K2
             KK4=KK3+K2
             KK5=KK4+K2
             KK6=KK5+K2
             KK7=KK6+K2
             DO L=1,LA
               IADD2=-ISKIP3+J
*vdir nodep             
CDIR$ IVDEP
CDEC$ IVDEP
               DO NU=0,NUMBER-1
                 IADD2=IADD2+ISKIP3
                 IADD1=NU+I
                 T1REAL=A(IADD1)+A(IADD1+IE)
                 T1IMAG=B(IADD1)+B(IADD1+IE)
                 T2REAL=A(IADD1)-A(IADD1+IE)
                 T2IMAG=B(IADD1)-B(IADD1+IE)
                 T3REAL=A(IADD1+IB)+A(IADD1+IF)
                 T3IMAG=B(IADD1+IB)+B(IADD1+IF)
                 T4REAL=A(IADD1+IB)-A(IADD1+IF)
                 T4IMAG=B(IADD1+IB)-B(IADD1+IF)
                 T5REAL=A(IADD1+IC)+A(IADD1+IG)
                 T5IMAG=B(IADD1+IC)+B(IADD1+IG)
                 T6REAL=A(IADD1+IC)-A(IADD1+IG)
                 T6IMAG=B(IADD1+IC)-B(IADD1+IG)
                 T7REAL=A(IADD1+ID)+A(IADD1+IH)
                 T7IMAG=B(IADD1+ID)+B(IADD1+IH)
                 T8REAL=A(IADD1+ID)-A(IADD1+IH)
                 T8IMAG=B(IADD1+ID)-B(IADD1+IH)
                 B1REAL=T1REAL+T5REAL
                 B1IMAG=T1IMAG+T5IMAG
                 B2REAL=T1REAL-T5REAL
                 B2IMAG=T1IMAG-T5IMAG
                 B3REAL=T3REAL+T7REAL
                 B3IMAG=T3IMAG+T7IMAG
                 B4REAL=T7IMAG-T3IMAG
                 B4IMAG=T3REAL-T7REAL
                 B5REAL=T2REAL-T6IMAG
                 B5IMAG=T2IMAG+T6REAL
                 B6REAL=T2REAL+T6IMAG
                 B6IMAG=T2IMAG-T6REAL
                 T1    =SQ2*(T4REAL+T4IMAG)
                 T2    =SQ2*(T4REAL-T4IMAG)
                 T3    =SQ2*(T8REAL+T8IMAG)
                 T4    =SQ2*(T8REAL-T8IMAG)
                 B7REAL=T2-T3
                 B7IMAG=T1+T4
                 B8REAL=T1-T4
                 B8IMAG=-(T2+T3)
                 X1REAL=B6REAL+B8REAL
                 X1IMAG=B6IMAG+B8IMAG
                 X2REAL=B2REAL-B4REAL
                 X2IMAG=B2IMAG-B4IMAG
                 X3REAL=B5REAL-B7REAL
                 X3IMAG=B5IMAG-B7IMAG
                 X4REAL=B1REAL-B3REAL
                 X4IMAG=B1IMAG-B3IMAG
                 X5REAL=B6REAL-B8REAL
                 X5IMAG=B6IMAG-B8IMAG
                 X6REAL=B2REAL+B4REAL
                 X6IMAG=B2IMAG+B4IMAG
                 X7REAL=B5REAL+B7REAL
                 X7IMAG=B5IMAG+B7IMAG
                 C(IADD2)   =B1REAL+B3REAL
                 D(IADD2)   =B1IMAG+B3IMAG
                 C(IADD2+JB)= TRIGS(KK1)*X1REAL+TRIGI(KK1)*X1IMAG
                 D(IADD2+JB)=-TRIGI(KK1)*X1REAL+TRIGS(KK1)*X1IMAG
                 C(IADD2+JC)= TRIGS(KK2)*X2REAL+TRIGI(KK2)*X2IMAG
                 D(IADD2+JC)=-TRIGI(KK2)*X2REAL+TRIGS(KK2)*X2IMAG
                 C(IADD2+JD)= TRIGS(KK3)*X3REAL+TRIGI(KK3)*X3IMAG
                 D(IADD2+JD)=-TRIGI(KK3)*X3REAL+TRIGS(KK3)*X3IMAG
                 C(IADD2+JE)= TRIGS(KK4)*X4REAL+TRIGI(KK4)*X4IMAG
                 D(IADD2+JE)=-TRIGI(KK4)*X4REAL+TRIGS(KK4)*X4IMAG
                 C(IADD2+JF)= TRIGS(KK5)*X5REAL+TRIGI(KK5)*X5IMAG
                 D(IADD2+JF)=-TRIGI(KK5)*X5REAL+TRIGS(KK5)*X5IMAG
                 C(IADD2+JG)= TRIGS(KK6)*X6REAL+TRIGI(KK6)*X6IMAG
                 D(IADD2+JG)=-TRIGI(KK6)*X6REAL+TRIGS(KK6)*X6IMAG
                 C(IADD2+JH)= TRIGS(KK7)*X7REAL+TRIGI(KK7)*X7IMAG
                 D(IADD2+JH)=-TRIGI(KK7)*X7REAL+TRIGS(KK7)*X7IMAG
               ENDDO
               I=I+ISKIP
               J=J+ISKIP1
             ENDDO
             J=J+JUMP
           ENDDO
         ENDIF
       ELSEIF (IFAC .EQ. 6) THEN
C
C        FACTOR 6
C
         IC=IB+MSKIP
         ID=IC+MSKIP
         IE=ID+MSKIP
         IF=IE+MSKIP
         JC=JB+LASKIP
         JD=JC+LASKIP
         JE=JD+LASKIP
         JF=JE+LASKIP
         IF(ISIGN.LT.0) THEN
           DO L=1,LA
             IADD2=-ISKIP3+J
*vdir nodep             
CDIR$ IVDEP
CDEC$ IVDEP
             DO NU=0,NUMBER-1
               IADD2=IADD2+ISKIP3
               IADD1=NU+I
               T1REAL=A(IADD1)+A(IADD1+ID)
               T1IMAG=B(IADD1)+B(IADD1+ID)
               T2REAL=A(IADD1)-A(IADD1+ID)
               T2IMAG=B(IADD1)-B(IADD1+ID)
               T3REAL=A(IADD1+IB)+A(IADD1+IE)
               T3IMAG=B(IADD1+IB)+B(IADD1+IE)
               T4REAL=A(IADD1+IB)-A(IADD1+IE)
               T4IMAG=B(IADD1+IB)-B(IADD1+IE)
               T5REAL=A(IADD1+IC)+A(IADD1+IF)
               T5IMAG=B(IADD1+IC)+B(IADD1+IF)
               T6REAL=A(IADD1+IC)-A(IADD1+IF)
               T6IMAG=B(IADD1+IC)-B(IADD1+IF)
               B1REAL=T3REAL+T5REAL
               B1IMAG=T3IMAG+T5IMAG
               B1REAL5=PO5*B1REAL
               B1IMAG5=PO5*B1IMAG
               B2REAL=(T3REAL-T5REAL)*SIN60
               B2IMAG=(T3IMAG-T5IMAG)*SIN60
               B3REAL=T4REAL-T6REAL
               B3IMAG=T4IMAG-T6IMAG
               B3REAL5=PO5*B3REAL
               B3IMAG5=PO5*B3IMAG
               B4REAL=(T4REAL+T6REAL)*SIN60
               B4IMAG=(T4IMAG+T6IMAG)*SIN60
               C(IADD2)   =T1REAL+B1REAL
               D(IADD2)   =T1IMAG+B1IMAG
               C(IADD2+JB)=T2REAL+B3REAL5-B4IMAG
               D(IADD2+JB)=T2IMAG+B3IMAG5+B4REAL
               C(IADD2+JC)=T1REAL-B1REAL5-B2IMAG
               D(IADD2+JC)=T1IMAG-B1IMAG5+B2REAL
               C(IADD2+JD)=T2REAL-B3REAL
               D(IADD2+JD)=T2IMAG-B3IMAG
               C(IADD2+JE)=T1REAL-B1REAL5+B2IMAG
               D(IADD2+JE)=T1IMAG-B1IMAG5-B2REAL
               C(IADD2+JF)=T2REAL+B3REAL5+B4IMAG
               D(IADD2+JF)=T2IMAG+B3IMAG5-B4REAL
             ENDDO
             I=I+ISKIP
             J=J+ISKIP1
           ENDDO
           J=J+JUMP
           DO K=LA,M-LA,LA
             K2=K+K
             KK=K2+1
             KK2=KK+K2
             KK4=KK2+K2
             KK6=KK4+K2
             KK8=KK6+K2
             DO L=1,LA
               IADD2=-ISKIP3+J
*vdir nodep             
CDIR$ IVDEP
               DO NU=0,NUMBER-1
                 IADD2=IADD2+ISKIP3
                 IADD1=NU+I
                 T1REAL=A(IADD1)+A(IADD1+ID)
                 T1IMAG=B(IADD1)+B(IADD1+ID)
                 T2REAL=A(IADD1)-A(IADD1+ID)
                 T2IMAG=B(IADD1)-B(IADD1+ID)
                 T3REAL=A(IADD1+IB)+A(IADD1+IE)
                 T3IMAG=B(IADD1+IB)+B(IADD1+IE)
                 T4REAL=A(IADD1+IB)-A(IADD1+IE)
                 T4IMAG=B(IADD1+IB)-B(IADD1+IE)
                 T5REAL=A(IADD1+IC)+A(IADD1+IF)
                 T5IMAG=B(IADD1+IC)+B(IADD1+IF)
                 T6REAL=A(IADD1+IC)-A(IADD1+IF)
                 T6IMAG=B(IADD1+IC)-B(IADD1+IF)
                 B1REAL=T3REAL+T5REAL
                 B1IMAG=T3IMAG+T5IMAG
                 B1REAL5=PO5*B1REAL
                 B1IMAG5=PO5*B1IMAG
                 B2REAL=(T3REAL-T5REAL)*SIN60
                 B2IMAG=(T3IMAG-T5IMAG)*SIN60
                 B3REAL=T4REAL-T6REAL
                 B3IMAG=T4IMAG-T6IMAG
                 B3REAL5=PO5*B3REAL
                 B3IMAG5=PO5*B3IMAG
                 B4REAL=(T4REAL+T6REAL)*SIN60
                 B4IMAG=(T4IMAG+T6IMAG)*SIN60
                 X2REAL=T2REAL+B3REAL5-B4IMAG
                 X2IMAG=T2IMAG+B3IMAG5+B4REAL
                 X3REAL=T1REAL-B1REAL5-B2IMAG
                 X3IMAG=T1IMAG-B1IMAG5+B2REAL
                 X4REAL=T2REAL-B3REAL
                 X4IMAG=T2IMAG-B3IMAG
                 X5REAL=T1REAL-B1REAL5+B2IMAG
                 X5IMAG=T1IMAG-B1IMAG5-B2REAL
                 X6REAL=T2REAL+B3REAL5+B4IMAG
                 X6IMAG=T2IMAG+B3IMAG5-B4REAL
                 C(IADD2)   =T1REAL+B1REAL
                 D(IADD2)   =T1IMAG+B1IMAG
                 C(IADD2+JB)=TRIGS(KK)*X2REAL-TRIGI(KK)*X2IMAG
                 D(IADD2+JB)=TRIGI(KK)*X2REAL+TRIGS(KK)*X2IMAG
                 C(IADD2+JC)=TRIGS(KK2)*X3REAL-TRIGI(KK2)*X3IMAG
                 D(IADD2+JC)=TRIGI(KK2)*X3REAL+TRIGS(KK2)*X3IMAG
                 C(IADD2+JD)=TRIGS(KK4)*X4REAL-TRIGI(KK4)*X4IMAG
                 D(IADD2+JD)=TRIGI(KK4)*X4REAL+TRIGS(KK4)*X4IMAG
                 C(IADD2+JE)=TRIGS(KK6)*X5REAL-TRIGI(KK6)*X5IMAG
                 D(IADD2+JE)=TRIGI(KK6)*X5REAL+TRIGS(KK6)*X5IMAG
                 C(IADD2+JF)=TRIGS(KK8)*X6REAL-TRIGI(KK8)*X6IMAG
                 D(IADD2+JF)=TRIGI(KK8)*X6REAL+TRIGS(KK8)*X6IMAG
               ENDDO
               I=I+ISKIP
               J=J+ISKIP1
             ENDDO
             J=J+JUMP
           ENDDO
         ELSE
           DO L=1,LA
             IADD2=-ISKIP3+J
*vdir nodep             
CDIR$ IVDEP
CDEC$ IVDEP
             DO NU=0,NUMBER-1
               IADD2=IADD2+ISKIP3
               IADD1=NU+I
               T1REAL=A(IADD1)+A(IADD1+ID)
               T1IMAG=B(IADD1)+B(IADD1+ID)
               T2REAL=A(IADD1)-A(IADD1+ID)
               T2IMAG=B(IADD1)-B(IADD1+ID)
               T3REAL=A(IADD1+IB)+A(IADD1+IE)
               T3IMAG=B(IADD1+IB)+B(IADD1+IE)
               T4REAL=A(IADD1+IB)-A(IADD1+IE)
               T4IMAG=B(IADD1+IB)-B(IADD1+IE)
               T5REAL=A(IADD1+IC)+A(IADD1+IF)
               T5IMAG=B(IADD1+IC)+B(IADD1+IF)
               T6REAL=A(IADD1+IC)-A(IADD1+IF)
               T6IMAG=B(IADD1+IC)-B(IADD1+IF)
               B1REAL=T3REAL+T5REAL
               B1IMAG=T3IMAG+T5IMAG
               B1REAL5=PO5*B1REAL
               B1IMAG5=PO5*B1IMAG
               B2REAL=(T3REAL-T5REAL)*SIN60
               B2IMAG=(T3IMAG-T5IMAG)*SIN60
               B3REAL=T4REAL-T6REAL
               B3IMAG=T4IMAG-T6IMAG
               B3REAL5=PO5*B3REAL
               B3IMAG5=PO5*B3IMAG
               B4REAL=(T4REAL+T6REAL)*SIN60
               B4IMAG=(T4IMAG+T6IMAG)*SIN60
               C(IADD2)   =T1REAL+B1REAL
               D(IADD2)   =T1IMAG+B1IMAG
               C(IADD2+JB)=T2REAL+B3REAL5+B4IMAG
               D(IADD2+JB)=T2IMAG+B3IMAG5-B4REAL
               C(IADD2+JC)=T1REAL-B1REAL5+B2IMAG
               D(IADD2+JC)=T1IMAG-B1IMAG5-B2REAL
               C(IADD2+JD)=T2REAL-B3REAL
               D(IADD2+JD)=T2IMAG-B3IMAG
               C(IADD2+JE)=T1REAL-B1REAL5-B2IMAG
               D(IADD2+JE)=T1IMAG-B1IMAG5+B2REAL
               C(IADD2+JF)=T2REAL+B3REAL5-B4IMAG
               D(IADD2+JF)=T2IMAG+B3IMAG5+B4REAL
             ENDDO
             I=I+ISKIP
             J=J+ISKIP1
           ENDDO
           J=J+JUMP
           DO K=LA,M-LA,LA
             K2=K+K
             KK=K2+1
             KK2=KK+K2
             KK4=KK2+K2
             KK6=KK4+K2
             KK8=KK6+K2
             DO L=1,LA
               IADD2=-ISKIP3+J
*vdir nodep             
CDIR$ IVDEP
CDEC$ IVDEP
               DO NU=0,NUMBER-1
                 IADD2=IADD2+ISKIP3
                 IADD1=NU+I
                 T1REAL=A(IADD1)+A(IADD1+ID)
                 T1IMAG=B(IADD1)+B(IADD1+ID)
                 T2REAL=A(IADD1)-A(IADD1+ID)
                 T2IMAG=B(IADD1)-B(IADD1+ID)
                 T3REAL=A(IADD1+IB)+A(IADD1+IE)
                 T3IMAG=B(IADD1+IB)+B(IADD1+IE)
                 T4REAL=A(IADD1+IB)-A(IADD1+IE)
                 T4IMAG=B(IADD1+IB)-B(IADD1+IE)
                 T5REAL=A(IADD1+IC)+A(IADD1+IF)
                 T5IMAG=B(IADD1+IC)+B(IADD1+IF)
                 T6REAL=A(IADD1+IC)-A(IADD1+IF)
                 T6IMAG=B(IADD1+IC)-B(IADD1+IF)
                 B1REAL=T3REAL+T5REAL
                 B1IMAG=T3IMAG+T5IMAG
                 B1REAL5=PO5*B1REAL
                 B1IMAG5=PO5*B1IMAG
                 B2REAL=(T3REAL-T5REAL)*SIN60
                 B2IMAG=(T3IMAG-T5IMAG)*SIN60
                 B3REAL=T4REAL-T6REAL
                 B3IMAG=T4IMAG-T6IMAG
                 B3REAL5=PO5*B3REAL
                 B3IMAG5=PO5*B3IMAG
                 B4REAL=(T4REAL+T6REAL)*SIN60
                 B4IMAG=(T4IMAG+T6IMAG)*SIN60
                 X2REAL=T2REAL+B3REAL5+B4IMAG
                 X2IMAG=T2IMAG+B3IMAG5-B4REAL
                 X3REAL=T1REAL-B1REAL5+B2IMAG
                 X3IMAG=T1IMAG-B1IMAG5-B2REAL
                 X4REAL=T2REAL-B3REAL
                 X4IMAG=T2IMAG-B3IMAG
                 X5REAL=T1REAL-B1REAL5-B2IMAG
                 X5IMAG=T1IMAG-B1IMAG5+B2REAL
                 X6REAL=T2REAL+B3REAL5-B4IMAG
                 X6IMAG=T2IMAG+B3IMAG5+B4REAL
                 C(IADD2)   = T1REAL+B1REAL
                 D(IADD2)   = T1IMAG+B1IMAG
                 C(IADD2+JB)= TRIGS(KK)*X2REAL+TRIGI(KK)*X2IMAG
                 D(IADD2+JB)=-TRIGI(KK)*X2REAL+TRIGS(KK)*X2IMAG
                 C(IADD2+JC)= TRIGS(KK2)*X3REAL+TRIGI(KK2)*X3IMAG
                 D(IADD2+JC)=-TRIGI(KK2)*X3REAL+TRIGS(KK2)*X3IMAG
                 C(IADD2+JD)= TRIGS(KK4)*X4REAL+TRIGI(KK4)*X4IMAG
                 D(IADD2+JD)=-TRIGI(KK4)*X4REAL+TRIGS(KK4)*X4IMAG
                 C(IADD2+JE)= TRIGS(KK6)*X5REAL+TRIGI(KK6)*X5IMAG
                 D(IADD2+JE)=-TRIGI(KK6)*X5REAL+TRIGS(KK6)*X5IMAG
                 C(IADD2+JF)= TRIGS(KK8)*X6REAL+TRIGI(KK8)*X6IMAG
                 D(IADD2+JF)=-TRIGI(KK8)*X6REAL+TRIGS(KK8)*X6IMAG
               ENDDO
               I=I+ISKIP
               J=J+ISKIP1
             ENDDO
             J=J+JUMP
           ENDDO
         ENDIF
       ELSEIF (IFAC .EQ. 5 ) THEN
C
C        FACTOR 5
C
         IC=IB+MSKIP
         ID=IC+MSKIP
         IE=ID+MSKIP
         JC=JB+LASKIP
         JD=JC+LASKIP
         JE=JD+LASKIP
         IF(ISIGN.LT.0) THEN
           DO L=1,LA
             IADD1=-ISKIP3+J
*vdir nodep             
CDIR$ IVDEP
CDEC$ IVDEP
             DO NU=0,NUMBER-1
               IADD1=IADD1+ISKIP3
               IADD2=NU+I
               T1REAL=A(IADD2+IB)+A(IADD2+IE)
               T1IMAG=B(IADD2+IB)+B(IADD2+IE)
               T2REAL=A(IADD2+IC)+A(IADD2+ID)
               T2IMAG=B(IADD2+IC)+B(IADD2+ID)
               T3REAL=A(IADD2+IB)-A(IADD2+IE)
               T3IMAG=B(IADD2+IB)-B(IADD2+IE)
               T4REAL=A(IADD2+IC)-A(IADD2+ID)
               T4IMAG=B(IADD2+IC)-B(IADD2+ID)
               T5REAL=T1REAL+T2REAL
               T5IMAG=T1IMAG+T2IMAG
               T6REAL=SQ54*(T1REAL-T2REAL)
               T6IMAG=SQ54*(T1IMAG-T2IMAG)
               T7REAL=A(IADD2)-T5REAL*PO25
               T7IMAG=B(IADD2)-T5IMAG*PO25
               T8REAL=T7REAL+T6REAL
               T8IMAG=T7IMAG+T6IMAG
               T9REAL=T7REAL-T6REAL
               T9IMAG=T7IMAG-T6IMAG
               T10REAL=SIN72*T3REAL+SIN36*T4REAL
               T10IMAG=SIN72*T3IMAG+SIN36*T4IMAG
               T11REAL=SIN36*T3REAL-SIN72*T4REAL
               T11IMAG=SIN36*T3IMAG-SIN72*T4IMAG
               C(IADD1)=A(IADD2)+T5REAL
               D(IADD1)=B(IADD2)+T5IMAG
               C(IADD1+JB)=T8REAL-T10IMAG
               D(IADD1+JB)=T8IMAG+T10REAL
               C(IADD1+JC)=T9REAL-T11IMAG
               D(IADD1+JC)=T9IMAG+T11REAL
               C(IADD1+JD)=T9REAL+T11IMAG
               D(IADD1+JD)=T9IMAG-T11REAL
               C(IADD1+JE)=T8REAL+T10IMAG
               D(IADD1+JE)=T8IMAG-T10REAL
             ENDDO
             I=I+ISKIP
             J=J+ISKIP1
           ENDDO
           J=J+JUMP
           DO K=LA,M-LA,LA
             K2=K+K
             KK=K2+1
             KK2=KK+K2
             KK4=KK2+K2
             KK6=KK4+K2
             DO L=1,LA
               IADD1=-ISKIP3+J
*vdir nodep             
CDIR$ IVDEP
CDEC$ IVDEP
               DO NU=0,NUMBER-1
                 IADD1=IADD1+ISKIP3
                 IADD2=NU+I
                 T1REAL=A(IADD2+IB)+A(IADD2+IE)
                 T1IMAG=B(IADD2+IB)+B(IADD2+IE)
                 T2REAL=A(IADD2+IC)+A(IADD2+ID)
                 T2IMAG=B(IADD2+IC)+B(IADD2+ID)
                 T3REAL=A(IADD2+IB)-A(IADD2+IE)
                 T3IMAG=B(IADD2+IB)-B(IADD2+IE)
                 T4REAL=A(IADD2+IC)-A(IADD2+ID)
                 T4IMAG=B(IADD2+IC)-B(IADD2+ID)
                 T5REAL=T1REAL+T2REAL
                 T5IMAG=T1IMAG+T2IMAG
                 T6REAL=SQ54*(T1REAL-T2REAL)
                 T6IMAG=SQ54*(T1IMAG-T2IMAG)
                 T7REAL=A(IADD2)-T5REAL*PO25
                 T7IMAG=B(IADD2)-T5IMAG*PO25
                 T8REAL=T7REAL+T6REAL
                 T8IMAG=T7IMAG+T6IMAG
                 T9REAL=T7REAL-T6REAL
                 T9IMAG=T7IMAG-T6IMAG
                 T10REAL=SIN72*T3REAL+SIN36*T4REAL
                 T10IMAG=SIN72*T3IMAG+SIN36*T4IMAG
                 T11REAL=SIN36*T3REAL-SIN72*T4REAL
                 T11IMAG=SIN36*T3IMAG-SIN72*T4IMAG
                 X1REAL=T8REAL-T10IMAG
                 X1IMAG=T8IMAG+T10REAL
                 X2REAL=T9REAL-T11IMAG
                 X2IMAG=T9IMAG+T11REAL
                 X3REAL=T9REAL+T11IMAG
                 X3IMAG=T9IMAG-T11REAL
                 X4REAL=T8REAL+T10IMAG
                 X4IMAG=T8IMAG-T10REAL
                 C(IADD1)   =A(IADD2)+T5REAL
                 D(IADD1)   =B(IADD2)+T5IMAG
                 C(IADD1+JB)=TRIGS(KK) *X1REAL-TRIGI(KK)*X1IMAG
                 D(IADD1+JB)=TRIGI(KK)*X1REAL+TRIGS(KK) *X1IMAG
                 C(IADD1+JC)=TRIGS(KK2)*X2REAL-TRIGI(KK2)*X2IMAG
                 D(IADD1+JC)=TRIGI(KK2)*X2REAL+TRIGS(KK2)*X2IMAG
                 C(IADD1+JD)=TRIGS(KK4)*X3REAL-TRIGI(KK4)*X3IMAG
                 D(IADD1+JD)=TRIGI(KK4)*X3REAL+TRIGS(KK4)*X3IMAG
                 C(IADD1+JE)=TRIGS(KK6)*X4REAL-TRIGI(KK6)*X4IMAG
                 D(IADD1+JE)=TRIGI(KK6)*X4REAL+TRIGS(KK6)*X4IMAG
               ENDDO
               I=I+ISKIP
               J=J+ISKIP1
             ENDDO
             J=J+JUMP
           ENDDO
         ELSE
           DO L=1,LA
             IADD1=-ISKIP3+J
*vdir nodep             
CDIR$ IVDEP
CDEC$ IVDEP
             DO NU=0,NUMBER-1
               IADD1=IADD1+ISKIP3
               IADD2=NU+I
               T1REAL=A(IADD2+IB)+A(IADD2+IE)
               T1IMAG=B(IADD2+IB)+B(IADD2+IE)
               T2REAL=A(IADD2+IC)+A(IADD2+ID)
               T2IMAG=B(IADD2+IC)+B(IADD2+ID)
               T3REAL=A(IADD2+IB)-A(IADD2+IE)
               T3IMAG=B(IADD2+IB)-B(IADD2+IE)
               T4REAL=A(IADD2+IC)-A(IADD2+ID)
               T4IMAG=B(IADD2+IC)-B(IADD2+ID)
               T5REAL=T1REAL+T2REAL
               T5IMAG=T1IMAG+T2IMAG
               T6REAL=SQ54*(T1REAL-T2REAL)
               T6IMAG=SQ54*(T1IMAG-T2IMAG)
               T7REAL=A(IADD2)-T5REAL*PO25
               T7IMAG=B(IADD2)-T5IMAG*PO25
               T8REAL=T7REAL+T6REAL
               T8IMAG=T7IMAG+T6IMAG
               T9REAL=T7REAL-T6REAL
               T9IMAG=T7IMAG-T6IMAG
               T10REAL=SIN72*T3REAL+SIN36*T4REAL
               T10IMAG=SIN72*T3IMAG+SIN36*T4IMAG
               T11REAL=SIN36*T3REAL-SIN72*T4REAL
               T11IMAG=SIN36*T3IMAG-SIN72*T4IMAG
               C(IADD1)   =A(IADD2)+T5REAL
               D(IADD1)   =B(IADD2)+T5IMAG
               C(IADD1+JB)=T8REAL+T10IMAG
               D(IADD1+JB)=T8IMAG-T10REAL
               C(IADD1+JC)=T9REAL+T11IMAG
               D(IADD1+JC)=T9IMAG-T11REAL
               C(IADD1+JD)=T9REAL-T11IMAG
               D(IADD1+JD)=T9IMAG+T11REAL
               C(IADD1+JE)=T8REAL-T10IMAG
               D(IADD1+JE)=T8IMAG+T10REAL
             ENDDO
             I=I+ISKIP
             J=J+ISKIP1
           ENDDO
           J=J+JUMP
           DO K=LA,M-LA,LA
             K2=K+K
             KK=K2+1
             KK2=KK+K2
             KK4=KK2+K2
             KK6=KK4+K2
             DO L=1,LA
               IADD1=-ISKIP3+J
*vdir nodep             
CDIR$ IVDEP
CDEC$ IVDEP
               DO NU=0,NUMBER-1
                 IADD1=IADD1+ISKIP3
                 IADD2=NU+I
                 T1REAL=A(IADD2+IB)+A(IADD2+IE)
                 T1IMAG=B(IADD2+IB)+B(IADD2+IE)
                 T2REAL=A(IADD2+IC)+A(IADD2+ID)
                 T2IMAG=B(IADD2+IC)+B(IADD2+ID)
                 T3REAL=A(IADD2+IB)-A(IADD2+IE)
                 T3IMAG=B(IADD2+IB)-B(IADD2+IE)
                 T4REAL=A(IADD2+IC)-A(IADD2+ID)
                 T4IMAG=B(IADD2+IC)-B(IADD2+ID)
                 T5REAL=T1REAL+T2REAL
                 T5IMAG=T1IMAG+T2IMAG
                 T6REAL=SQ54*(T1REAL-T2REAL)
                 T6IMAG=SQ54*(T1IMAG-T2IMAG)
                 T7REAL=A(IADD2)-T5REAL*PO25
                 T7IMAG=B(IADD2)-T5IMAG*PO25
                 T8REAL=T7REAL+T6REAL
                 T8IMAG=T7IMAG+T6IMAG
                 T9REAL=T7REAL-T6REAL
                 T9IMAG=T7IMAG-T6IMAG
                 T10REAL=SIN72*T3REAL+SIN36*T4REAL
                 T10IMAG=SIN72*T3IMAG+SIN36*T4IMAG
                 T11REAL=SIN36*T3REAL-SIN72*T4REAL
                 T11IMAG=SIN36*T3IMAG-SIN72*T4IMAG
                 X1REAL=T8REAL+T10IMAG
                 X1IMAG=T8IMAG-T10REAL
                 X2REAL=T9REAL+T11IMAG
                 X2IMAG=T9IMAG-T11REAL
                 X3REAL=T9REAL-T11IMAG
                 X3IMAG=T9IMAG+T11REAL
                 X4REAL=T8REAL-T10IMAG
                 X4IMAG=T8IMAG+T10REAL
                 C(IADD1)   =A(IADD2)+T5REAL
                 D(IADD1)   =B(IADD2)+T5IMAG
                 C(IADD1+JB)= TRIGS(KK) *X1REAL+TRIGI(KK) *X1IMAG
                 D(IADD1+JB)=-TRIGI(KK) *X1REAL+TRIGS(KK) *X1IMAG
                 C(IADD1+JC)= TRIGS(KK2)*X2REAL+TRIGI(KK2)*X2IMAG
                 D(IADD1+JC)=-TRIGI(KK2)*X2REAL+TRIGS(KK2)*X2IMAG
                 C(IADD1+JD)= TRIGS(KK4)*X3REAL+TRIGI(KK4)*X3IMAG
                 D(IADD1+JD)=-TRIGI(KK4)*X3REAL+TRIGS(KK4)*X3IMAG
                 C(IADD1+JE)= TRIGS(KK6)*X4REAL+TRIGI(KK6)*X4IMAG
                 D(IADD1+JE)=-TRIGI(KK6)*X4REAL+TRIGS(KK6)*X4IMAG
               ENDDO
               I=I+ISKIP
               J=J+ISKIP1
             ENDDO
             J=J+JUMP
           ENDDO
         ENDIF
       ELSEIF (IFAC .EQ. 4) THEN
C
C        FACTOR 4 
C
         IC=IB+MSKIP
         ID=IC+MSKIP
         JC=JB+LASKIP
         JD=JC+LASKIP
         IF(ISIGN.LT.0) THEN
           DO L=1,LA
             IADD1=-ISKIP3+J
*vdir nodep             
CDIR$ IVDEP
CDEC$ IVDEP
             DO NU=0,NUMBER-1
               IADD1=IADD1+ISKIP3
               IADD2=NU+I
               T1REAL=A(IADD2)+A(IADD2+IC)
               T1IMAG=B(IADD2)+B(IADD2+IC)
               T2REAL=A(IADD2+IB)+A(IADD2+ID)
               T2IMAG=B(IADD2+IB)+B(IADD2+ID)
               T3REAL=A(IADD2)-A(IADD2+IC)
               T3IMAG=B(IADD2)-B(IADD2+IC)
               T4REAL=A(IADD2+IB)-A(IADD2+ID)
               T4IMAG=B(IADD2+IB)-B(IADD2+ID)
               C(IADD1)   =T1REAL+T2REAL
               D(IADD1)   =T1IMAG+T2IMAG
               C(IADD1+JB)=T3REAL-T4IMAG
               D(IADD1+JB)=T3IMAG+T4REAL
               C(IADD1+JC)=T1REAL-T2REAL
               D(IADD1+JC)=T1IMAG-T2IMAG
               C(IADD1+JD)=T3REAL+T4IMAG
               D(IADD1+JD)=T3IMAG-T4REAL
             ENDDO
             I=I+ISKIP
             J=J+ISKIP1
           ENDDO
           J=J+JUMP
           DO K=LA,M-LA,LA
             K2=K+K
             KK=K2+1
             KK2=KK+K2
             KK4=KK2+K2
             DO L=1,LA
               IADD1=-ISKIP3+J
*vdir nodep             
CDIR$ IVDEP
CDEC$ IVDEP
               DO NU=0,NUMBER-1
                 IADD1=IADD1+ISKIP3
                 IADD2=NU+I
                 T1REAL=A(IADD2)+A(IADD2+IC)
                 T1IMAG=B(IADD2)+B(IADD2+IC)
                 T2REAL=A(IADD2+IB)+A(IADD2+ID)
                 T2IMAG=B(IADD2+IB)+B(IADD2+ID)
                 T3REAL=A(IADD2)-A(IADD2+IC)
                 T3IMAG=B(IADD2)-B(IADD2+IC)
                 T4REAL=A(IADD2+IB)-A(IADD2+ID)
                 T4IMAG=B(IADD2+IB)-B(IADD2+ID)
                 X1REAL=T3REAL-T4IMAG
                 X1IMAG=T3IMAG+T4REAL
                 X2REAL=T1REAL-T2REAL
                 X2IMAG=T1IMAG-T2IMAG
                 X3REAL=T3REAL+T4IMAG
                 X3IMAG=T3IMAG-T4REAL
                 C(IADD1)   =T1REAL+T2REAL
                 D(IADD1)   =T1IMAG+T2IMAG
                 C(IADD1+JB)=TRIGS(KK)*X1REAL-TRIGI(KK)*X1IMAG
                 D(IADD1+JB)=TRIGI(KK)*X1REAL+TRIGS(KK)*X1IMAG
                 C(IADD1+JC)=TRIGS(KK2)*X2REAL-TRIGI(KK2)*X2IMAG
                 D(IADD1+JC)=TRIGI(KK2)*X2REAL+TRIGS(KK2)*X2IMAG
                 C(IADD1+JD)=TRIGS(KK4)*X3REAL-TRIGI(KK4)*X3IMAG
                 D(IADD1+JD)=TRIGI(KK4)*X3REAL+TRIGS(KK4)*X3IMAG
               ENDDO
               I=I+ISKIP
               J=J+ISKIP1
             ENDDO
             J=J+JUMP
           ENDDO
         ELSE
           DO L=1,LA
             IADD1=-ISKIP3+J
*vdir nodep             
CDIR$ IVDEP
CDEC$ IVDEP
             DO NU=0,NUMBER-1
               IADD1=IADD1+ISKIP3
               IADD2=NU+I
               T1REAL=A(IADD2)+A(IADD2+IC)
               T1IMAG=B(IADD2)+B(IADD2+IC)
               T2REAL=A(IADD2+IB)+A(IADD2+ID)
               T2IMAG=B(IADD2+IB)+B(IADD2+ID)
               T3REAL=A(IADD2)-A(IADD2+IC)
               T3IMAG=B(IADD2)-B(IADD2+IC)
               T4REAL=A(IADD2+IB)-A(IADD2+ID)
               T4IMAG=B(IADD2+IB)-B(IADD2+ID)
               C(IADD1)   =T1REAL+T2REAL
               D(IADD1)   =T1IMAG+T2IMAG
               C(IADD1+JB)=T3REAL+T4IMAG
               D(IADD1+JB)=T3IMAG-T4REAL
               C(IADD1+JC)=T1REAL-T2REAL
               D(IADD1+JC)=T1IMAG-T2IMAG
               C(IADD1+JD)=T3REAL-T4IMAG
               D(IADD1+JD)=T3IMAG+T4REAL
             ENDDO
             I=I+ISKIP
             J=J+ISKIP1
           ENDDO
           J=J+JUMP
           DO K=LA,M-LA,LA
             K2=K+K
             KK=K2+1
             KK2=KK+K2
             KK4=KK2+K2
             DO L=1,LA
               IADD1=-ISKIP3+J
*vdir nodep             
CDIR$ IVDEP
CDEC$ IVDEP
               DO NU=0,NUMBER-1
                 IADD1=IADD1+ISKIP3
                 IADD2=NU+I
                 T1REAL=A(IADD2)+A(IADD2+IC)
                 T1IMAG=B(IADD2)+B(IADD2+IC)
                 T2REAL=A(IADD2+IB)+A(IADD2+ID)
                 T2IMAG=B(IADD2+IB)+B(IADD2+ID)
                 T3REAL=A(IADD2)-A(IADD2+IC)
                 T3IMAG=B(IADD2)-B(IADD2+IC)
                 T4REAL=A(IADD2+IB)-A(IADD2+ID)
                 T4IMAG=B(IADD2+IB)-B(IADD2+ID)
                 X1REAL=T3REAL+T4IMAG
                 X1IMAG=T3IMAG-T4REAL
                 X2REAL=T1REAL-T2REAL
                 X2IMAG=T1IMAG-T2IMAG
                 X3REAL=T3REAL-T4IMAG
                 X3IMAG=T3IMAG+T4REAL
                 C(IADD1)   = T1REAL+T2REAL
                 D(IADD1)   = T1IMAG+T2IMAG
                 C(IADD1+JB)= TRIGS(KK) *X1REAL+TRIGI(KK)*X1IMAG
                 D(IADD1+JB)=-TRIGI(KK)*X1REAL+TRIGS(KK) *X1IMAG
                 C(IADD1+JC)= TRIGS(KK2)*X2REAL+TRIGI(KK2)*X2IMAG
                 D(IADD1+JC)=-TRIGI(KK2)*X2REAL+TRIGS(KK2)*X2IMAG
                 C(IADD1+JD)= TRIGS(KK4)*X3REAL+TRIGI(KK4)*X3IMAG
                 D(IADD1+JD)=-TRIGI(KK4)*X3REAL+TRIGS(KK4)*X3IMAG
               ENDDO
               I=I+ISKIP
               J=J+ISKIP1
             ENDDO
             J=J+JUMP
           ENDDO
         ENDIF
       ELSEIF (IFAC .EQ. 3) THEN
C
C        FACTOR THREE 
C
         IC=IB+MSKIP
         JC=JB+LASKIP
         IF(ISIGN.LT.0) THEN
           DO L=1,LA
             IADD1=-ISKIP3+J
*vdir nodep             
CDIR$ IVDEP
CDEC$ IVDEP
             DO NU=0,NUMBER-1
               IADD1=IADD1+ISKIP3
               IADD2=NU+I
               T1REAL=A(IADD2+IB)+A(IADD2+IC)
               T1IMAG=B(IADD2+IB)+B(IADD2+IC)
               T2REAL=A(IADD2)-PO5*T1REAL
               T2IMAG=B(IADD2)-PO5*T1IMAG
               T3REAL=SIN60*(A(IADD2+IB)-A(IADD2+IC))
               T3IMAG=SIN60*(B(IADD2+IB)-B(IADD2+IC))
               C(IADD1)   =A(IADD2)+T1REAL
               D(IADD1)   =B(IADD2)+T1IMAG
               C(IADD1+JB)=T2REAL-T3IMAG
               D(IADD1+JB)=T2IMAG+T3REAL
               C(IADD1+JC)=T2REAL+T3IMAG
               D(IADD1+JC)=T2IMAG-T3REAL
             ENDDO
             I=I+ISKIP
             J=J+ISKIP1
           ENDDO
           J=J+JUMP
           DO K=LA,M-LA,LA
             K2=K+K
             KK=K2+1
             KK2=KK+K2
             DO L=1,LA
               IADD1=-ISKIP3+J
*vdir nodep             
CDIR$ IVDEP
CDEC$ IVDEP
               DO NU=0,NUMBER-1
                 IADD1=IADD1+ISKIP3
                 IADD2=NU+I
                 T1REAL=A(IADD2+IB)+A(IADD2+IC)
                 T1IMAG=B(IADD2+IB)+B(IADD2+IC)
                 T2REAL=A(IADD2)-PO5*T1REAL
                 T2IMAG=B(IADD2)-PO5*T1IMAG
                 T3REAL=SIN60*(A(IADD2+IB)-A(IADD2+IC))
                 T3IMAG=SIN60*(B(IADD2+IB)-B(IADD2+IC))
                 X1REAL=T2REAL-T3IMAG
                 X1IMAG=T2IMAG+T3REAL
                 X2REAL=T2REAL+T3IMAG
                 X2IMAG=T2IMAG-T3REAL
                 C(IADD1)   =A(IADD2)+T1REAL
                 D(IADD1)   =B(IADD2)+T1IMAG
                 C(IADD1+JB)=TRIGS(KK)*X1REAL-TRIGI(KK)*X1IMAG
                 D(IADD1+JB)=TRIGI(KK)*X1REAL+TRIGS(KK)*X1IMAG
                 C(IADD1+JC)=TRIGS(KK2)*X2REAL-TRIGI(KK2)*X2IMAG
                 D(IADD1+JC)=TRIGI(KK2)*X2REAL+TRIGS(KK2)*X2IMAG
               ENDDO
               I=I+ISKIP
               J=J+ISKIP1
             ENDDO
             J=J+JUMP
           ENDDO
         ELSE
           DO L=1,LA
             IADD1=-ISKIP3+J
*vdir nodep             
CDIR$ IVDEP
CDEC$ IVDEP
             DO NU=0,NUMBER-1
               IADD1=IADD1+ISKIP3
               IADD2=NU+I
               T1REAL=A(IADD2+IB)+A(IADD2+IC)
               T1IMAG=B(IADD2+IB)+B(IADD2+IC)
               T2REAL=A(IADD2)-PO5*T1REAL
               T2IMAG=B(IADD2)-PO5*T1IMAG
               T3REAL=SIN60*(A(IADD2+IB)-A(IADD2+IC))
               T3IMAG=SIN60*(B(IADD2+IB)-B(IADD2+IC))
               C(IADD1)   =A(IADD2)+T1REAL
               D(IADD1)   =B(IADD2)+T1IMAG
               C(IADD1+JB)=T2REAL+T3IMAG
               D(IADD1+JB)=T2IMAG-T3REAL
               C(IADD1+JC)=T2REAL-T3IMAG
               D(IADD1+JC)=T2IMAG+T3REAL
             ENDDO
             I=I+ISKIP
             J=J+ISKIP1
           ENDDO
           J=J+JUMP
           DO K=LA,M-LA,LA
             K2=K+K
             KK=K2+1
             KK2=KK+K2
             DO L=1,LA
               IADD1=-ISKIP3+J
*vdir nodep             
CDIR$ IVDEP
CDEC$ IVDEP
               DO NU=0,NUMBER-1
                 IADD1=IADD1+ISKIP3
                 IADD2=NU+I
                 T1REAL=A(IADD2+IB)+A(IADD2+IC)
                 T1IMAG=B(IADD2+IB)+B(IADD2+IC)
                 T2REAL=A(IADD2)-PO5*T1REAL
                 T2IMAG=B(IADD2)-PO5*T1IMAG
                 T3REAL=SIN60*(A(IADD2+IB)-A(IADD2+IC))
                 T3IMAG=SIN60*(B(IADD2+IB)-B(IADD2+IC))
                 X1REAL=T2REAL+T3IMAG
                 X1IMAG=T2IMAG-T3REAL
                 X2REAL=T2REAL-T3IMAG
                 X2IMAG=T2IMAG+T3REAL
                 C(IADD1)   =A(IADD2)+T1REAL
                 D(IADD1)   =B(IADD2)+T1IMAG
                 C(IADD1+JB)= TRIGS(KK)*X1REAL+TRIGI(KK)*X1IMAG
                 D(IADD1+JB)=-TRIGI(KK)*X1REAL+TRIGS(KK)*X1IMAG
                 C(IADD1+JC)= TRIGS(KK2)*X2REAL+TRIGI(KK2)*X2IMAG
                 D(IADD1+JC)=-TRIGI(KK2)*X2REAL+TRIGS(KK2)*X2IMAG
               ENDDO
               I=I+ISKIP
               J=J+ISKIP1
             ENDDO
             J=J+JUMP
           ENDDO
         ENDIF
       ELSEIF (IFAC .EQ. 2) THEN
C
C        FACTOR TWO 
C
         DO L=1,LA
           IADD1=-ISKIP3+J
*vdir nodep             
CDIR$ IVDEP
CDEC$ IVDEP
           DO NU=0,NUMBER-1
             IADD1=IADD1+ISKIP3
             IADD2=NU+I
             C(IADD1)   =A(IADD2)+A(IADD2+IB)
             D(IADD1)   =B(IADD2)+B(IADD2+IB)
             C(IADD1+JB)=A(IADD2)-A(IADD2+IB)
             D(IADD1+JB)=B(IADD2)-B(IADD2+IB)
           ENDDO
           I=I+ISKIP
           J=J+ISKIP1
         ENDDO
         J=J+JUMP
         DO K=LA,M-LA,LA
           KK=2*K+1    
           TIMAG=TRIGI(KK)
           IF(ISIGN.GT.0) TIMAG=-TIMAG
           DO L=1,LA
             IADD1=-ISKIP3+J
*vdir nodep             
CDIR$ IVDEP
CDEC$ IVDEP
             DO NU=0,NUMBER-1
               IADD1=IADD1+ISKIP3
               IADD2=NU+I
               X1REAL=A(IADD2)-A(IADD2+IB)
               X1IMAG=B(IADD2)-B(IADD2+IB)
               C(IADD1)   =A(IADD2)+A(IADD2+IB)
               D(IADD1)   =B(IADD2)+B(IADD2+IB)
               C(IADD1+JB)=TRIGS(KK)*X1REAL-TIMAG*X1IMAG
               D(IADD1+JB)=TIMAG*X1REAL+TRIGS(KK)*X1IMAG
             ENDDO
             I=I+ISKIP
             J=J+ISKIP1
            ENDDO
           J=J+JUMP
         ENDDO
       ENDIF
       RETURN
       END
