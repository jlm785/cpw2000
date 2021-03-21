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

!>     Reads the crystal structure from a data file
!>     compatible with Berkeley/Sverre data files.
!>     Hence the complicated fixed format structure.

       subroutine read_data(adot,ntype,natom,nameat,rat,atmass,alatt,    &
     & mxdtyp,mxdatm)

!      Adapted from Sverre Froyen plane wave program
!      Written 15 October 1993. jlm
!      Modified 15 January 1999. jlm
!      Modified 2013/2014. CLR
!      Modified,f90, near_rational, 30 May 2014. JLM
!      Modified, near_rational bug trap, 21 June 2014. JLM
!      Modified, ntry, printing, STRUCT, December 8, 2016. JLM
!      Modified, documentation, August 2019. JLM
!      copyright INESC-MN/Jose Luis Martins/Sverre Froyen/Carlos Loia


!      version 4.94


       implicit none
       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdtyp                     !<  array dimension of types of atoms
       integer, intent(in)                ::  mxdatm                     !<  array dimension of types of atoms

!      output

       real(REAL64), intent(out)          ::  adot(3,3)                  !<  metric in direct space
       integer, intent(out)               ::  ntype                      !<  number of types of atoms
       integer, intent(out)               ::  natom(mxdtyp)              !<  number of atoms of type i
       character(len=2), intent(out)      ::  nameat(mxdtyp)             !<  chemical symbol for the type i
       real(REAL64), intent(out)          ::  rat(3,mxdatm,mxdtyp)       !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i
       real(REAL64), intent(out)          ::  atmass(mxdtyp)             !<  atomic mass (in a.u.) of atom of type i 

       real(REAL64), intent(out)          ::  alatt                      !<  lattice constant

!      other variables

       real(REAL64)      ::  avec(3,3)

       character(len=1)  ::  iop(3,3),jop(3)
       real(REAL64)      ::  car(3)
       real(REAL64)      ::  sgn
       real(REAL64)      ::  tmp
       real(REAL64)      ::  difmax
       
       logical           ::  lf                     !  an approximation was found
       logical           ::  lsq                    !  it is the square of a rati onal
       integer           ::  nnum, nden                 !  x ~ nnum/ndem or x**2 ~ nnum/nden and sign of nnum is the same as x
       integer           ::  ntry                       !  tries denominators up to ntry


!      constants

       real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
       real(REAL64), parameter  :: AMU = 1822.888485_REAL64
       real(REAL64), parameter  :: EPS = 2.0E-6                          !  Criteria for rounding the values (depends on the f12.6 in input)
       real(REAL64), parameter  :: EPSMALL = 1.0E-14                     !  Criteria for worthwhile rounding

!      counters

       integer    ::  i, j, ja, jmax, nt, ntt

!      hardcoded value of ntry

       ntry = 36

!      read scale factor

       read(5,'(f12.6)') alatt

       write(6,*)
       write(6,*)
       write(6,*)'  Crystal Structure:'
       write(6,*)
       write(6,*)
       write(6,'(4x,f14.6,"   reference lattice constant")') alatt
       write(6,*)

!      read the basis vectors for the lattice (cartesian coordinates and
!      atomic units divided by ascale). modify according to iop...
!      multiply by alatt.

       do j=1,3
         read(5,'(3(2x,a1,f12.6))') (iop(i,j),avec(i,j),i=1,3)
       enddo

       do j=1,3
       do i=1,3
         sgn = UM
         if (avec(i,j) < ZERO) sgn = -sgn
         avec(i,j) = abs(avec(i,j))
         if (iop(i,j)  ==  'S') avec(i,j) = sqrt(avec(i,j))
         if (iop(i,j)  ==  'C') avec(i,j) = avec(i,j)**(UM/3)
         if (iop(i,j)  ==  'T') avec(i,j) = avec(i,j)/3
         if (iop(i,j)  ==  'H') avec(i,j) = avec(i,j)/sqrt(3*UM)
         avec(i,j) = sgn*avec(i,j)
       enddo
       enddo

 !      checks if it is close to a rational number or the sqrt of a rational

       do i=1,3
       do j=1,3
         
         call near_rational(avec(j,i),lf,lsq,nnum,nden,ntry,5.0*EPS)

         if(lf) then
           if(lsq) then
             tmp = sign(sqrt(UM*nnum/(UM*nden)),avec(j,i))
           else
             tmp = UM*nnum/(UM*nden)
           endif
           if(abs(tmp - avec(j,i)) < EPS .and.                           &
     &        abs(tmp - avec(j,i)) > EPSMALL) then
             write(6,*)
             write(6,'("   WARNING:   avec(",i3,",",i3,")")') j,i
             write(6,'("   changed from ",f20.12," to ",f20.12)')        &
     &                avec(j,i),tmp
             write(6,*)
             avec(j,i) = tmp
           endif
         endif
         
       enddo
       enddo

       do j=1,3
       do i=1,3
         avec(i,j) = alatt*avec(i,j)
       enddo
       enddo

       write(6,*)
       write(6,*) '  Primitive Translation Vectors'
       write(6,'(24x,"in a.u.",28x,"in lattice units")')

       do j=1,3
         write(6,'("  a",i1,"=",3(1x,e13.6),5x,3(2x,f7.3))')             &
     &                j,(avec(i,j),i=1,3),(avec(i,j)/alatt,i=1,3)
       enddo

!      read the number of atoms of each type and their names.
!      read atomic positions in unit cell (in units of
!      basis vectors). modify according to iop.

       read(5,'(i5)') ntype

       if(ntype > mxdtyp .or. ntype <= 0) then
         write(6,*)
         write(6,'("    STOPPED  in read_data:   ",                      &
     &        "  increase mxdtyp to ",i6)') ntype

         stop

       endif

       do nt=1,ntype

         read(5,'(i5,3x,a2)') natom(nt),nameat(nt)

         if(natom(nt) > mxdatm .or. natom(nt) <= 0) then
           write(6,*)
           write(6,'("    STOPPED in read_data:   ",                     &
     &        "  increase mxdatm to ",i6)') natom(nt)

           stop

         endif

         call p_tbl_mass(nameat(nt),atmass(nt))

         jmax = natom(nt)
         do ja=1,jmax
           read(5,'(3(2x,a1,f12.6))') (jop(i),rat(i,ja,nt),i=1,3)
           do i=1,3
             sgn = UM
             if (rat(i,ja,nt) < ZERO) sgn = -sgn
             rat(i,ja,nt) = abs(rat(i,ja,nt))
             if (jop(i)  ==  'S') rat(i,ja,nt) = sqrt(rat(i,ja,nt))
             if (jop(i)  ==  'T')  rat(i,ja,nt) = rat(i,ja,nt)/3
             rat(i,ja,nt) = sgn*rat(i,ja,nt)
           enddo
         enddo
       enddo

!      checks if it is close to a rational number or the sqrt of a rational

       difmax = ZERO
  
       do nt=1,ntype
       do ja=1,natom(nt)
       do i=1,3

          call near_rational(rat(i,ja,nt),lf,lsq,nnum,nden,ntry,5.0*EPS)

         if(lf) then
           if(lsq) then
             tmp = sign(sqrt(UM*nnum/(UM*nden)),rat(i,ja,nt))
           else
             tmp = UM*nnum/(UM*nden)
           endif
           if(abs(tmp - rat(i,ja,nt)) < EPS .and.                        &
     &        abs(tmp - rat(i,ja,nt)) > EPSMALL) then
!             write(6,*)
!             write(6,'("   WARNING:    rat(",i3,",",i3,",",i3,")")')     &
!     &                i,ja,nt
!             write(6,'("   changed from ",f20.12," to ",f20.12)')        &
!     &              rat(i,ja,nt),tmp
!             write(6,*)
             difmax = max(difmax,abs(tmp - rat(i,ja,nt)))
             rat(i,ja,nt) = tmp
           endif
         endif
 
       enddo
       enddo
       enddo
       
       if(difmax > sqrt(EPS*EPSMALL)) then
             write(6,*)
             write(6,'("   WARNING:    atomic positions changed in ",    &
     &           "read_data")')
             write(6,'("   Maximum correction ",f20.12)') difmax
             write(6,*)
       endif

!      printout

       write(6,*)
       write(6,'("  No. Type     Position(lattice coord.)",10x,          &
     &          "Position(cartesian coord.)",15x,"Mass")')
       write(6,*)

       ntt = 0
       do nt=1,ntype
         jmax = natom(nt)
         do ja=1,jmax
           ntt = ntt + 1
           car(1) = avec(1,1)*rat(1,ja,nt) + avec(1,2)*rat(2,ja,nt) +    &
     &              avec(1,3)*rat(3,ja,nt)
           car(2) = avec(2,1)*rat(1,ja,nt) + avec(2,2)*rat(2,ja,nt) +    &
     &              avec(2,3)*rat(3,ja,nt)
           car(3) = avec(3,1)*rat(1,ja,nt) + avec(3,2)*rat(2,ja,nt) +    &
     &              avec(3,3)*rat(3,ja,nt)

           write(6,'(1x,i3,3x,a2,2x,3(2x,f7.3),3x,3(1x,e13.6),5x,f8.3)') &
     &               ntt,nameat(nt),(rat(i,ja,nt),i=1,3),                &
     &               (car(i),i=1,3),atmass(nt)

         enddo
       enddo

!      convert mass

       do nt=1,ntype
         atmass(nt) = AMU*atmass(nt)
       enddo

!      calculate metric

       adot(1,1) = avec(1,1)*avec(1,1) + avec(2,1)*avec(2,1) +           &
     &             avec(3,1)*avec(3,1)
       adot(2,2) = avec(1,2)*avec(1,2) + avec(2,2)*avec(2,2) +           &
     &             avec(3,2)*avec(3,2)
       adot(3,3) = avec(1,3)*avec(1,3) + avec(2,3)*avec(2,3) +           &
     &             avec(3,3)*avec(3,3)
       adot(1,2) = avec(1,1)*avec(1,2) + avec(2,1)*avec(2,2) +           &
     &             avec(3,1)*avec(3,2)
       adot(1,3) = avec(1,1)*avec(1,3) + avec(2,1)*avec(2,3) +           &
     &             avec(3,1)*avec(3,3)
       adot(2,3) = avec(1,2)*avec(1,3) + avec(2,2)*avec(2,3) +           &
     &             avec(3,2)*avec(3,3)
       adot(2,1) = adot(1,2)
       adot(3,1) = adot(1,3)
       adot(3,2) = adot(2,3)
       
!       lfile = .FALSE.
!       inquire(file='STRUCT.DAT',exist=lfile)
!       if (lfile) then
!         open(unit=222,file='STRUCT.DAT',form='unformatted')
!         write(6,*) ' In read_data str lfile=', lfile
!       endif
!       if(lfile) then
!         write(6,*) ' Found struct.dat file!' 
!         write(6,*) ' Metric and atom positions will be overrriden.'
!         write(6,*)
!         do i=1,3
!           read(222) (adot(ja,i),ja=1,3)
!         enddo
!         read(222) ntype
!         do nt=1,ntype
!           read(222) natom(nt)
!           do i=1,natom(nt)
!             read(222)(rat(k,i,nt),k=1,3)
!           enddo
!         enddo       
!         close(unit = 222)
!       endif
       

       return

       end subroutine read_data


