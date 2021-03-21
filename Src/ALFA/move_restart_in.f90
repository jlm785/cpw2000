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

!>     Reads the restart file.
!>     For molecular dynamics trajectory should be the same.
!>     for l_bfgs minimization, algorithm is restarted.

       subroutine move_restart_in(flgcal,                                     &
     & ntype,natom,nameat,atmass,rat,vat,adot,vadot,                     &
     & istmd,tstep,beta,tempk,iseed,strext,press,celmas,                 &
     & rat1,frc1,adot1,frcel1,                                           &
     & mxdatm,mxdtyp)

!      Written 26 may 99. jlm
!      Rewritten 19 november 2002. jlm
!      Modified 14 october 2003. jlm
!      Modified (formats) 17 august 2004. jlm
!      Modified 6 January 2017, f90. JLM
!      Modified, documentation, August 2019. JLM
!      Copyright inesc-mn/Jose Luis Martins/Benedito Costa Cabral

!      version 4.94 of pw


       implicit none

       integer, parameter :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdtyp                     !<  array dimension of types of atoms
       integer, intent(in)                ::  mxdatm                     !<  array dimension of number of atoms of a given type

       integer, intent(in)                ::  ntype                      !<  number of types of atoms
       integer, intent(in)                ::  natom(mxdtyp)              !<  number of atoms of type i
       character(len=2), intent(in)       ::  nameat(mxdtyp)             !<  chemical symbol for the type i

!      input and output

       character(len=6), intent(inout)    ::  flgcal

!      output

       real(REAL64), intent(out)          ::  atmass(mxdtyp)             !<  atomic mass of atoms of type i

       real(REAL64), intent(out)          ::  rat(3,mxdatm,mxdtyp)       !<  lattice coordinates of atom j of type i
       real(REAL64), intent(out)          ::  vat(3,mxdatm,mxdtyp)       !<  d rat / d t  velocity in lattice coordinates of atom j of type i
       real(REAL64), intent(out)          ::  adot(3,3)                  !<  metric in real space
       real(REAL64), intent(out)          ::  vadot(3,3)                 !<  d adot / d t  rate of change of metric

       integer,intent(out)                ::  istmd                      !<  md step. Equal to 1 in first step of molecular dynamics
       real(REAL64), intent(out)          ::  tstep                      !<  molecular dynamics time step (in a.u.)

       real(REAL64), intent(out)          ::  tempk                      !<  ionic temperature (in Kelvin)
       real(REAL64), intent(out)          ::  beta                       !<  friction coefficient/mass (in a.u.)

       integer, intent(out)               ::  iseed                      !<  seed for random number generator

       real(REAL64), intent(out)          ::  strext(3,3)                !<  external applied stress
       real(REAL64), intent(out)          ::  press                      !<  external pressure
       real(real64), intent(out)          ::  celmas                     !<  fictitious cell mass       

       real(REAL64), intent(out)          ::  rat1(3,mxdatm,mxdtyp)      !<  previous value of lattice coordinates of atom j of type i
       real(REAL64), intent(out)          ::  frc1(3,mxdatm,mxdtyp)      !<  previous value of force on atom j of type i
       real(REAL64), intent(out)          ::  adot1(3,3)                 !<  previous value of adot
       real(REAL64), intent(out)          ::  frcel1(3,3)                !<  previous cell "force" (covariant components)

!      local variables

       integer              ::  ntyold, natold
       character(len=6)     ::  flgold
       character(len=2)     ::  namold

!      counters

       integer       ::  nt, i, j

       open(unit = 7,file = 'RESTART.DAT',status = 'unknown',            &
     &      form = 'formatted')

       read(7,'(2x,a6)') flgold

       read(7,'(9(2x,e24.16))') ((adot(i,j),i=1,3),j=1,3)
       read(7,'(2x,i10)') ntyold

       if(ntyold /= ntype) then
         write(6,'("   STOPPED in restart_in:    old and new number",    &
     &        " of types of atoms are",2i7)') ntyold, ntype

         stop

       endif

       do nt = 1,ntype
         read(7,'(2x,i10,2x,a2,2x,e24.16)') natold,namold,atmass(nt)

         if(natold /= natom(nt) .or. namold /= nameat(nt)) then
           write(6,'("   STOPPED in restart_in:  for atom number ",i4,   &
     &        " the old and new number of atoms and symbols are",        &
     &        2i7,2x,a2,2x,a2)') nt,natold,natom(nt),namold,nameat(nt)

           stop

         endif

         do i = 1,natom(nt)
           read(7,'(3(2x,e24.16))') (rat(j,i,nt),j=1,3)
         enddo
       enddo

       if(flgcal == 'RSTRT ') flgcal = flgold

       if(flgcal == 'MICRO ') then

         read(7,'(2x,i10,2x,e24.16)') istmd,tstep
         do nt=1,ntype
           do i=1,natom(nt)
             read(7,'(9(2x,e24.16))') (vat(j,i,nt),j=1,3),               &
     &                   (rat1(j,i,nt),j=1,3),                           &
     &                   (frc1(j,i,nt),j=1,3)
           enddo
         enddo

       elseif(flgcal == 'LANG  ') then

         read(7,'(2x,i10,3(2x,e24.16),2x,i10)') istmd,tstep,beta,        &
     &           tempk,iseed
         do nt = 1,ntype
           do i = 1,natom(nt)
             read(7,'(11(2x,e24.16))') (vat(j,i,nt),j=1,3),              &
     &                   (rat1(j,i,nt),j=1,3),(frc1(j,i,nt),j=1,3)
           enddo
         enddo

       elseif(flgcal == 'VCSLNG' .or. flgcal == 'EPILNG') then

         read(7,'(2x,i10,3(2x,e24.16),2x,i10)') istmd,tstep,beta,        &
     &           tempk,iseed
         read(7,'(11(2x,e24.16))') press,((strext(i,j),i=1,3),j=1,3),    &
     &           celmas
         do nt = 1,ntype
           do i = 1,natom(nt)
             read(7,'(9(2x,e24.16))') (vat(j,i,nt),j=1,3),               &
     &                   (rat1(j,i,nt),j=1,3),(frc1(j,i,nt),j=1,3)
           enddo
         enddo
         read(7,'(9(2x,e24.16))') ((vadot(i,j),i=1,3),j=1,3)
         read(7,'(9(2x,e24.16))') ((adot1(i,j),i=1,3),j=1,3)
         read(7,'(9(2x,e24.16))') ((frcel1(i,j),i=1,3),j=1,3)

       elseif(flgcal == 'VCSMIC') then

         read(7,'(2x,i10,2x,e24.16)') istmd,tstep
         read(7,'(11(2x,e24.16))') press,((strext(i,j),i=1,3),j=1,3),    &
     &           celmas
         do nt = 1,ntype
           do i = 1,natom(nt)
             read(7,'(9(2x,e24.16))') (vat(j,i,nt),j=1,3),               &
     &                   (rat1(j,i,nt),j=1,3),(frc1(j,i,nt),j=1,3)
           enddo
         enddo
         read(7,'(9(2x,e24.16))') ((vadot(i,j),i=1,3),j=1,3)
         read(7,'(9(2x,e24.16))') ((adot1(i,j),i=1,3),j=1,3)
         read(7,'(9(2x,e24.16))') ((frcel1(i,j),i=1,3),j=1,3)

       elseif(flgcal == 'VCSLBF' .or. flgcal == 'EPILBF') then

         read(7,'(11(2x,e24.16))') press,((strext(i,j),i=1,3),j=1,3)

       elseif(flgcal == 'LBFSYM') then

!      do nothing

       else
         write(6,'("   STOPPED in restart_in:  do not know how to",      &
     &        " restart.  flgcal = ",a6)') flgcal

         stop

       endif

       close(unit=7)

       return
       end subroutine move_restart_in
