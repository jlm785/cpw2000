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

!>     Reads the parameters of the calculation and descritpion of the crystal

       subroutine pw_rho_v_in_crystal_calc(io,                           &
     &        pwline,title,subtitle,meta_cpw2000,                        &
     &        author,flgscf,flgdal,emax,teleck,                          &
     &        nx,ny,nz,sx,sy,sz,nband,                                   &
     &        ng,ns,                                                     &
     &        ntrans, mtrx, tnp,                                         &
     &        alatt,adot,ntype,natom,nameat,rat,                         &
     &        mxdtyp,mxdatm,mxdgve,mxdnst,mxdlqp)

!      This sunx,ny,nzbroutines reads the first part of file io,
!      prints sx,sy,szcalculation flags other information
!      and returns the crystal geometry.
!      If allocated space is not correct stops the calculation.

!      Written May 26 2014. JLM
!      Splitted from code written January 12, 2014. JLM
!      Modified pwline,title,subtitle, 1 August 2014. JLM
!      Modified, meta_cpw2000, author,nx,etc, January 10, 2017. JLM
!      Modified, documentation, spacegroup, 31 December 2020. JLM
!      copyright  Jose Luis Martins/INESC-MN

!      version 4.99

       implicit none
       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdtyp                     !<  array dimension of types of atoms
       integer, intent(in)                ::  mxdatm                     !<  array dimension of types of atoms
       integer, intent(in)                ::  mxdgve                     !<  array dimension for g-space vectors
       integer, intent(in)                ::  mxdnst                     !<  array dimension for g-space stars
       integer, intent(in)                ::  mxdlqp                     !<  array dimension for local potential

       integer, intent(in)                ::  io                         !<  number of tape to which the pseudo is added.

!      output

       character(len=60), intent(out)     ::  pwline                     !<  identification of the calculation
       character(len=50), intent(out)     ::  title                      !<  title for plots
       character(len=140), intent(out)    ::  subtitle                   !<  title for plots
       character(len=250), intent(out)    ::  meta_cpw2000               !<  metadata from cpw2000

       integer, intent(out)               ::  ntrans                     !<  number of symmetry operations in the factor group
       integer, intent(out)               ::  mtrx(3,3,48)               !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group
       real(REAL64), intent(out)          ::  tnp(3,48)                  !<  2*pi* i-th component (in lattice coordinates) of the fractional translation vector associated with the k-th symmetry operation of the factor group

       real(REAL64), intent(out)          ::  adot(3,3)                  !<  metric in direct space
       integer, intent(out)               ::  ntype                      !<  number of types of atoms
       integer, intent(out)               ::  natom(mxdtyp)              !<  number of atoms of type i
       character(len=2), intent(out)      ::  nameat(mxdtyp)             !<  chemical symbol for the type i
       real(REAL64), intent(out)          ::  rat(3,mxdatm,mxdtyp)       !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i

       integer, intent(out)               ::  ng                         !<  total number of g-vectors with length less than gmax
       integer, intent(out)               ::  ns                         !<  number os stars with length less than gmax

!      information about the calculation

       character(len=3), intent(out)      ::  author                     !<  type of xc wanted (CA=PZ , PW92 , PBE)
       character(len=6), intent(out)      ::  flgscf                     !<  type of self consistent field and diagonalization
       character(len=4), intent(out)      ::  flgdal                     !<  whether the dual approximation is used
       real(REAL64), intent(out)          ::  emax                       !<  kinetic energy cutoff of plane wave expansion (Hartree).
       real(REAL64), intent(out)          ::  teleck                     !<  electronic temperature (in Kelvin)

       integer, intent(out)               ::  nband                      !<  target for number of bands
       integer, intent(out)               ::  nx,ny,nz                   !<  divisions of Brillouin zone for integration (Monkhorst-Pack)
       real(REAL64), intent(out)          ::  sx,sy,sz                   !<  shift of points in division of Brillouin zone for integration (Monkhorst-Pack)
       real(REAL64), intent(out)          ::  alatt                      !<  lattice constant

!      other local variables

       character(len=9 )   ::  bdate
       character(len=8)    ::  btime
       integer             ::  mxdl
       character(len=140)  ::  line140
       character(len=250)  ::  line250
       integer             ::  ioerr, ioerr2

!      counters

       integer    ::  i, j, k


!      reads the first record from the file

       read(io,iostat = ioerr) ntype,ng,ns,mxdl, ntrans
       if(ioerr == 0) then
         if(ntrans > 48 .or. ntrans < 1) then
           write(6,*)
           write(6,'("  STOPPED in in_v_rho")')
           write(6,'("  ntrans = ",i8)') ntrans

           stop

         endif
       else
         backspace(io)
         read(io) ntype,ng,ns,mxdl
         ntrans = 0
       endif

       if(ntype > mxdtyp) then
         write(6,*)
         write(6,'("  STOPPED in in_v_rho")')
         write(6,'("  ntype = ",i8," is greater than mxdtyp =",i8)')     &
     &        ntype,mxdtyp

         stop

       endif

       if(ng > mxdgve) then
         write(6,*)
         write(6,'("  STOPPED in in_v_rho")')
         write(6,'("  ng = ",i8," is greater than mxdgve =",i8)')        &
     &        ng,mxdgve

         stop

       endif

       if(ns > mxdnst) then
         write(6,*)
         write(6,'("  STOPPED in in_v_rho")')
         write(6,'("  ns = ",i8," is greater than mxdnst =",i8)')        &
     &        ns,mxdnst

         stop

       endif

       if(mxdl > mxdlqp) then
         write(6,*)
         write(6,'("  STOPPED in in_v_rho")')
         write(6,'("  mxdl = ",i8," is greater than mxdlqp =",i8)')      &
     &        ns,mxdnst

         stop

       endif

       read(io) (natom(i),i=1,ntype)
       do i=1,ntype
         if(natom(i) > mxdatm) then
           write(6,*)
           write(6,'("  STOPPED in in_v_rho")')
           write(6,'("  for atom type ",i4," natom = ",i8," is ",        &
     &       "greater than mxdatm =",i8)') i,natom(i),mxdatm

           stop

         endif
       enddo

       read(io) bdate,btime
       write(6,'("  The data file was written on ",a9," at ",a8)')       &
     &             bdate,btime

       read(io) author,flgscf,flgdal
       if(author == 'ca ' .or. author == 'CA ' .or. author == 'pz '      &
     &         .or. author == 'PZ ') then
         write(6,*)
         write(6,'("  The potential was calculated in the local ",       &
     &     "density aproximation using Ceperley and Alder ",             &
     &     "correlation")')
         write(6,'("  (as parametrized by Perdew and Zunger)")')
       elseif(author == 'pbe' .or. author == 'PBE' ) then
         write(6,*)
         write(6,'("  The potential was calculated in the generalized",  &
     &     " gradient aproximation as parametrized by Perdew, Burke ",   &
     &     "and Ernzerhof")')
       else
         write(6,*)
         write(6,'("  The XC flag is:   ",a3)') author
       endif
       if(flgscf(5:6) /= 'PW') then
         write(6,*)
         write(6,'("  The SCF flag is:   ",a6)') flgscf
         write(6,*)
         write(6,*) '  The calculation may have used a restricted basis'
       endif
       if(flgdal /= '    ') then
         if(flgdal == 'DUAL') then
           write(6,*)
           write(6,*) '  The full dual approximation was used'        
         else
           write(6,*)
           write(6,'("  The partial dual approximation was used ",a4)')  &
     &           flgdal
         endif
       endif

       read(io,iostat = ioerr) emax,teleck,nx,ny,nz,sx,sy,sz,nband,alatt
       
       if(ioerr /= 0) then
         nband = 0
         alatt = -1.0
         backspace(io)
         read(io) emax,teleck,nx,ny,nz,sx,sy,sz
       endif


       write(6,*)
       write(6,'("  The energy cutoff for wave-function kinetic ",       &
     &       "energy was ",f12.3," Hartree")') emax
       if(teleck > 0.1) then
         write(6,*)
         write(6,'("  The electronic temperature was ",f12.3," K")')     &
     &        teleck
       endif
       write(6,*)
       write(6,'("  The Brillouin zone integrations were done on a ",    &
     &     "Monkhorst-Pack grid of",i5," x ",i5," x ",i5,                &
     &     "    points")') nx,ny,nz
       write(6,'("  with a shift of",3f10.3)') sx,sy,sz

       read(io,iostat = ioerr) line250,meta_cpw2000
       
       if(ioerr == 0) then
         pwline = line250(1:60)
         title = line250(61:110)
         subtitle = line250(111:250)
         write(6,*)
         write(6,'("  The metadata of the crystal structure file ",      &
     &          "was:")')
         write(6,*) pwline,title,subtitle
         write(6,'("  The metadata of the cpw calculation was:")')
         write(6,*) meta_cpw2000
         write(6,*)
       else
         do i = 1,250
           meta_cpw2000(i:i) = ' '
         enddo
         do i = 1,140
           subtitle(i:i) = ' '
         enddo
         do i = 1,50
           title(i:i) = ' '
         enddo
         do i = 1,60
           pwline(i:i) = ' '
         enddo
         backspace(io)
         read(io,iostat = ioerr2) line140
         if(ioerr2 == 0) then
           pwline = line140(1:60)
           title = line140(61:110)
           subtitle = line140(111:140)
           write(6,*)
           write(6,'("  The metadata of the crystal structure file ",    &
     &          "was:")')
           write(6,*) pwline,title,subtitle
           write(6,*)
         else
           backspace(io)
           read(io) pwline
           write(6,*)
           write(6,'("  The old identifier of the input file was:")')
           write(6,*) pwline
         endif
       endif

!      old files do not have spacegroup information

       if(ntrans == 0) then
         read(io) ((adot(i,j),i=1,3),j=1,3)
       else
         read(io) ( (adot(i,j),i=1,3),j=1,3 ),                           &
     &            ( ((mtrx(i,j,k),j=1,3),i=1,3), k=1,ntrans ),           &
     &            ( (tnp(i,k),i=1,3), k=1,ntrans )

       endif

       read(io) (nameat(i),i=1,ntype)

       do i=1,ntype
         read(io) ((rat(j,k,i),j=1,3),k=1,natom(i))
       enddo
       
       return
       end subroutine pw_rho_v_in_crystal_calc
