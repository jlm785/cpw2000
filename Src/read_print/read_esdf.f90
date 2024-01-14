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

!>     reads the parameters of the calculation

       subroutine read_esdf(fname,vdriv,flgcal,flgkeat,                  &
     & flgpsd,flgscf,flgdal,flgmix,                                      &
     & adot,ntype,natom,nameat,rat,atmass,alatt,lgeom,                   &
     & emax,nbandin,nx,ny,nz,sx,sy,sz,lbz,                               &
     & meta_cpw2000,                                                     &
     & symkip,symtol,                                                    &
     & author,tblaha,iprglob,itmax,epscv,epscvao,epspsi,                 &
     & tempk,teleck,tempinik,nstep,tstep,beta,iseed,pgtol,dxmax,         &
     & press,strext,celmas,flgkplusg,epskplusg,                          &
     & mxdtyp,mxdatm)

!      Written August 5, 2002. jlm
!      Modified 18 September 2002. jlm
!      Modified 6 February 2008
!      Modified to use esdf as interface. CLR
!      Modified, f90, comments, 6 October 2015. JLM
!      Modified to read geometry and get rid of PW.DAT, June 2017. JLM
!      modified (iargc) 27 June 2017.  JLM
!      Modified, documentation, kplusg,August 10 2019. JLM
!      copyright inesc-mn/Jose Luis Martins/Carlos Loia Reis


!      version 4.94

       use esdf

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdtyp                     !<  array dimension of types of atoms
       integer, intent(in)                ::  mxdatm                     !<  array dimension of types of atoms

       character(len=*), intent(in)       ::  fname                      !<  filename with input data
       character(len=4), intent(in)       ::  vdriv                      !<  version of the calling program

!      output

       character(len=6), intent(out)      ::  flgcal                     !<  type of calculation
       character(len=6), intent(out)      ::  flgkeat                    !<  adds keating force field correction to energy

       character(len=6), intent(out)      ::  flgpsd                     !<  type of pseudopotential
       character(len=6), intent(out)      ::  flgscf                     !<  type of self consistent field and diagonalization
       character(len=4), intent(out)      ::  flgdal                     !<  whether the dual approximation is used
       character(len=6), intent(out)      ::  flgmix                     !<  choice of potential mixing

       real(REAL64), intent(out)          ::  adot(3,3)                  !<  metric in direct space
       integer, intent(out)               ::  ntype                      !<  number of types of atoms
       integer, intent(out)               ::  natom(mxdtyp)              !<  number of atoms of type i
       character(len=2), intent(out)      ::  nameat(mxdtyp)             !<  chemical symbol for the type i
       real(REAL64), intent(out)          ::  rat(3,mxdatm,mxdtyp)       !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i
       real(REAL64), intent(out)          ::  atmass(mxdtyp)             !<  atomic mass (in a.u.) of atom of type i

       real(REAL64), intent(out)          ::  alatt                      !<  lattice constant
       logical, intent(out)               ::  lgeom                      !<  indicates if geometry was successfully read.

       real(REAL64), intent(out)          ::  emax                       !<  kinetic energy cutoff of plane wave expansion (Hartree).

       integer, intent(out)               ::  nx, ny, nz                 !<  size of the integration mesh in k-space (nx*ny*nz)
       real(REAL64), intent(out)          ::  sx, sy, sz                 !<  offset of the integration mesh (usually 0.5)

       integer, intent(out)               ::  nbandin                    !<  target for number of bands
       logical, intent(out)               ::  lbz                        !<  indicates if Brillouin Zone data was successfully read.

       character(len=250), intent(out)    ::  meta_cpw2000               !<  metadata coded in cpw2000

       logical, intent(out)               ::  symkip                     !<  whether symmetry should be conserved
       real(REAL64), intent(out)          ::  symtol                     !<  tolerance for symmetry recognition subroutines

       character(len=4), intent(out)      ::  author                     !<  type of xc wanted (ca=pz , pw92 , pbe,...)
       real(REAL64), intent(out)          ::  tblaha                     !<  Tran-Blaha constant

       integer, intent(out)               ::  iprglob                    !<  controls the amount of printing by subroutines
       integer, intent(out)               ::  itmax                      !<  maximum number of self consistency cycles
       real(REAL64), intent(out)          ::  epscv                      !<  convergence criteria for potential self consistency
       real(REAL64), intent(out)          ::  epscvao                    !<  convergence criteria for potential self consistency in atomic orbitals
       real(REAL64), intent(out)          ::  epspsi                     !<  convergence criteria for iterative diagonalization

       real(REAL64), intent(out)          ::  tempk                      !<  ionic temperature (in Kelvin)
       real(REAL64), intent(out)          ::  teleck                     !<  electronic temperature (in Kelvin)
       real(REAL64), intent(out)          ::  tempinik                   !<  initial ionic temperature (in Kelvin)
       integer, intent(out)               ::  nstep                      !<  number of MD or minimization steps
       real(REAL64), intent(out)          ::  tstep                      !<  time step for MD (in atomic units)
       real(REAL64), intent(out)          ::  beta                       !<  friction coefficient/mass (in a.u.)
       integer, intent(out)               ::  iseed                      !<  initial seed for random number generator
       real(REAL64), intent(out)          ::  pgtol                      !<  criteria for force and stress convergence
       real(REAL64), intent(out)          ::  dxmax                      !<  maximum step size in forces and stresses

       real(REAL64), intent(out)          ::  press                      !<  external pressure
       real(REAL64), intent(out)          ::  strext(3,3)                !<  external applied stress
       real(real64), intent(out)          ::  celmas                     !<  fictitious cell mass

       logical, intent(out)               ::  flgkplusg                  !<  finish cell minimization with fixed k+G
       real(REAL64), intent(out)          ::  epskplusg                  !<  criteria for switching to fixed k+G

!      local variables

       character(len=4)           ::  vdrloc
       logical                    ::  ldevel                             !  development branch, minor version may be incompatible
       integer                    ::  ipr
       logical                    ::  lrede

!      parameters

       real(REAL64), parameter ::  UM = 1.0_REAL64


!      checks consistency between cpw2000 and this subroutine


       call version(vdrloc,ldevel)

       if(ldevel) then

         if(vdriv  /=  vdrloc) then
           write(6,*)
           write(6,'("   Stopped in read_esdf, version ",a4,             &
     &     " because it was called from a cpw version ",a4)')            &
     &         vdrloc,vdriv

           stop

         endif

       else

         if(vdriv(1:3)  /=  vdrloc(1:3)) then
           write(6,*)
           write(6,'("   Stopped in read_esdf, version ",a4,             &
     &     " because it was called from a cpw version ",a4)')            &
     &         vdrloc,vdriv

           stop

         endif

         if(vdriv(4:4) /= vdrloc(4:4)) then

           write(6,*)
           write(6,'("     WARNING    WARNING     in read_esdf")')
           write(6,*)
           write(6,'("  library version: ",a4," cpw driver version: ",   &
     &                a4)') vdriv,vdrloc
           write(6,*)

         endif

       endif

!      esdf initialization

       write(6,*)
       write(6,*) '  input read from ',fname
       write(6,*)

       call esdf_init(fname)

!      type of calculation


!      ONE:    just do one scf calculation
!      MICRO:  do a microcanonical molecular dynamical simulation
!      LANG:   do a langevin molecular dynamical simulation
!      LBFSYM: optimization of internal coordinates with symmetry conservation
!      VCSLNG: variable cell shape with langevin
!      VCSMIC: variable cell shape with microcanonic
!      VCSLBF: variable cell shape optimization
!      EPILNG: variable cell shape with langevin and epitaxial constraint
!      RSTRT:  restarts an old calculation from the point where it was interrupted

       flgcal = 'ONE   '
!      flgcal = 'MICRO '
!      flgcal = 'LANG  '
!      flgcal = 'LBFSYM'
!      flgcal = 'VCSLNG'
!      flgcal = 'VCSMIC'
!      flgcal = 'VCSLBF'
!      flgcal = 'EPILBF'
!      flgcal = 'EPILNG'
!      flgcal = 'RSTRT '

       flgcal = esdf_string('MD.TypeOfRun',flgcal)
       call chrcap(flgcal,6)

       if (command_argument_count() >=1) then
!       if (iargc() >=1) then
!         call getarg(1,flgcal)
         call get_command_argument(1,flgcal)
       endif

       if(flgcal == 'ONE   ') nstep = 1


!      choice of exchange and correlation


       author = 'CA  '
!      author = 'PBE '

       author = esdf_string('XC.Authors',author)
       call chrcap(author,4)

       if (command_argument_count() >=2) then
!       if (iargc() >=2) then
!         call getarg(2,author)
         call get_command_argument(2,author)
       endif

!      and Tran-Blaha constant

       tblaha = UM
       tblaha = esdf_double('Xc.TBL.C',tblaha)


!      global printing flag, the highest the more detail is printed


       if(flgcal == 'ONE   ') then
         iprglob = 3
       else
         iprglob = 1
       endif

       iprglob = esdf_integer('PrintingLevel',iprglob)


!      type of pseudopotential
!      currently only one option is available


       flgpsd = 'PSEUKB'

       flgpsd = esdf_string('TypeOfPseudopotential',flgpsd)
       call chrcap(flgpsd,6)


       ipr = iprglob

       call read_esdf_crystal(ipr,                                       &
     & adot,ntype,natom,nameat,rat,                                      &
     & atmass,alatt,lgeom,                                               &
     & mxdtyp,mxdatm)

       call read_esdf_bz(ipr,                                            &
     & emax,nbandin,nx,ny,nz,sx,sy,sz,lbz)

       call read_esdf_rede(ipr,meta_cpw2000,lrede)

       if(.NOT. lrede) then
         meta_cpw2000 = esdf_string('SystemLabel',meta_cpw2000)
       endif

       call read_esdf_md(ipr,                                            &
     & tempk,tempinik,nstep,tstep,flgkeat,                               &
     & beta,iseed,pgtol,dxmax,press,strext,celmas,flgkplusg,epskplusg)

       call read_esdf_options(flgcal,ipr,                                &
     & flgscf,flgdal,flgmix,teleck,                                      &
     & itmax,epscv,epscvao,epspsi,symkip,symtol)


!      deallocates the esdf arrays

       call esdf_close


!      prints the results

       if(ipr > 2) then
         write(6,*)
         write(6,*)
         write(6,*)  '  The values set by read_esdf are:'
         write(6,*)  '  (but not on called subroutines)'
         write(6,*)
         write(6,'("   The value of flgcal is: ",a6)') flgcal
         write(6,'("   The value of author is: ",a3)') author
         write(6,'("   The value of tblaha is: ",f14.6)') tblaha
         write(6,'("   The value of iprglob is: ",i5)') iprglob
         write(6,'("   The value of flgpsd is: ",a6)') flgpsd
         write(6,*)
         write(6,*)
       endif



       return
       end subroutine read_esdf




