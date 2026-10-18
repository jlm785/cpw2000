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

!>  Reads the parameters of the calculation
!>
!>  \author       Jose Luis Martins and Carlos Loia Reis
!>  \version      5.13
!>  \date         5 August 2002. 18 October 2026.
!>  \copyright    GNU Public License v2

subroutine read_esdf(fname, vdriv, flgcal, flgkeat,                      &
    flgpsd, flgscf, flgdal, flgmix,                                      &
    adot, ntype, natom, nameat, rat, atmass, alatt,                      &
    emax, esmear, nbandin, nx, ny, nz, sx, sy, sz,                       &
    meta_cpw2000,                                                        &
    symkip, symtol,                                                      &
    author, tblaha, iprglob, itmax, icdiagmax, epscv, epscvao, epspsi,   &
    tempk, teleck, tempinik, nstep, tstep, beta, iseed, pgtol, dxmax,    &
    press, strext, celmas, flgkplusg, epskplusg,                         &
    pseudo_path, pseudo_suffix, itape_pseudo,                            &
    save_psi_path, itape_save_psi,                                       &
    file_restart, iotape_restart,                                        &
    jbset,                                                               &
    mxdtyp, mxdatm)

! Written August 5, 2002. jlm
! Modified 18 September 2002. jlm
! Modified 6 February 2008
! Modified to use esdf as interface. CLR
! Modified, f90, comments, 6 October 2015. JLM
! Modified to read geometry and get rid of PW.DAT, June 2017. JLM
! modified (iargc) 27 June 2017.  JLM
! Modified, documentation, kplusg,August 10 2019. JLM
! Modified, icdiagmax, indentation. 14 May 2025. JLM
! Modified, filenames for pseudo and saving psi to disk. 10 October 2025. JLM
! Documentation, one argument per declaration. 28 September 2026. JLM+claude
! esmear, smearing of the kinetic energy cutoff (MD.PWCutoffSmearing). 15 October 2026. JLM+claude WARNING NEW API
! RSTRT default of esmear from the restart file, file_restart, iotape_restart, EPS. 16 October 2026. JLM+claude WARNING NEW API
! Old PW.DAT no longer accepted, lgeom, lbz removed, jbset (read_esdf_atorb_choice). 18 October 2026. JLM+claude WARNING NEW API


  use esdf

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of types of atoms

  character(len=*), intent(in)       ::  fname                           !<  filename with input data
  character(len=4), intent(in)       ::  vdriv                           !<  version of the calling program

  character(len=*), intent(in)       ::  file_restart                    !<  name of the file with the restart information (to find the type of run if flgcal = RSTRT)
  integer, intent(in)                ::  iotape_restart                  !<  tape number of the restart file

! output

  character(len=6), intent(out)      ::  flgcal                          !<  type of calculation
  character(len=6), intent(out)      ::  flgkeat                         !<  adds keating force field correction to energy

  character(len=6), intent(out)      ::  flgpsd                          !<  type of pseudopotential
  character(len=6), intent(out)      ::  flgscf                          !<  type of self consistent field and diagonalization
  character(len=4), intent(out)      ::  flgdal                          !<  whether the dual approximation is used
  character(len=6), intent(out)      ::  flgmix                          !<  choice of potential mixing

  real(REAL64), intent(out)          ::  adot(3,3)                       !<  metric in direct space
  integer, intent(out)               ::  ntype                           !<  number of types of atoms
  integer, intent(out)               ::  natom(mxdtyp)                   !<  number of atoms of type i
  character(len=2), intent(out)      ::  nameat(mxdtyp)                  !<  chemical symbol for the type i
  real(REAL64), intent(out)          ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i
  real(REAL64), intent(out)          ::  atmass(mxdtyp)                  !<  atomic mass (in a.u.) of atom of type i

  real(REAL64), intent(out)          ::  alatt                           !<  lattice constant

  real(REAL64), intent(out)          ::  emax                            !<  kinetic energy cutoff of plane wave expansion (Hartree).
  real(REAL64), intent(out)          ::  esmear                          !<  smearing of the kinetic energy cutoff (Hartree), no smearing if zero

  integer, intent(out)               ::  nx                              !<  size of the integration mesh in k-space (nx*ny*nz), direction 1
  integer, intent(out)               ::  ny                              !<  size of the integration mesh in k-space (nx*ny*nz), direction 2
  integer, intent(out)               ::  nz                              !<  size of the integration mesh in k-space (nx*ny*nz), direction 3
  real(REAL64), intent(out)          ::  sx                              !<  offset of the integration mesh (usually 0.5), direction 1
  real(REAL64), intent(out)          ::  sy                              !<  offset of the integration mesh (usually 0.5), direction 2
  real(REAL64), intent(out)          ::  sz                              !<  offset of the integration mesh (usually 0.5), direction 3

  integer, intent(out)               ::  nbandin                         !<  target for number of bands

  character(len=250), intent(out)    ::  meta_cpw2000                    !<  metadata coded in cpw2000

  logical, intent(out)               ::  symkip                          !<  whether symmetry should be conserved
  real(REAL64), intent(out)          ::  symtol                          !<  tolerance for symmetry recognition subroutines

  character(len=4), intent(out)      ::  author                          !<  type of xc wanted (ca=pz , pw92 , pbe,...)
  real(REAL64), intent(out)          ::  tblaha                          !<  Tran-Blaha constant

  integer, intent(out)               ::  iprglob                         !<  controls the amount of printing by subroutines
  integer, intent(out)               ::  itmax                           !<  maximum number of self consistency cycles
  integer, intent(out)               ::  icdiagmax                       !<  maximum number of diagonalization  outer iterations
  real(REAL64), intent(out)          ::  epscv                           !<  convergence criteria for potential self consistency
  real(REAL64), intent(out)          ::  epscvao                         !<  convergence criteria for potential self consistency in atomic orbitals
  real(REAL64), intent(out)          ::  epspsi                          !<  convergence criteria for iterative diagonalization

  real(REAL64), intent(out)          ::  tempk                           !<  ionic temperature (in Kelvin)
  real(REAL64), intent(out)          ::  teleck                          !<  electronic temperature (in Kelvin)
  real(REAL64), intent(out)          ::  tempinik                        !<  initial ionic temperature (in Kelvin)
  integer, intent(out)               ::  nstep                           !<  number of MD or minimization steps
  real(REAL64), intent(out)          ::  tstep                           !<  time step for MD (in atomic units)
  real(REAL64), intent(out)          ::  beta                            !<  friction coefficient/mass (in a.u.)
  integer, intent(out)               ::  iseed                           !<  initial seed for random number generator
  real(REAL64), intent(out)          ::  pgtol                           !<  criteria for force and stress convergence
  real(REAL64), intent(out)          ::  dxmax                           !<  maximum step size in forces and stresses

  real(REAL64), intent(out)          ::  press                           !<  external pressure
  real(REAL64), intent(out)          ::  strext(3,3)                     !<  external applied stress
  real(real64), intent(out)          ::  celmas                          !<  fictitious cell mass

  logical, intent(out)               ::  flgkplusg                       !<  finish cell minimization with fixed k+G
  real(REAL64), intent(out)          ::  epskplusg                       !<  criteria for switching to fixed k+G

  character(len=200), intent(out)    ::  pseudo_path                     !<  path to pseudopotentials
  character(len=50), intent(out)     ::  pseudo_suffix                   !<  suffix for the pseudopotentials
  integer, intent(out)               ::  itape_pseudo                    !<  tape number to read pseudo

  character(len=200), intent(out)    ::  save_psi_path                   !<  path to save psi files to disk
  integer, intent(out)               ::  itape_save_psi                  !<  tape number to read and write psi.  If < 10 (default 0) do not use.

  integer, intent(out)               ::  jbset(mxdtyp)                   !<  atomic basis set chosen for atom type nt

! local variables

  logical                    ::  lgeom                                   !  indicates if geometry was successfully read.
  logical                    ::  lbz                                     !  indicates if Brillouin Zone data was successfully read.

  character(len=4)           ::  vdrloc
  logical                    ::  ldevel                                  !  development branch, minor version may be incompatible
  integer                    ::  ipr
  logical                    ::  lrede

  character(len=6)           ::  flgold                                  !  type of calculation in the restart file
  integer                    ::  ios

! parameters

  real(REAL64), parameter ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter ::  EPS = 1.0E-7_REAL64                       !  smaller values of esmear are set to zero


! checks consistency between cpw2000 and this subroutine


  call version(vdrloc,ldevel)

  if(ldevel) then

    if(vdriv  /=  vdrloc) then
      write(6,*)
      write(6,'("   Stopped in read_esdf, version ",a4,                  &
          & " because it was called from a cpw version ",a4)')           &
            vdrloc,vdriv

      stop

    endif

  else

    if(vdriv(1:3)  /=  vdrloc(1:3)) then
      write(6,*)
      write(6,'("   Stopped in read_esdf, version ",a4,                  &
          & " because it was called from a cpw version ",a4)')           &
            vdrloc,vdriv

      stop

    endif

    if(vdriv(4:4) /= vdrloc(4:4)) then

      write(6,*)
      write(6,'("     WARNING    WARNING     in read_esdf")')
      write(6,*)
      write(6,'("  library version: ",a4," cpw driver version: ",a4)') vdriv,vdrloc
      write(6,*)

    endif

  endif

! esdf initialization

  write(6,*)
  write(6,*) '  input read from ',fname
  write(6,*)

  call esdf_init(fname)

! type of calculation


! ONE:    just do one scf calculation
! MICRO:  do a microcanonical molecular dynamical simulation
! LANG:   do a langevin molecular dynamical simulation
! LBFSYM: optimization of internal coordinates with symmetry conservation
! VCSLNG: variable cell shape with langevin
! VCSMIC: variable cell shape with microcanonic
! VCSLBF: variable cell shape optimization
! EPILNG: variable cell shape with langevin and epitaxial constraint
! RSTRT:  restarts an old calculation from the point where it was interrupted

  flgcal = 'ONE   '
! flgcal = 'MICRO '
! flgcal = 'LANG  '
! flgcal = 'LBFSYM'
! flgcal = 'VCSLNG'
! flgcal = 'VCSMIC'
! flgcal = 'VCSLBF'
! flgcal = 'EPILBF'
! flgcal = 'EPILNG'
! flgcal = 'RSTRT '

  flgcal = esdf_string('MD.TypeOfRun',flgcal)
  call chrcap(flgcal,6)

  if (command_argument_count() >=1) then
!  if (iargc() >=1) then
!    call getarg(1,flgcal)
    call get_command_argument(1,flgcal)
  endif

  if(flgcal == 'ONE   ') nstep = 1


! choice of exchange and correlation

  author = 'CA  '
! author = 'PBE '

  author = esdf_string('XC.Authors', author)

  if (command_argument_count() >=2) then
!  if (iargc() >=2) then
!    call getarg(2,author)
    call get_command_argument(2,author)
  endif
  call chrcap(author,4)


! and Tran-Blaha constant

  tblaha = UM
  tblaha = esdf_double('Xc.TBL.C', tblaha)


! global printing flag, the highest the more detail is printed


  if(flgcal == 'ONE   ') then
    iprglob = 3
  else
    iprglob = 1
  endif

  iprglob = esdf_integer('PrintingLevel', iprglob)


! type of pseudopotential
! currently only one option is available


  flgpsd = 'PSEUKB'

  flgpsd = esdf_string('TypeOfPseudopotential', flgpsd)
  call chrcap(flgpsd,6)


  ipr = iprglob

  call read_esdf_crystal(ipr,                                            &
      adot, ntype, natom, nameat, rat,                                   &
      atmass, alatt, lgeom,                                              &
      mxdtyp, mxdatm)

  if(.NOT. lgeom) then
    write(6,*)
    write(6,'("   STOPPED in read_esdf:  the crystal structure was not ",&
       &      "found in ",a)') trim(fname)
    write(6,'("   The old PW.DAT format is no longer accepted, ",        &
       &      "convert it with Tools/pwdat2cpwin.f90")')
    write(6,*)

    stop

  endif

  call read_esdf_bz(ipr,                                                 &
      emax, nbandin, nx, ny, nz, sx, sy, sz, lbz)

! smearing of the kinetic energy cutoff (see kinetic_smear).  By default
! only for runs with variable cell, with the value recommended by ABINIT
! for ecutsm (0.5 Hartree), but not more than 5% of the cutoff.
! A restart uses the default of the type of run in the restart file.

  flgold = flgcal
  if(flgcal == 'RSTRT ') then
    open(unit = iotape_restart, file = adjustl(trim(file_restart)),      &
         status = 'old', form = 'formatted', iostat = ios)
    if(ios == 0) read(iotape_restart,'(2x,a6)', iostat = ios) flgold
    if(ios /= 0) then
      write(6,*)
      write(6,'("   STOPPED in read_esdf:  could not read the type of ", &
         &      "calculation from the restart file ",a)')                &
                adjustl(trim(file_restart))
      write(6,*)

      stop

    endif
    close(unit = iotape_restart)
  endif

  esmear = ZERO
  if(flgold(1:3) == 'VCS' .or. flgold(1:3) == 'EPI' .or.                 &
     flgold(1:3) == 'EOS') then
    esmear = min(UM/2, emax/20)
  endif

  if(esdf_defined('MD.PWCutoffSmearing')) then
    esmear = esdf_physical('MD.PWCutoffSmearing', esmear, 'hartree')
  endif

  if(esmear >= ZERO .and. esmear < EPS) esmear = ZERO

  if(esmear < ZERO .or. esmear >= emax) then
    write(6,*)
    write(6,'("   STOPPED in read_esdf:  MD.PWCutoffSmearing = ",        &
       &      f12.4," should be positive and smaller than the cutoff")') &
              esmear
    write(6,*)

    stop

  endif

  if(ipr > 0 .and. esmear > ZERO) then
    write(6,*)
    write(6,'("   The kinetic energy cutoff is smeared by ",f10.4,       &
       &      " hartree")') esmear
    write(6,*)
  endif

  call read_esdf_rede(ipr, meta_cpw2000, lrede)

  if(.NOT. lrede) then
    meta_cpw2000 = esdf_string('SystemLabel', meta_cpw2000)
  endif

  call read_esdf_md(ipr,                                                 &
      tempk, tempinik, nstep, tstep, flgkeat,                            &
      beta, iseed, pgtol, dxmax, press, strext, celmas,                  &
      flgkplusg, epskplusg)

  call read_esdf_options(flgcal, ipr,                                    &
      flgscf, flgdal, flgmix, teleck,                                    &
      itmax, icdiagmax, epscv, epscvao, epspsi, symkip, symtol)

  call read_esdf_filenames(ipr,                                          &
      pseudo_path, pseudo_suffix, itape_pseudo,                          &
      save_psi_path, itape_save_psi)

  call read_esdf_atorb_choice(ipr, ntype, nameat, jbset,                 &
      mxdtyp)


! deallocates the esdf arrays

  call esdf_close


! prints the results

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




