!>    Keyword file for ESDF.  Original Author and copyright: Chris J. Pickard (c)

!     ===============================================================
!     Copyright (C) 2011 CLR & JLM, INESC-MN
!     Keywords for the cpw input package. (experimental)
!     Adapted from the parsec (1.1) package esdf_key_mod.f
!     ===============================================================

!     Modified for f90. Jose Luis Martins, December 2016.
!     Keywords added in June 2017.
!     Keywords added in May 2025.
!     Keyword MD.PWCutoffSmearing added 15 October 2026. JLM+claude
!     Keyword ChoiceOfDefaultAtomicBasisSet added 18 October 2026. JLM+claude
!     Keyword Read.PW.DAT removed (old PW.DAT format), gap at 29 closed, renumbered. 18 October 2026. JLM+claude

! Module to hold keyword list. This must be updated as
! new keywords are brought into existence.
!
! The 'label' is the label as used in calling the esdf routines
! 'typ' defines the type, with the following syntax. It is 3 characters
! long.
! The first indicates:
!  I - integer
!  S - single
!  D - double
!  P - physical
!  T - string (text)
!  E - defined (exists)
!  L - boolean (logical)
!  B - block
! The second is always a colon (:)
! The third indicates the "level" of the keyword
!  B - Basic
!  I - Intermediate
!  E - Expert
!  D - Dummy
!
! 'dscrpt' is a description of the variable. It should contain a (short) title
! enclosed between *! ... !*, and then a more detailed description of the
! variable.

Module esdf_key

  Implicit None

  Type kw_type
     Character(50)   :: label
     Character(3)    :: typ
     Character(3000) :: dscrpt
  End Type kw_type

  Integer, Parameter :: numkw=63
  Type(kw_type) :: kw(numkw)

! Now define the keywords


!     FLGCAL     TYPE OF CALCULATION (MD,...)
      data kw(1)%label     /'MD.TypeOfRun'/
      data kw(1)%typ       /'T:B' /
      data kw(1)%dscrpt                                                  &
     & /'FLGCAL (ONE,MICRO,LANG,LBFSYM,VCSLNG,VCSMIC,EPILBF,VCSLBF,      &
     &  EPILNG,RSTRT)'/

!     FLGPSD     TYPE OF PSEUDOPOTENTIAL
      data kw(2)%label     /'TypeOfPseudopotential'/
      data kw(2)%typ       /'T:B' /
      data kw(2)%dscrpt                                                  &
     & /'FLGPSD (PSEUKB) ONLY ONE OPTION HAS BEEN THOROUGHLY TESTED'/

!      FLGSCF     TYPE OF SELF CONSISTENT FIELD AND DIAGONALIZATION
      data kw(3)%label     /'TypeOfScfDiag'/
      data kw(3)%typ       /'T:B' /
      data kw(3)%dscrpt    /'FLGSCF (PW,AO ,AOJC,AOJCPW)'/

!      This should be changed to a boolean !
!      FLGDAL     WHETHER THE DUAL APPROXIMATION IS USED
      data kw(4)%label     /'DualApproximation'/
      data kw(4)%typ       /'L:B' /
      data kw(4)%dscrpt    /'FLGDAL (true, false)'/

!      FLGMIX     CHOICE OF POTENTIAL MIXING
      data kw(5)%label     /'TypeOfPotentialMixing'/
      data kw(5)%typ       /'T:B' /
      data kw(5)%dscrpt    /'FLGMIX (BFGS,BROYD1)'/

!      FLGRHO     WRITES FINAL CHARGE DENSITY TO A FILE
      data kw(6)%label     /'WriteChargeDensity'/
      data kw(6)%typ       /'L:B' /
      data kw(6)%dscrpt    /'FLGRHO (true,false)'/

!      FLGPSI     WRITES FINAL WAVEFUNCTIONS TO A FILE
      data kw(7)%label     /'WriteWaveFunctions'/
      data kw(7)%typ       /'L:B' /
      data kw(7)%dscrpt    /'FLGPSI (true,false)'/

!      FLGDOS     WRITES DATA FOR DENSITY OF STATES PLOT TO OUTPUT
      data kw(8)%label     /'WriteDOS'/
      data kw(8)%typ       /'L:B' /
      data kw(8)%dscrpt    /'FLGDOS (true,false)'/

!      FLGBND     WRITES DATA FOR BAND STRUCTURE PLOT TO OUTPUT
      data kw(9)%label     /'WriteBands'/
      data kw(9)%typ       /'L:B' /
      data kw(9)%dscrpt    /'FLGBND (true,false)'/

!      SYMKIP     WHETHER SYMMETRY SHOULD BE CONSERVED
      data kw(10)%label     /'UseSymmetry'/
      data kw(10)%typ       /'L:B' /
      data kw(10)%dscrpt    /'SYMKIP (true,false)'/

!      SYMTOL     TOLERANCE FOR SYMMETRY RECOGNITION SUBROUTINES
      data kw(11)%label     /'SymmTolerance'/
      data kw(11)%typ       /'D:B' /
      data kw(11)%dscrpt    /'SYMTOL (1.0D-5)'/

!      AUTHOR     TYPE OF XC WANTED (CA=PZ , PW92 , PBE)
      data kw(12)%label     /'XC.Authors'/
      data kw(12)%typ       /'T:B' /
      data kw(12)%dscrpt    /'AUTHOR (CA=PZ , PW92 , PBE)'/

!      IPRGLOB    AMOUNT OF PRINTING DESIRED
      data kw(13)%label     /'PrintingLevel'/
      data kw(13)%typ       /'I:B' /
      data kw(13)%dscrpt    /'IPRGLOB (0,1,2,3)'/

!      ITMAX      MAXIMUM NUMBER OF SELF CONSISTENCY CYCLES
      data kw(14)%label     /'MaxSCFIterations'/
      data kw(14)%typ       /'I:B' /
      data kw(14)%dscrpt    /'ITMAX (30)'/

!      EPSCV      CONVERGENCE CRITERIA FOR SELF CONSISTENCY
      data kw(15)%label     /'ScfTolerance'/
      data kw(15)%typ       /'D:B' /
      data kw(15)%dscrpt    /'EPSCV (0.00005)'/


!      EPSPSI     CONVERGENCE CRITERIA FOR DIAGONALIZATION
      data kw(16)%label     /'DiagTolerance'/
      data kw(16)%typ       /'D:B' /
      data kw(16)%dscrpt    /'EPSPSI (0.0001)'/

!      TEMPK      IONIC TEMPERATURE (IN KELVIN)
      data kw(17)%label     /'MD.TargetTemperature'/
      data kw(17)%typ       /'P:B' /
      data kw(17)%dscrpt    /'TEMPK (500.0 K)'/

!      TELECK     ELECTRONIC TEMPERATURE (IN KELVIN)
      data kw(18)%label     /'ElectronicTemperature'/
      data kw(18)%typ       /'P:B' /
      data kw(18)%dscrpt    /'TELECK (0.0 K)'/

!      TEMPINIK   INITIAL TEMPERATURE (IN KELVIN)
      data kw(19)%label     /'MD.InitialTemperature'/
      data kw(19)%typ       /'P:B' /
      data kw(19)%dscrpt    /'TEMPINIK (50.0 K)'/

!      NSTEP      NUMBER OF STEPS
      data kw(20)%label     /'MD.NumberOfSteps'/
      data kw(20)%typ       /'I:B' /
      data kw(20)%dscrpt    /'NSTEP (1000)'/

!      TSTEP      TIME STEP (IN ATOMIC UNITS)
      data kw(21)%label     /'MD.LengthTimeStep'/
      data kw(21)%typ       /'P:B' /
      data kw(21)%dscrpt                                                 &
     & /'TSTEP (100) (IN ATOMIC UNITS, 1 AU = 2.4 10**-17 S)'/

!      ISEED      INITIAL SEED FOR RANDOM NUMBERS
      data kw(22)%label     /'MD.Seed'/
      data kw(22)%typ       /'I:B' /
      data kw(22)%dscrpt    /'ISEED (876978) '/

!      PGTOL       CRITERIA FOR FORCE AND STRESS CONVERGENCE
      data kw(23)%label     /'MD.CG.Tolerance'/
      data kw(23)%typ       /'P:B' /
      data kw(23)%dscrpt                                                 &
     & /'PGTOL (1.0D-4) TOLERANCE FOR L-BFGS OPTIMIZATION'/

!      DXMAX       MAXIMUM STEP SIZE IN FORCES AND STRESSES
      data kw(24)%label     /'MD.CG.StepMax'/
      data kw(24)%typ       /'P:B' /
      data kw(24)%dscrpt                                                 &
     & /'DXMAX (0.1D0) MAXIMUM STEP FOR L-BFGS OPTIMIZATION'/

!      PRESS      EXTERNAL PRESSURE
      data kw(25)%label     /'MD.TargetPressure'/
      data kw(25)%typ       /'P:B' /
      data kw(25)%dscrpt    /'PRESS (0.0) EXTERNAL PRESSURE'/

!      STREXT     EXTERNAL STRESS
      data kw(26)%label     /'MD.TargetStress'/
      data kw(26)%typ       /'B:B' /
      data kw(26)%dscrpt                                                 &
     & /'STREXT (BLOCK 3x3=0.) EXTERNAL STRESS TENSOR COMPONENTS'/

!      CELMAS     FICTITIOUS CELL MASS
      data kw(27)%label     /'MD.CellMass'/
      data kw(27)%typ       /'D:B' /
      data kw(27)%dscrpt                                                 &
     & /'CELMAS (10.0D0) FICTITIOUS CELL MASS FOR VCS MD'/

      data kw(28)%label     /'ECut'/
      data kw(28)%typ       /'P:B' /
      data kw(28)%dscrpt                                                 &
     &   /'EMAX    PLANE WAVE KINETIC ENERGY CUTOFF (HARTREE)  '/

      data kw(29)%label     /'KPointGrid'/
      data kw(29)%typ       /'B:B' /
      data kw(29)%dscrpt    /'     '/

      data kw(30)%label     /'LatticeConstant'/
      data kw(30)%typ       /'P:B' /
      data kw(30)%dscrpt    /'LATTICE CONSTANT     '/

      data kw(31)%label     /'LatticeVectors'/
      data kw(31)%typ       /'B:B' /
      data kw(31)%dscrpt                                                 &
     &     /'LATTICE VECTORS IN UNITS OF LATTICE CONSTANT'/

      data kw(32)%label     /'NumberOfSpecies'/
      data kw(32)%typ       /'I:B' /
      data kw(32)%dscrpt    /'NUMBER OF DIFFERENT ATOMIC SPECIES'/

      data kw(33)%label     /'NumberOfAtoms'/
      data kw(33)%typ       /'I:B' /
      data kw(33)%dscrpt    /'TOTAL NUMBER OF ATOMS'/

      data kw(34)%label     /'Chemical_Species_Label'/
      data kw(34)%typ       /'B:B' /
      data kw(34)%dscrpt                                                 &
     & /'INDEX OF CHEMICAL SPECIES, ATOMIC NUMBER, CHEMICAL SYMBOL'/

      data kw(35)%label     /'AtomicCoordinatesAndAtomicSpecies'/
      data kw(35)%typ       /'B:B' /
      data kw(35)%dscrpt                                                 &
     & /'Coordinates of the atoms and index of species'/

      data kw(36)%label     /'SystemName'/
      data kw(36)%typ       /'T:B' /
      data kw(36)%dscrpt    /'     '/

      data kw(37)%label     /'MD.UseKeatingCorrections'/
      data kw(37)%typ       /'L:B' /
      data kw(37)%dscrpt    /'     '/

      data kw(38)%label     /'Xc.TBL.C'/
      data kw(38)%typ       /'D:B' /
      data kw(38)%dscrpt    /'        '/

      data kw(39)%label     /'NumberOfBands'/
      data kw(39)%typ       /'I:B' /
      data kw(39)%dscrpt    /'     '/

      data kw(40)%label     /'AtomSpecies'/
      data kw(40)%typ       /'B:B' /
      data kw(40)%dscrpt    /'     '/

      data kw(41)%label     /'PWEnergyCutoff'/
      data kw(41)%typ       /'P:B' /
      data kw(41)%dscrpt    /'Plane-wave energy cutoff     '/

      data kw(42)%label     /'NumberOfEigenStates'/
      data kw(42)%typ       /'I:B' /
      data kw(42)%dscrpt    /'Number of desired eigenstates     '/

      data kw(43)%label     /'kgrid_Monkhorst_Pack'/
      data kw(43)%typ       /'B:B' /
      data kw(43)%dscrpt                                                 &
     &       /'Grid for Brillouin Zone Gauss-Fourier integration  '/

      data kw(44)%label     /'Rede.NumberOfLatticePlanes'/
      data kw(44)%typ       /'I:I' /
      data kw(44)%dscrpt                                                 &
     & /'Number of lattice planes in the superlattice (rede.f90)'/

      data kw(45)%label     /'Rede.Version'/
      data kw(45)%typ       /'T:I' /
      data kw(45)%dscrpt    /'Version of rede.f90 that generated SL'/

      data kw(46)%label     /'Rede.Title'/
      data kw(46)%typ       /'T:I' /
      data kw(46)%dscrpt    /'Title for the calculation (rede.f90)'/

      data kw(47)%label     /'Rede.Date'/
      data kw(47)%typ       /'T:I' /
      data kw(47)%dscrpt    /'Day when rede.f90 was run'/

      data kw(48)%label     /'Rede.Time'/
      data kw(48)%typ       /'T:I' /
      data kw(48)%dscrpt    /'Time when rede.f90 was run'/

      data kw(49)%label     /'Rede.Name'/
      data kw(49)%typ       /'T:I' /
      data kw(49)%dscrpt    /'Name generated by rede.f90 from input'/

      data kw(50)%label     /'Rede.Superlattice'/
      data kw(50)%typ       /'B:I' /
      data kw(50)%dscrpt                                                 &
     &    /'Information about the SL vectors (rede.f90)'/

      data kw(51)%label     /'MD.FrictionFracInvTimeStep'/
      data kw(51)%typ       /'D:I' /
      data kw(51)%dscrpt                                                 &
     &    /'Friction coefficient as a function of time step,    beta =   &
     &    (1/ Delta t)/coeff. Recommended:  coeff = 20.0'/

      data kw(52)%label     /'SystemLabel'/
      data kw(52)%typ       /'T:B' /
      data kw(52)%dscrpt    /'Label that can be used to open files'/

!      LKPLUSG     FINISH CELL MINIMIZATION WITH FIXED k+G
      data kw(53)%label     /'MD.CG.UseFixedkplusG'/
      data kw(53)%typ       /'L:I' /
      data kw(53)%dscrpt                                                 &
     & /'LKPLUSG  (.FALSE.)  FINISH VCS/EPI L-BFGS OPT. WITH FIXED k+G'/

!      EPSKPLUSG     CRITERIA FOR SWITCHING TO FIXED k+G
      data kw(54)%label     /'MD.CG.FixedkplusGTol'/
      data kw(54)%typ       /'P:I' /
      data kw(54)%dscrpt                                                 &
     & /'EPSKPLUSG  (1.0D-2)   CRITERIA FOR SWITCHING TO FIXED k+G'/

!     FLGMOD     TYPE OF EMPIRICAL POTENTIAL MODEL
      data kw(55)%label     /'MD.PotentialModel'/
      data kw(55)%typ       /'T:B' /
      data kw(55)%dscrpt                                                  &
     & /'FLGMOD (LENJON,LJCLST,KEATNG)'/

!      ITDIAGMAX      MAXIMUM NUMBER OF SELF CONSISTENCY CYCLES
      data kw(56)%label     /'MaxDiagIterations'/
      data kw(56)%typ       /'I:B' /
      data kw(56)%dscrpt    /'ITMAX (40)'/

!     PSEUDO_PATH    PATH TO PSEUDOPOTENTIAL FILES
      data kw(57)%label     /'PathToPseudos'/
      data kw(57)%typ       /'T:I' /
      data kw(57)%dscrpt    /'PSEUDO_PATH  (empty) '/

!     PSEUDO_SUFFIX    SUFFIX OF THE PSEUDO FILES (ADDED TO CHEMICAL SYMBOL
      data kw(58)%label     /'PseudoSuffix'/
      data kw(58)%typ       /'T:I' /
      data kw(58)%dscrpt    /'PSEUDO_SUFFIX   _POTKB_F.DAT'/

!     ITAPE_PSEUDO     TAPE NUMBER FOR READING PSEUDOPOTENTIAL FILES
      data kw(59)%label     /'TapeToReadPseudo'/
      data kw(59)%typ       /'I:I' /
      data kw(59)%dscrpt    /'ITAPE_PSEUDO     40'/

!     ITAPE_SAVE_PSI     TAPE NUMBER FOR SAVING WAVE-FUNCTIONS TO DISK (IF > 10)
      data kw(60)%label     /'TapeToSavePsi'/
      data kw(60)%typ       /'I:I' /
      data kw(60)%dscrpt    /'ITAPE_SAVE_PSI     0'/

!     SAVE_PSI_PATH     PATH FOR EVENTUAL SAVING WAVE-FUNCTIONS TO DISK
      data kw(61)%label     /'PathToSavePsi'/
      data kw(61)%typ       /'T:I' /
      data kw(61)%dscrpt    /'SAVE_PSI_PATH  (empty) '/

!     ESMEAR     SMEARING OF THE KINETIC ENERGY CUTOFF
      data kw(62)%label     /'MD.PWCutoffSmearing'/
      data kw(62)%typ       /'P:I' /
      data kw(62)%dscrpt    /'ESMEAR (min(0.5 Ha, 5% cutoff) or 0)'/

!     IBSET      CHOICE OF ATOMIC BASIS SET FOR EACH SPECIES
      data kw(63)%label     /'ChoiceOfDefaultAtomicBasisSet'/
      data kw(63)%typ       /'B:I' /
      data kw(63)%dscrpt    /'IBSET (species index, basis set; default 1)'/


      end module esdf_key
