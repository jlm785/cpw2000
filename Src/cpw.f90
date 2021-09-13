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


!>     Total Energy Plane Wave Calculation
!>
!>  \author       Jose Luis Martins and many others
!>  \version      5.02
!>  \date         13 September 2021
!>  \copyright    GNU Public License v2

program cpw2000


  use cpw_variables

  implicit none


  type(dims_t)                       ::  dims_                           !<  array dimensions

  type(flags_t)                      ::  flags_                          !<  computational flags

! SHOULD BE CHECKED INSIDE CHOICE_OF_PARAMS

  type(acc_t)                        ::  acc_                            !<  accuracy parameters

  type(crys_t)                       ::  crys_                           !<  crystal structure

  type(pwexp_t)                      ::  pwexp_                          !<  plane-wave expansion choices

  type(moldyn_t)                     ::  moldyn_                         !<  molecular dynamics variables

  type(vcsdyn_t)                     ::  vcsdyn_                         !<  variational cell shape molecular dynamics variables

  type(spaceg_t)                     ::  spaceg_                         !<  space group information

  type(recip_t)                      ::  recip_                          !<  reciprocal space information

  type(strfac_t)                     ::  strfac_                         !<  structure factors

  type(chdens_t)                     ::  chdens_                         !<  charge densities

  type(vcomp_t)                      ::  vcomp_                          !<  Componemts of local potential

  type(pseudo_t)                     ::  pseudo_                         !<  pseudo-potential (Kleinman-Bylander)

  type(atorb_t)                      ::  atorb_                          !<  atomic orbitals in G-space

  type(enfrst_t)                     ::  total_                          !<  Total energy force stress

  type(enfrst_t)                     ::  ewald_                          !<  Ewald energy force stress

  type(kpoint_t)                     ::  kpoint_                         !<  k-point data

  type(xc_t)                         ::  xc_                             !<  exchange and correlation choice

  type(hamallk_t)                    ::  hamallk_                        !<  hamiltonian size and indexation for all k-points

  type(psiallk_t)                    ::  psiallk_                        !<  psi for all k-points

  character(len=4)                   ::  vdriv                           !<  version of this program
  character(len=250)                 ::  meta_cpw2000                    !<  metadata from cpw.in
  character(len=250)                 ::  meta_pwdat                      !<  metadata from PW.DAT
  integer                            ::  iprglob                         !<  level of detail of printout

  real(REAL64)                       ::  ealpha                          !<  G=0 contribution to the total energy (Hartree)

  integer                            ::  icmax                           !<  maximum value of outer iteration
  logical                            ::  lsafescf                        !<  is true if self-consistency was reached.

  real(REAL64)                       ::  exc                             !<  Exchange and correlation energy (Hartree)
  real(REAL64)                       ::  rhovxc                          !<  Integral of rho times vxc (Hartree)
  real(REAL64)                       ::  strxc(3,3)                      !<  d exc / d adot,  exchange and correlatio contribution to stress tensor (contravariant components)

! potential in the FFT mesh

  integer                            ::  kmscr(7)                        !<  max value of kgv(i,n) used for the potential FFT mesh (DUAL APPROXIMATION TYPE)

! dimensions, eigenvalues, eigenvectors, occupations for all k-points

  real(REAL64)                       ::  efermi                          !<  eigenvalue of highest occupied state (T=0) or fermi energy (T/=0), Hartree
  real(REAL64)                       ::  elects                          !<  TS, electronic temperature*entropy (Hartree)


! saving density

  complex(REAL64), allocatable        ::  chdsave(:,:,:)                 !<   saved density for charge extrapolation
  integer                             ::  nsave(3)                       !<   array dimensionsfor charge extrapolation


  integer, parameter                  ::  mxdlbf = 30                    !<  array dimension for l-bfgs optimization

  character(len=2)                    ::  icorr                          !<  eXchange-Correlation choice

  integer              ::  iguess
  real(REAL64)         ::  t0,tin,tout,tinit,tfinal
  real(REAL64)         ::  telin,telout,tel0,telfin
  real(REAL64)         ::  symtol
  logical              ::  symkip
  logical                             ::  newcalc                        !  indicates that it is a new calculation (equivalent iter = ist0+1)
  logical                             ::  newcel                         !  indicates that adot has changed on output
  integer              ::  ist0,istep,iconv

  logical              ::  lkpg                                          !  If true use the previous G-vectors (same mtxd and isort)

  logical              ::  lnewnrk

  real(REAL64)          ::  deltentpy                                    !<  enthalpy difference from last iteration
  real(REAL64)          ::  errfrc                                       !<  maximum error in force (cartesian coordiantes) Hartree/Bohr
  real(REAL64)          ::  errstr                                       !<  maximum error in stress



! Driver program version

  vdriv = '5.02'

! timing

  call zesec(tinit)
  call zeelap(telin)

! prints top of output file

  call tpage(vdriv)



! reads the input file and does the respective allocations

  call cpw_read_data('cpw.in', 'PW.DAT', 5, vdriv,                       &
     iprglob, meta_cpw2000, meta_pwdat, symkip, symtol,                  &
     flags_, crys_, pwexp_, kpoint_, xc_, acc_, moldyn_, vcsdyn_,        &
     dims_, total_, ewald_)




! reads the pseudopotential data

  icorr = xc_%author

  call cpw_read_pseudo(iprglob, icorr,                                   &
       crys_, pseudo_, atorb_, dims_)





! starts or restarts the molecular dynamics, prints crystal structure
! and molecular dynamics parameters and finds symmetry operations

  ist0 = 0

  call cpw_init_mov_print_sym(ist0, iprglob, symkip, symtol,             &
       crys_, flags_, dims_, spaceg_, pwexp_, acc_, xc_,                 &
       moldyn_, vcsdyn_)



! constructs reciprocal space and allocates recip_,strfac_,pseudo_,chdens_,vcomp_


  call cpw_gspace(iprglob, kmscr,                                        &
     dims_, crys_, spaceg_, pwexp_, recip_, strfac_, pseudo_,chdens_,    &
     vcomp_, flags_)


! allow some extra space for degeneracy

  dims_%mxdbnd = pwexp_%nbandin + 5

! calculates integration k-points

  call cpw_bzint(iprglob, lnewnrk,                                       &
     dims_, kpoint_, crys_, spaceg_, pwexp_)

!
! hamiltonian matrix dimensions

  call cpw_size_alloc_hampsi(lnewnrk,                                    &
     dims_, crys_, kpoint_, pwexp_, recip_, atorb_, hamallk_, psiallk_)


! array for saving charge density (unsymmetrized)

  nsave(1) = recip_%kmax(1)
  nsave(2) = recip_%kmax(2)
  nsave(3) = recip_%kmax(3)
  allocate(chdsave(-nsave(1):nsave(1),-nsave(2):nsave(2),           &
                         -nsave(3):nsave(3)))



! Molecular Dynamics / Minimization loop

  newcel = .FALSE.
  newcalc = .TRUE.

  lkpg = .FALSE.

  do istep =  ist0+1,moldyn_%nstep

    call zesec(t0)
    call zeelap(tel0)


!-----------------------clr----------------------------------------

    if(xc_%author == "TBL") then
      if((flags_%flgcal ==  'ONE   ') .or. (flags_%flgcal == 'ONEVRD')) then
        write(6,*)
        write(6,*) 'Using TBL exchange and correlation'
        write(6,*)
      else
        write(6,*)
        write(6,'("do not use TBL to move atoms")')

        stop

      endif
    endif

!-----------------------clr----------------------------------------


!   recalculates, deallocates and reallocates for new cell

    if(newcel) then

     call cpw_newcell(lkpg, symkip, iprglob, symtol, kmscr, lnewnrk,     &
     dims_, recip_, crys_, spaceg_, pwexp_, strfac_, pseudo_,            &
     chdens_, vcomp_, flags_, kpoint_, atorb_, hamallk_, psiallk_)


    endif

!   End of newcel

    call cpw_scf_prepare(ealpha, iprglob, newcalc,                       &
       nsave, chdsave,                                                   &
       exc, strxc, rhovxc,                                               &
       dims_, crys_, recip_, strfac_, pseudo_, chdens_, vcomp_, flags_,  &
       ewald_, xc_)

    deallocate(chdsave)

    call zesec(tin)
    write(6,*)
    write(6,'("  Computing time for starting (s):  ",3x,f10.2)')         &
               tin - t0


!   Self-Consistency


    if(newcalc .or. newcel) then
      iguess = 0
    else
      iguess = 1
    endif

    icmax = 40


    if((flags_%flgscf == 'AO    ' .or.                                   &
        flags_%flgscf == 'AOJC  ' .or.                                   &
        flags_%flgscf == 'AOJCPW') .and. atorb_%latorb) then


      call cpw_scf('AO', iprglob, icmax, iguess, kmscr,                  &
      efermi, elects, exc, strxc, ealpha, lkpg, lsafescf,                &
      dims_, crys_, flags_, pwexp_, recip_, acc_, xc_, strfac_,          &
      vcomp_, pseudo_, atorb_, kpoint_, hamallk_, psiallk_,              &
      total_, ewald_, chdens_)


      iguess = 1

    endif

    if(flags_%flgscf == '    PW' .or. flags_%flgscf == 'AOJCPW') then

      call cpw_scf('PW', iprglob, icmax, iguess, kmscr,                  &
      efermi, elects, exc, strxc, ealpha, lkpg, lsafescf,                &
      dims_, crys_, flags_, pwexp_, recip_, acc_, xc_, strfac_,          &
      vcomp_, pseudo_, atorb_, kpoint_, hamallk_, psiallk_,              &
      total_, ewald_, chdens_)


    endif



!   Forces and Stresses


    call cpw_force(iprglob,strxc, ealpha, deltentpy, errfrc,             &
    errstr,                                                              &
    dims_, crys_, total_, ewald_, flags_, spaceg_, recip_, pseudo_,      &
    vcomp_, chdens_, hamallk_, psiallk_, kpoint_, vcsdyn_)


    if(flags_%flgcal /= 'ONE   ') then

      nsave(1) = recip_%kmax(1)
      nsave(2) = recip_%kmax(2)
      nsave(3) = recip_%kmax(3)
      allocate(chdsave(-nsave(1):nsave(1),-nsave(2):nsave(2),            &
                         -nsave(3):nsave(3)))

    endif

    if(flags_%flgcal == 'ONE   ') then
      iconv = 1
    else

      call cpw_move(newcel, newcalc, iconv, istep ,mxdlbf,               &
      nsave, chdsave,                                                    &
      dims_, crys_, flags_, vcsdyn_, total_, moldyn_, spaceg_,           &
      chdens_, recip_)

      newcalc = .FALSE.

    endif

    if(flags_%flgcal /= 'ONE   ') then

      if(flags_%flgcal == 'EPILBF' .or. flags_%flgcal == 'VCSLBF') then
        write(6,*)
        write(6,'(5x,2f14.6,5x,"errfrc,errstr")') errfrc,errstr
        write(6,*)
      endif

      if(flags_%flgcal == 'LBFSYM') then
        write(6,*)
        write(6,'(5x,f14.6,5x,"errfrc")') errfrc
        write(6,*)
      endif

      if((flags_%flgcal == 'EPILBF' .or. flags_%flgcal == 'VCSLBF')      &
            .and. .not. lkpg) then
        if(pwexp_%lkplusg .and. errfrc < pwexp_%epskplusg .and.          &
          errstr < pwexp_%epskplusg) then
          lkpg = .TRUE.
          write(6,*)
          write(6,*)   '  Switched on the fixed k+G mode'
          write(6,*)
        endif
      endif



    endif

    call cpw_print_crystal_energy(iprglob, iconv, istep,                 &
     elects,                                                             &
     dims_, crys_, flags_, total_, moldyn_, vcsdyn_)



    if(flags_%flgcal /= 'ONE   ') then

      call zesec(tfinal)
      call zeelap(telfin)


      write(6,*)
      write(6,'("  Computing time for md step (s):",f10.2,               &
      "    elapsed time (s):",3x,f10.2)') tfinal-t0,telfin-tel0
      write(6,*)

    endif

    if(iconv == 1) exit

  enddo


! End of molecular dynamics loop

  if(lsafescf) then

    call cpw_finish('PW_RHO_V.DAT', 21, meta_pwdat, meta_cpw2000,        &
      dims_, crys_, spaceg_, xc_, flags_, pwexp_, kpoint_, recip_,       &
      vcomp_, chdens_)

  endif



  call zesec(tout)
  call zeelap(telout)

  write(6,*)
  write(6,'("  Total computing time (s):   ",f10.2,                      &
     "    Elapsed time (s):",3x,f10.2)') tout-tinit,telout-telin
  write(6,*)


  stop

end program cpw2000
