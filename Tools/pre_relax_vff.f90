!>     pre-relaxation with valence force field

       program pre_relax_vff

!      written 15 january 1999. jlm
!      Modified, f90, new subroutines 3 June 2019. JLM
!      Specialized for VFF relaxation
!      copyright inesc-mn/jose luis martins

       implicit none
       integer, parameter          :: REAL64 = selected_real_kind(12)

!      main variables

       integer                            ::  mxdtyp                     !< array dimension of types of atoms
       integer                            ::  mxdatm                     !< array dimension of types of atoms

       real(REAL64)                       ::  adot(3,3)                  !< metric in direct space
       integer                            ::  ntype                      !< number of types of atoms
       integer, allocatable               ::  natom(:)                   !< number of atoms of type i
       character(len=2), allocatable      ::  nameat(:)                  !< chemical symbol for the type i
       real(REAL64), allocatable          ::  rat(:,:,:)                 !< k-th component (in lattice coordinates) of the position of the n-th atom of type i
       real(REAL64), allocatable          ::  atmass(:)                  !< atomic mass (in a.u.) of atom of type i 

       real(REAL64)                       ::  alatt                      !< lattice constant

       real(REAL64), allocatable          ::  force(:,:,:)               !< d energy / d rat,  force on the n-th atom of type i (contravarian components, hartree/bohr)
       real(REAL64)                       ::  stress(3,3)                !< d energy / d adot,  stress tensor (contravariant)
       real(REAL64)                       ::  energy                     !< energy

       real(REAL64), allocatable          ::  vat(:,:,:)                 !< d rat / d t  velocity in lattice coordinates of atom j of type i
       real(REAL64)                       ::  vadot(3,3)                 !< d adot / d t  rate of change of metric

       logical                            ::  newcel                     !< indicates that adot has changed on output

       real(REAL64), allocatable          ::  rat1(:,:,:)                !< previous value of lattice coordinates of atom j of type i
       real(REAL64), allocatable          ::  frc1(:,:,:)                !< previous value of force on atom j of type i
       real(REAL64)                       ::  adot1(3,3)                 !< previous metric
       real(REAL64)                       ::  frcel1(3,3)                !< previous fictitious cell force

       real(REAL64)                       ::  ekin                       !< kinetic energy of the nuclei (atomic units)

!       real(REAL64)                       ::  acc
       integer                            ::  ipr
       character(len=50)                  ::  title
       integer                            ::  io
       character(len=20)                  ::  filename
!       character(len=20)                  ::  filepw

       logical                            ::  lgeom

!      md variables

       real(REAL64)                       ::  tempk                      !<  ionic temperature (in Kelvin)
       real(REAL64)                       ::  tempinik                   !<  initial ionic temperature (in Kelvin)
       integer                            ::  nstep                      !<  number of MD or minimization steps
       real(REAL64)                       ::  tstep                      !<  time step for MD (in atomic units)
       real(REAL64)                       ::  beta                       !<  friction coefficient/mass (in a.u.)
       integer                            ::  iseed                      !<  initial seed for random number generator        
       real(REAL64)                       ::  pgtol                      !<  criteria for force and stress convergence
       real(REAL64)                       ::  dxmax                      !<  maximum step size in forces and stresses
       
       real(REAL64)                       ::  press                      !<  external pressure
       real(REAL64)                       ::  strext(3,3)                !<  external applied stress
       real(real64)                       ::  celmas                     !<  fictitious cell mass

       real(REAL64)                       ::  ekcell                     !< fictitious cell kinetic energy

       character(len=6)                   ::  flgcal                     !<  type of calculation
       character(len=6)                   ::  flgmod                     !<  classical potential model
!       logical                            ::  lcluster                   !<  if true it is a cluster (not periodic)

!      pass through variables (not used by md)

       real(REAL64)                       ::  emax                       !<  kinetic energy cutoff of plane wave expansion (Hartree).
       integer                            ::  nbandin                    !<  target for number of bands      
       integer                            ::  nx, ny, nz                 !<  size of the integration mesh in k-space (nx*ny*nz)
       real(REAL64)                       ::  sx, sy, sz                 !<  offset of the integration mesh (usually 0.5)
       character(len=250)                 ::  meta_pwdat                 !<  metadata from cpw_in or PW.DAT

!      other variables

       integer                            ::  iconv
       integer                            ::  iotape

!      constants

       real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
       integer, parameter                 ::  mxdlbf = 20                !< array dimension for l-bfgs optimization

!      counters

       integer         ::  nt1, nt2, istep

!       lcluster = .FALSE. 
!       acc = 8.0
       iconv = 0

       ipr = 3
       filename = 'cpw.in'
!       filepw = 'LJ.DAT'
       call size_mxdtyp_mxdatm_esdf(ipr,filename,lgeom,mxdtyp,mxdatm)
!        if(.NOT. lgeom) then
!          io = 20
!          call size_mxdtyp_mxdatm(filepw,io,mxdtyp,mxdatm)
!        endif

       allocate(natom(mxdtyp))
       allocate(nameat(mxdtyp))
       allocate(rat(3,mxdatm,mxdtyp))
       allocate(atmass(mxdtyp))
       allocate(force(3,mxdatm,mxdtyp))

       allocate(vat(3,mxdatm,mxdtyp))

       allocate(rat1(3,mxdatm,mxdtyp))
       allocate(frc1(3,mxdatm,mxdtyp))

!       if(lgeom) then

         call read_md_esdf(ipr,filename,flgcal,flgmod,                   &
     &   adot,ntype,natom,nameat,rat,atmass,alatt,lgeom,                 &
     &   emax,nbandin,nx,ny,nz,sx,sy,sz,                                 &
     &   meta_pwdat,                                                     &
     &   tempk,tempinik,nstep,tstep,beta,iseed,pgtol,dxmax,              &
     &   press,strext,celmas,                                            &
     &   mxdtyp,mxdatm)

!        use higher accuracy

         pgtol = 0.001*pgtol

!        else
! 
!        read data from old file
! 
!          write(6,*)
!          write(6,*) ' WARNING '
!          write(6,*) ' Using file ',filepw
!          write(6,*)
! 
!          open(unit=5,file=filepw,form='formatted',status='old')
!          read(5,'(a50)') title
! 
!          call read_data(adot,ntype,natom,nameat,rat,atmass,alatt,        &
!      &   mxdtyp,mxdatm)
! 
!        endif


!      Random initial velocities

       call move_random_v_atom(iseed,tempinik,vat,vadot,                 &
     & ntype,natom,atmass,adot,                                          &
     & mxdtyp,mxdatm)

!       newcel = .TRUE.
!       newcalc = .TRUE.

       call print_md_parameters(flgcal,                                  &
     & tempk,tempinik,nstep,tstep,beta,iseed,                            &
     & press,strext,celmas)

       do istep =  1,nstep

         ipr = 1
         call print_crystal(ipr,adot,ntype,natom,nameat,rat,             &
     &   mxdtyp,mxdatm)

!          if(flgmod == 'LENJON' .or. flgmod == 'LJCLST') then
! 
!            lcluster = .FALSE.
!            if(flgmod == 'LJCLST') lcluster = .TRUE.
! 
!            call lennard_jones(lcluster,acc,                              &
!      &     ntype,natom,nameat,rat,adot,                                  &
!      &     force,energy,stress,                                          &
!      &     mxdtyp,mxdatm)
! 
!          else

           call vff_keating(ipr, 6,                                      &
     &       ntype, natom, nameat, rat, adot,                            &
     &       force, energy, stress,                                      &
     &       mxdtyp, mxdatm)

!         endif

         ipr = 1
         call print_energy(ipr,'potential',energy,force,stress,          &
     &   adot,ntype,natom,nameat,                                        &
     &   mxdtyp,mxdatm)

         call move_md_atom_cell(flgcal,newcel,iconv,istep,               &
     &   rat,vat,adot,vadot,energy,force,stress,                         &
     &   tstep,ekin,ekcell,strext,press,celmas,                          &
     &   beta,tempk,iseed,                                               &
     &   mxdlbf,pgtol,dxmax,                                             &
     &   rat1,frc1,adot1,frcel1,                                         &
     &   ntype,natom,nameat,atmass,                                      &
     &   mxdtyp,mxdatm)

         if(adjustl(trim(flgcal)) == 'ONE') then

           exit

         endif

         call move_print_mdenergy(flgcal,istep,energy,ekin,ekcell,ZERO)

         if(iconv == 1) exit

       enddo

       if(adjustl(trim(flgcal)) /= 'ONE') then
         call print_crystal(ipr,adot,ntype,natom,nameat,rat,             &
     &   mxdtyp,mxdatm)
       endif

       filename = 'cpw.out'

       call write_md_cpwout(filename, meta_pwdat, flgcal,                &
     & adot, ntype, natom, nameat, rat, atmass, alatt,                   &
     & emax, nbandin, nx, ny, nz, sx, sy, sz,                            &
     & .TRUE., .FALSE.,                                                  &
     & mxdtyp, mxdatm)

       iotape = 26
       open(unit = iotape, file = 'vesta_pre-relaxed.xsf',               &
     &                  status='UNKNOWN', form='FORMATTED')
       
       call plot_xsf_crys(iotape, .TRUE.,                                &
     &     adot, ntype, natom, nameat, rat,                              &
     &     mxdtyp, mxdatm)

       close(unit = iotape)

       deallocate(natom)
       deallocate(nameat)
       deallocate(rat)
       deallocate(atmass)
       deallocate(force)

       stop

       end program pre_relax_vff


!!>     Dummy subroutine, just calls adot_to_avec 
!!>     to avoid linking to symmetry package
!
!       subroutine adot_to_avec_sym(adot,avec,bvec)
!
!!      Written by Ivo Souza
!!      Modified 15 January 1999. jlm
!!      Converted to f90. JLM
!!      copyright inesc-mn/Jose Luis Martins/Ivo Souza
!
!!      Version 4.94 of cpw
!!      Version 1.5 of md
!
!       implicit none
!       integer, parameter          :: REAL64 = selected_real_kind(12)
!
!!      input
!
!       real(REAL64), intent(in)           ::  adot(3,3)                  !< metric in direct space in atomic units (Bohr radius)
!
!!      output
!
!       real(REAL64), intent(out)          ::  avec(3,3)                  !< primitive lattice vectors that generate adot in canonical orientation
!       real(REAL64), intent(out)          ::  bvec(3,3)                  !< reciprocal lattice vectors
!
!       call adot_to_avec(adot,avec,bvec)
!
!       return
!       end subroutine adot_to_avec_sym


!>     performs one step of molecular dynamics (md version)

       subroutine move_md_atom_cell(flgcal,newcel,iconv,istmd,           &
     & rat,vat,adot,vadot,energy,force,stress,                           &
     & tstep,ekin,ekcell,strext,press,celmas,                            &
     & beta,tempk,iseed,                                                 &
     & mkeep,pgtol,dxmax,                                                &
     & rat1,frc1,adot1,frcel1,                                           &
     & ntype,natom,nameat,atmass,                                        &
     & mxdtyp,mxdatm)

!      modified June 5 2019.  JLM
!      added move_epi_langevin, August 2019. JLM
!      copyright INESC-MN/ Jose Luis Martins

       implicit none

       integer, parameter :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdtyp                     !<  array dimension of types of atoms
       integer, intent(in)                ::  mxdatm                     !<  array dimension of number of atoms of a given type

       character(len=6), intent(in)       ::  flgcal                     !<  type of md calculation

       real(REAL64), intent(in)           ::  energy                     !<  energy in hartree
       real(REAL64), intent(in)           ::  force(3,mxdatm,mxdtyp)     !<  d energy / d rat,  force on the n-th atom of type i (contravarian components, hartree/bohr)
       real(REAL64), intent(in)           ::  stress(3,3)                !<  d energy / d adot,  stress tensor (contravariant)

       integer,intent(in)                 ::  istmd                      !<  md step. Equal to 1 in first step of molecular dynamics
       real(REAL64), intent(in)           ::  tstep                      !<  molecular dynamics time step (in a.u.)

       real(REAL64), intent(in)           ::  tempk                      !<  ionic temperature (in Kelvin)
       real(REAL64), intent(in)           ::  beta                       !<  friction coefficient/mass (in a.u.)

       integer, intent(in)                ::  mkeep                      !<  number of iterations remembered by lbfgs
       real(REAL64), intent(in)           ::  pgtol                      !<  tolerance for minimization convergence
       real(REAL64), intent(in)           ::  dxmax                      !<  maximum step size

       real(REAL64), intent(in)           ::  strext(3,3)                !<  external applied stress
       real(REAL64), intent(in)           ::  press                      !<  external pressure
       real(real64), intent(in)           ::  celmas                     !<  fictitious cell mass       

       integer, intent(in)                ::  ntype                      !<  number of types of atoms
       integer, intent(in)                ::  natom(mxdtyp)              !<  number of atoms of type i
       character(len=2), intent(in)       ::  nameat(mxdtyp)             !<  chemical symbol for the type i
       real(REAL64), intent(in)           ::  atmass(mxdtyp)             !<  atomic mass of atoms of type i

!      input and output
       
       real(REAL64), intent(inout)        ::  rat(3,mxdatm,mxdtyp)       !<  lattice coordinates of atom j of type i
       real(REAL64), intent(inout)        ::  vat(3,mxdatm,mxdtyp)       !<  d rat / d t  velocity in lattice coordinates of atom j of type i

       real(REAL64), intent(inout)        ::  adot(3,3)                  !<  metric in real space
       real(REAL64), intent(inout)        ::  vadot(3,3)                 !<  d adot / d t  rate of change of metric

       integer, intent(inout)             ::  iseed                      !<  seed for rndom number generator

       real(REAL64), intent(inout)        ::  rat1(3,mxdatm,mxdtyp)      !<  previous value of lattice coordinates of atom j of type i
       real(REAL64), intent(inout)        ::  frc1(3,mxdatm,mxdtyp)      !<  previous value of force on atom j of type i
       real(REAL64), intent(inout)        ::  adot1(3,3)                 !<  previous value of adot
       real(REAL64), intent(inout)        ::  frcel1(3,3)                !<  previous cell "force" (covariant components)

!      output

       real(REAL64), intent(out)          ::  ekin                       !<  kinetic energy of the nuclei (atomic units)
       real(REAL64), intent(out)          ::  ekcell                     !<  kinetic energy in Hartree of cell (fictitious)C

       logical, intent(out)               ::  newcel                     !<  indicates that adot has changed on output
       integer, intent(out)               ::  iconv                      !<  iconv = 1, lbfgs converged of flgcal = one

!      local variables

       integer               ::  minrat                     !  if =1 minimize with respect to atomic positions
       integer               ::  minstr                     !  if =1 minimize with respect to all adot variables, if =2 minimize with respect to adot(3,3)
       integer, save         ::  ilbfgs = 0                 !  minimization step in lbfgs           





       if(flgcal == 'ONE   ') then
         iconv = 1
       else

         iconv = 0
         ekcell = 0.0

         if(flgcal == 'MICRO ') then
       
           call move_md_micro(rat,vat,force,istmd,tstep,ekin,            &
     &     rat1,frc1,                                                    &
     &     ntype,natom,atmass,adot,                                      &
     &     mxdatm,mxdtyp)

           newcel = .FALSE.

         elseif(flgcal == 'LANG  ') then

           call move_md_langevin(rat,vat,force,istmd,tstep,ekin,         &
     &     beta,tempk,iseed,                                             &
     &     rat1,frc1,                                                    &
     &     ntype,natom,atmass,adot,                                      &
     &     mxdatm,mxdtyp)

           newcel = .FALSE.

         elseif(flgcal == 'LBFSYM' .or. flgcal == 'VCSLBF' .or.          &
     &          flgcal == 'EPILBF') then

           minrat = 1
           if(flgcal == 'LBFSYM') then
             minstr = 0
             newcel = .FALSE.
           elseif(flgcal == 'VCSLBF') then
             minstr = 1
             newcel = .TRUE.
           elseif(flgcal == 'EPILBF') then
             minstr = 2
             newcel = .TRUE.
           endif
           
!           if(istmd == ist0 + 1) ilbfgs = 0 

           call move_lbfgs_call(energy,rat,force,adot,stress,               &
     &     strext,press,ilbfgs,iconv,minrat,minstr,                         &
     &     ntype,natom,                                                     &
     &     mkeep,pgtol,dxmax,                                               &
     &     mxdtyp,mxdatm)
     
           if(iconv < 0) then
             ilbfgs = 0
           else
             ilbfgs = ilbfgs + 1
           endif

         elseif(flgcal == 'VCSMIC') then

           call move_vcs_micro(rat,vat,adot,vadot,force,stress,          &
     &     istmd,tstep,ekin,ekcell,strext,press,celmas,                  &
     &     rat1,frc1,adot1,frcel1,                                       &
     &     ntype,natom,atmass,                                           &
     &     mxdatm,mxdtyp)

         elseif(flgcal == 'VCSLNG') then

           call move_vcs_langevin(rat,vat,adot,vadot,force,stress,       &
     &     istmd,tstep,ekin,ekcell,strext,press,celmas,                  &
     &     beta,tempk,iseed,                                             &
     &     rat1,frc1,adot1,frcel1,                                       &
     &     ntype,natom,atmass,                                           &
     &     mxdatm,mxdtyp)

         elseif(flgcal == 'EPILNG') then

           call move_epi_langevin(rat,vat,adot,vadot,force,stress,       &
     &     istmd,tstep,ekin,ekcell,strext,press,celmas,                  &
     &     beta,tempk,iseed,                                             &
     &     rat1,frc1,adot1,frcel1,                                       &
     &     ntype,natom,atmass,                                           &
     &     mxdatm,mxdtyp)

         endif

!        prints velocities

         call move_print_velocity(flgcal,                                &
     &   adot,ntype,natom,nameat,vat,vadot,                              &
     &   mxdtyp,mxdatm)

       endif

       return
       end subroutine move_md_atom_cell


!>     prints the information about the parameters used in the calculation

       subroutine print_md_parameters(flgcal,                            &
     & tempk,tempinik,nstep,tstep,beta,iseed,                            &
     & press,strext,celmas)

!      written 20 september 2002. jlm
!      modified 6 february 2008. jlm
!      modified 18 february 2008. jlm
!      Modified f90, tblaha, 24 Outubro 2015. JLM
!      copyright INESC-MN/Jose Luis Martins

!      version 4.68

       implicit none
       
       integer, parameter          :: real64 = selected_real_kind(12)

!      intput

       character(len=6), intent(in)       ::  flgcal                     !<  type of calculation
       
       real(REAL64), intent(in)           ::  tempk                      !<  ionic temperature (in Kelvin)
       real(REAL64), intent(in)           ::  tempinik                   !<  initial ionic temperature (in Kelvin)
       integer, intent(in)                ::  nstep                      !<  number of MD or minimization steps
       real(REAL64), intent(in)           ::  tstep                      !<  time step for MD (in atomic units)
       real(REAL64), intent(in)           ::  beta                       !<  friction coefficient/mass (in a.u.)
       integer, intent(in)                ::  iseed                      !<  initial seed for random number generator        
       
       real(REAL64), intent(in)           ::  press                      !<  external pressure
       real(REAL64), intent(in)           ::  strext(3,3)                !<  external applied stress
       real(real64), intent(in)           ::  celmas                     !<  fictitious cell mass

!      parameters

       real(REAL64), parameter :: AUTOGPA = 29421.58_REAL64

!      counters

       integer   ::  i, j


       write(6,*)
       write(6,*)
       

       if(flgcal == 'MICRO ' .or. flgcal == 'LANG  ' .or.                &
     &    flgcal == 'VCSLNG' .or. flgcal == 'VCSMIC') then
         write(6,'("    The initial ion temperature is ",f13.3,          &
     &     " Kelvin")') tempinik
         write(6,'("    The time step is ",f12.3," au ",                 &
     & "      ( * 2.4 10**-17 s)")') tstep
         write(6,'("    The number of MD steps is : ",i10)') nstep
       endif
       
       if(flgcal == 'LANG  ' .or. flgcal == 'VCSLNG') then
         write(6,'("    The thermal bath temperature is ",f12.3,         &
     &     " Kelvin")') tempk
         write(6,'("    The relaxation time is",f10.2," time steps")')   &
     &        1.0/(beta*tstep)
         write(6,'("    The random number seed is : ",i10)') iseed
       endif

       if(flgcal == 'VCSLNG' .or. flgcal == 'VCSLBF'                     &
     &    .or. flgcal == 'EPILBF' .or. flgcal == 'VCSMIC') then
         write(6,'("    The applied pressure is ",f12.3," GPa.  ",       &
     &    "The stress tensor is",3f10.2,/,65x,3f10.2,/,65x,3f10.2)')     &
     &        press*AUTOGPA, ((strext(i,j)*AUTOGPA,i=1,3),j=1,3)
       endif
       if(flgcal == 'VCSLNG' .or. flgcal == 'VCSMIC') then
         write(6,'("    The fictitious cell mass is",f12.3)') celmas 
       endif
       
       if(flgcal == 'VCSLBF' .or. flgcal == 'LBFSYM'                     &
     &    .or. flgcal == 'EPILBF') then
         write(6,'("    Maximum number of steps is : ",i10)') nstep
       endif
       write(6,*)
       
       return
       end subroutine print_md_parameters


!>     reads the parameters for the md calculation

       subroutine read_md_esdf(ipr,fname,flgcal,flgmod,                  &
     & adot,ntype,natom,nameat,rat,atmass,alatt,lgeom,                   &
     & emax,nbandin,nx,ny,nz,sx,sy,sz,                                   &
     & meta_cpw2000,                                                     &
     & tempk,tempinik,nstep,tstep,beta,iseed,pgtol,dxmax,                &
     & press,strext,celmas,                                              &
     & mxdtyp,mxdatm)

!      Adapted for md June 3 2019, JLM
!
!      copyright inesc-mn/Jose Luis Martins/Carlos Loia Reis


!      version 1.5

       use esdf

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdtyp                     !<  array dimension of types of atoms
       integer, intent(in)                ::  mxdatm                     !<  array dimension of types of atoms


       integer, intent(in)                ::  ipr                        !<  controls print verbosity

       character(len=*), intent(in)       ::  fname                      !<  filename with input data

!      output

       character(len=6), intent(out)      ::  flgcal                     !<  type of calculation
       character(len=6), intent(out)      ::  flgmod                     !<  classical potential model

       real(REAL64), intent(out)          ::  adot(3,3)                  !<  metric in direct space
       integer, intent(out)               ::  ntype                      !<  number of types of atoms
       integer, intent(out)               ::  natom(mxdtyp)              !<  number of atoms of type i
       character(len=2), intent(out)      ::  nameat(mxdtyp)             !<  chemical symbol for the type i
       real(REAL64), intent(out)          ::  rat(3,mxdatm,mxdtyp)       !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i
       real(REAL64), intent(out)          ::  atmass(mxdtyp)             !<  atomic mass (in a.u.) of atom of type i 

       real(REAL64), intent(out)          ::  alatt                      !<  lattice constant
       logical, intent(out)               ::  lgeom                      !<  indicates if geometry was successfully read.

!      The next group of variables are not used im md, just passed to the next cpw.in

       real(REAL64), intent(out)          ::  emax                       !<  kinetic energy cutoff of plane wave expansion (Hartree).
       integer, intent(out)               ::  nx, ny, nz                 !<  size of the integration mesh in k-space (nx*ny*nz)
       real(REAL64), intent(out)          ::  sx, sy, sz                 !<  offset of the integration mesh (usually 0.5)
       integer, intent(out)               ::  nbandin                    !<  target for number of bands      
       character(len=250), intent(out)    ::  meta_cpw2000               !<  metadata coded in cpw2000
       
       real(REAL64), intent(out)          ::  tempk                      !<  ionic temperature (in Kelvin)
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

!      local variables

       logical                    ::  lrede
       logical                    ::  lbz
       character(len=6)           ::  flgkeat
       logical                    ::  flgkplusg                  !  unused
       real(REAL64)               ::  epskplusg                  !  unused

!      parameters

       real(REAL64), parameter ::  UM = 1.0_REAL64


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
!      RSTRT:  restarts an old calculation from the point where it was interrupted

       flgcal = 'ONE   '
!      flgcal = 'MICRO '
!      flgcal = 'LANG  '
!      flgcal = 'LBFSYM'
!      flgcal = 'VCSLNG'
!      flgcal = 'VCSMIC'
!      flgcal = 'VCSLBF'
!      flgcal = 'EPILBF'
!      flgcal = 'RSTRT '

       flgcal = esdf_string('MD.TypeOfRun',flgcal)
       call chrcap(flgcal,6)

!      may be overridden by command line

       if (command_argument_count() >=1) then
         call get_command_argument(1,flgcal)
       endif

       flgmod = 'LENJON'
       if(esdf_defined('MD.PotentialModel')) then
         flgmod = esdf_string('MD.PotentialModel',flgmod)
         call chrcap(flgmod,6)
       else
         write(6,*)
         write(6,*) '   WARNING'
         write(6,*)
         write(6,*) '   Assuming Lennard Jones Potential'    
       endif

       call read_esdf_crystal(ipr,                                       &
     & adot,ntype,natom,nameat,rat,                                      &
     & atmass,alatt,lgeom,                                               &
     & mxdtyp,mxdatm)

       call read_esdf_bz(ipr,                                            &
     & emax,nbandin,nx,ny,nz,sx,sy,sz,lbz)

       call read_esdf_rede(ipr,meta_cpw2000,lbz)

       call read_esdf_md(ipr,                                            &
     & tempk,tempinik,nstep,tstep,flgkeat,                               &
     & beta,iseed,pgtol,dxmax,press,strext,celmas,flgkplusg,epskplusg)

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
         write(6,*)
         write(6,'("   The value of flgmod is: ",a6)') flgmod
         write(6,*)
         write(6,*)
       endif

!      Overrides number of steps

       if(flgcal == 'ONE   ') nstep = 1

       return
       end subroutine read_md_esdf


!>     table of keating vff constants as tabulated by
!>     J.L.Martins and A.Zunger PhysRevB 30, 6217,1984.
!>     Distances in angstroms and Force constants in N/m from table.
!>     On output they are converted to atomic (Hartree/bohr) units.

!>     If the force constants are unavailable they are given a HUGE value
!>     and the interatomic distance is negative so it is obvious in the output
!>     in the case vff_check is not used.

       subroutine vff_constants(iprint, iowrite,                         &
     &                     ntype, nameat, alfa, dist, beta ,mxdtyp)

!      written January 2013, J.L.Martins
!      Modified, documentation, details, printing, December 2019. JLM
!      Added guessed constants for Pb, May 22, 2020.  JLM
!      and Arithmetic mean instead of geometric for carbides. JLM
!      Added constants from doi:10.1016/j.physe.2009.11.035
!      copyright INESC-MN/Jose Luis Martins/C.L. Reis

!      version 4.99 of pw
!      version 1.6.0 of md

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input variables

       integer, intent(in)         :: mxdtyp                             !<  maximum number of type of atoms

       integer, intent(in)         :: iprint                             !<  if iprint == 0 does not print details
       integer, intent(in)         :: iowrite                            !<  tape number

       integer, intent(in)         :: ntype                              !<  number of types of atoms
       character(len=2),intent(in) :: nameat(mxdtyp)                     !<  chemical symbol of type of atom

!      output variables

       real(REAL64), intent(out)   :: alfa((mxdtyp*(mxdtyp+1))/2)        !<  keating coefficient for bond
       real(REAL64), intent(out)   :: dist((mxdtyp*(mxdtyp+1))/2)        !<  equilibrium distance for bond
       real(REAL64), intent(out)   :: beta(mxdtyp,(mxdtyp*(mxdtyp+1))/2) !<  keating coefficient for angle

!      counters

       integer :: n,m,k,nn2,mm2,nm2,mk2,kk2

!      physical constants

       real(REAL64), parameter  :: BOHR = 0.52917721E-10
       real(REAL64), parameter  :: EV = 27.211385
       real(REAL64), parameter  :: HARTREE = EV * 1.6021765E-19
       real(REAL64), parameter  :: ANG = 1E-10

       real(REAL64), parameter  :: UM = 1.0_REAL64

       
!      same species alpha,dist,beta (group IV)

       do n=1,ntype
         nn2 = (n*(n+1))/2
         if(nameat(n) == 'C ' .or. nameat(n) == ' C') then
           alfa(nn2) = 129.33_REAL64
           dist(nn2) = 1.545_REAL64
           beta(n,nn2) = 0.655_REAL64*alfa(nn2)
         elseif(nameat(n) == 'Si') then
           alfa(nn2) = 48.50_REAL64
           dist(nn2) = 2.352_REAL64
           beta(n,nn2) = 0.285_REAL64*alfa(nn2)
         elseif(nameat(n) == 'Ge') then
           alfa(nn2) = 38.67_REAL64
           dist(nn2) = 2.450_REAL64
           beta(n,nn2) = 0.294_REAL64*alfa(nn2)
         elseif(nameat(n) == 'Sn') then
           alfa(nn2) = 25.45_REAL64
           dist(nn2) = 2.810_REAL64
           beta(n,nn2) = 0.253_REAL64*alfa(nn2)
         elseif(nameat(n) == 'Pb') then
           alfa(nn2) = 22_REAL64
           dist(nn2) = 2.99_REAL64
           beta(n,nn2) = 0.25_REAL64*alfa(nn2)
         else
           alfa(nn2) = 100000.0
           dist(nn2) = -1.0
           beta(n,nn2) = alfa(nn2)
         endif
       enddo
       
       if(ntype > 1) then
       
!        pairs of species alpha,dist,beta
         
         do n=1,ntype-1
         do m=n+1,ntype
           nm2 = (m*(m-1))/2+n
           nn2 = (n*(n+1))/2
           mm2 = (m*(m+1))/2

!          group IV

           if(((nameat(n) == 'C ' .or. nameat(n) == ' C') .and.           &
     &          nameat(m) == 'Si') .or. (nameat(n) == 'Si' .and.          &
     &          (nameat(m) == 'C ' .or. nameat(m) == ' C'))) then
             alfa(nm2) = 88.0_REAL64
             dist(nm2) = 1.888_REAL64
             beta(n,mm2) = 0.54_REAL64*alfa(nm2)
           elseif(((nameat(n) == 'C ' .or. nameat(n) == ' C') .and.       &
     &          nameat(m) == 'Ge') .or. (nameat(n) == 'Ge' .and.          &
     &          (nameat(m) == 'C ' .or. nameat(m) == ' C'))) then
             alfa(nm2) = (alfa(nn2)+alfa(mm2))/2
             dist(nm2) = (dist(nn2)*dist(mm2))/2
             beta(n,mm2) = (beta(n,nn2)+beta(m,mm2))/2
           elseif(((nameat(n) == 'C ' .or. nameat(n) == ' C') .and.       &
     &          nameat(m) == 'Sn') .or. (nameat(n) == 'Sn' .and.          &
     &          (nameat(m) == 'C ' .or. nameat(m) == ' C'))) then
             alfa(nm2) = (alfa(nn2)+alfa(mm2))/2
             dist(nm2) = (dist(nn2)*dist(mm2))/2
             beta(n,mm2) = (beta(n,nn2)+beta(m,mm2))/2
           elseif(((nameat(n) == 'C ' .or. nameat(n) == ' C') .and.       &
     &          nameat(m) == 'Pb') .or. (nameat(n) == 'Pb' .and.          &
     &          (nameat(m) == 'C ' .or. nameat(m) == ' C'))) then
             alfa(nm2) = (alfa(nn2)+alfa(mm2))/2
             dist(nm2) = (dist(nn2)*dist(mm2))/2
             beta(n,mm2) = (beta(n,nn2)+beta(m,mm2))/2
           elseif((nameat(n) == 'Si' .and. nameat(m) == 'Ge') .or.        &
     &          (nameat(n) == 'Ge' .and. nameat(m) == 'Si')) then
             alfa(nm2) = (alfa(nn2)+alfa(mm2))/2
             dist(nm2) = sqrt(dist(nn2)*dist(mm2))
             beta(n,mm2) = (beta(n,nn2)+beta(m,mm2))/2
           elseif((nameat(n) == 'Si' .and. nameat(m) == 'Sn') .or.        &
     &          (nameat(n) == 'Sn' .and. nameat(m) == 'Si')) then
             alfa(nm2) = (alfa(nn2)+alfa(mm2))/2
             dist(nm2) = sqrt(dist(nn2)*dist(mm2))
             beta(n,mm2) = (beta(n,nn2)+beta(m,mm2))/2
           elseif((nameat(n) == 'Si' .and. nameat(m) == 'Pb') .or.        &
     &          (nameat(n) == 'Pb' .and. nameat(m) == 'Si')) then
             alfa(nm2) = (alfa(nn2)+alfa(mm2))/2
             dist(nm2) = sqrt(dist(nn2)*dist(mm2))
             beta(n,mm2) = (beta(n,nn2)+beta(m,mm2))/2
           elseif((nameat(n) == 'Ge' .and. nameat(m) == 'Sn') .or.        &
     &          (nameat(n) == 'Sn' .and. nameat(m) == 'Ge')) then
             alfa(nm2) = (alfa(nn2)+alfa(mm2))/2
             dist(nm2) = sqrt(dist(nn2)*dist(mm2))
             beta(n,mm2) = (beta(n,nn2)+beta(m,mm2))/2
           elseif((nameat(n) == 'Ge' .and. nameat(m) == 'Pb') .or.        &
     &          (nameat(n) == 'Pb' .and. nameat(m) == 'Ge')) then
             alfa(nm2) = (alfa(nn2)+alfa(mm2))/2
             dist(nm2) = sqrt(dist(nn2)*dist(mm2))
             beta(n,mm2) = (beta(n,nn2)+beta(m,mm2))/2
           elseif((nameat(n) == 'Sn' .and. nameat(m) == 'Pb') .or.        &
     &          (nameat(n) == 'Pb' .and. nameat(m) == 'Sn')) then
             alfa(nm2) = (alfa(nn2)+alfa(mm2))/2
             dist(nm2) = sqrt(dist(nn2)*dist(mm2))
             beta(n,mm2) = (beta(n,nn2)+beta(m,mm2))/2

!          Group III-V

           elseif((nameat(n) == 'Al' .and. (nameat(m) == 'N ' .or.        &
     &           nameat(m) == ' N')) .or. ((nameat(n) == 'N ' .or.        &
     &           nameat(n) == ' N') .and. nameat(m) == 'Al')) then
             alfa(nm2) = 83.8_REAL64
             dist(nm2) = 1.89_REAL64
             beta(n,mm2) = 0.236_REAL64*alfa(nm2)
           elseif((nameat(n) == 'Al' .and. (nameat(m) == 'P ' .or.        &
     &           nameat(m) == ' P')) .or. ((nameat(n) == 'P ' .or.        &
     &           nameat(n) == ' P') .and. nameat(m) == 'Al')) then
             alfa(nm2) = 47.29_REAL64
             dist(nm2) = 2.367_REAL64
             beta(n,mm2) = 0.192_REAL64*alfa(nm2)
           elseif((nameat(n) == 'Al' .and. nameat(m) == 'As') .or.        &
     &          (nameat(n) == 'As' .and. nameat(m) == 'Al')) then
             alfa(nm2) = 43.05_REAL64
             dist(nm2) = 2.451_REAL64
             beta(n,mm2) = 0.229_REAL64*alfa(nm2)
           elseif((nameat(n) == 'Al' .and. nameat(m) == 'Sb') .or.        &
     &          (nameat(n) == 'Sb' .and. nameat(m) == 'Al')) then
             alfa(nm2) = 35.35_REAL64
             dist(nm2) = 2.656_REAL64
             beta(n,mm2) = 0.192_REAL64*alfa(nm2)

           elseif((nameat(n) == 'Ga' .and. (nameat(m) == 'N ' .or.        &
     &           nameat(m) == ' N')) .or. ((nameat(n) == 'N ' .or.        &
     &           nameat(n) == ' N') .and. nameat(m) == 'Ga')) then
             alfa(nm2) = 88.4_REAL64
             dist(nm2) = 1.95_REAL64
             beta(n,mm2) = 0.237_REAL64*alfa(nm2)
           elseif((nameat(n) == 'Ga' .and. (nameat(m) == 'P ' .or.        &
     &           nameat(m) == ' P')) .or. ((nameat(n) == 'P ' .or.        &
     &           nameat(n) == ' P') .and. nameat(m) == 'Ga')) then
             alfa(nm2) = 47.32_REAL64
             dist(nm2) = 2.360_REAL64
             beta(n,mm2) = 0.221_REAL64*alfa(nm2)
           elseif((nameat(n) == 'Ga' .and. nameat(m) == 'As') .or.        &
     &          (nameat(n) == 'As' .and. nameat(m) == 'Ga')) then
             alfa(nm2) = 41.19_REAL64
             dist(nm2) = 2.448_REAL64
             beta(n,mm2) = 0.217_REAL64*alfa(nm2)
           elseif((nameat(n) == 'Ga' .and. nameat(m) == 'Sb') .or.        &
     &          (nameat(n) == 'Sb' .and. nameat(m) == 'Ga')) then
             alfa(nm2) = 33.16_REAL64
             dist(nm2) = 2.640_REAL64
             beta(n,mm2) = 0.218_REAL64*alfa(nm2)

           elseif((nameat(n) == 'In' .and. (nameat(m) == 'N ' .or.        &
     &           nameat(m) == ' N')) .or. ((nameat(n) == 'N ' .or.        &
     &           nameat(n) == ' N') .and. nameat(m) == 'In')) then
             alfa(nm2) = 67.4_REAL64
             dist(nm2) = 2.16_REAL64
             beta(n,mm2) = 0.149_REAL64*alfa(nm2)
           elseif((nameat(n) == 'In' .and. (nameat(m) == 'P ' .or.        &
     &           nameat(m) == ' P')) .or. ((nameat(n) == 'P ' .or.        &
     &           nameat(n) == ' P') .and. nameat(m) == 'In')) then
             alfa(nm2) = 43.04_REAL64
             dist(nm2) = 2.541_REAL64
             beta(n,mm2) = 0.145_REAL64*alfa(nm2)
           elseif((nameat(n) == 'In' .and. nameat(m) == 'As') .or.        &
     &          (nameat(n) == 'As' .and. nameat(m) == 'In')) then
             alfa(nm2) = 35.18_REAL64
             dist(nm2) = 2.622_REAL64
             beta(n,mm2) = 0.156_REAL64*alfa(nm2)
           elseif((nameat(n) == 'In' .and. nameat(m) == 'Sb') .or.        &
     &          (nameat(n) == 'Sb' .and. nameat(m) == 'In')) then
             alfa(nm2) = 26.61_REAL64
             dist(nm2) = 2.805_REAL64
             beta(n,mm2) = 0.161_REAL64*alfa(nm2)

!          Group II-VI

           elseif((nameat(n) == 'Zn' .and. (nameat(m) == 'S ' .or.        &
     &           nameat(m) == ' S')) .or. ((nameat(n) == 'S ' .or.        &
     &           nameat(n) == ' S') .and. nameat(m) == 'Zn')) then
             alfa(nm2) = 44.92_REAL64
             dist(nm2) = 2.342_REAL64
             beta(n,mm2) = 0.107_REAL64*alfa(nm2)
           elseif((nameat(n) == 'Zn' .and. nameat(m) == 'Se') .or.        &
     &          (nameat(n) == 'Se' .and. nameat(m) == 'Zn')) then
             alfa(nm2) = 35.24_REAL64
             dist(nm2) = 2.454_REAL64
             beta(n,mm2) = 0.120_REAL64*alfa(nm2)
           elseif((nameat(n) == 'Zn' .and. nameat(m) == 'Te') .or.        &
     &          (nameat(n) == 'Te' .and. nameat(m) == 'Zn')) then
             alfa(nm2) = 31.35_REAL64
             dist(nm2) = 2.637_REAL64
             beta(n,mm2) = 0.142_REAL64*alfa(nm2)
             
           elseif((nameat(n) == 'Cd' .and. nameat(m) == 'Te') .or.        &
     &          (nameat(n) == 'Te' .and. nameat(m) == 'Cd')) then
             alfa(nm2) = 29.02_REAL64
             dist(nm2) = 2.806_REAL64
             beta(n,mm2) = 0.084_REAL64*alfa(nm2)

           elseif((nameat(n) == 'Hg' .and. (nameat(m) == 'S ' .or.        &
     &           nameat(m) == ' S')) .or. ((nameat(n) == 'S ' .or.        &
     &           nameat(n) == ' S') .and. nameat(m) == 'Hg')) then
             alfa(nm2) = 41.33_REAL64
             dist(nm2) = 2.534_REAL64
             beta(n,mm2) = 0.062_REAL64*alfa(nm2)
           elseif((nameat(n) == 'Hg' .and. nameat(m) == 'Se') .or.        &
     &          (nameat(n) == 'Se' .and. nameat(m) == 'Hg')) then
             alfa(nm2) = 36.35_REAL64
             dist(nm2) = 2.634_REAL64
             beta(n,mm2) = 0.065_REAL64*alfa(nm2)
           elseif((nameat(n) == 'Hg' .and. nameat(m) == 'Te') .or.        &
     &          (nameat(n) == 'Te' .and. nameat(m) == 'Hg')) then
             alfa(nm2) = 27.97_REAL64
             dist(nm2) = 2.798_REAL64
             beta(n,mm2) = 0.092_REAL64*alfa(nm2)
             
!          Group I-VII

           elseif((nameat(n) == 'Cu' .and. nameat(m) == 'Cl') .or.        &
     &          (nameat(n) == 'Cl' .and. nameat(m) == 'Cu')) then
             alfa(nm2) = 22.9_REAL64
             dist(nm2) = 2.341_REAL64
             beta(n,mm2) = 0.04_REAL64*alfa(nm2)
           elseif((nameat(n) == 'Cu' .and. nameat(m) == 'Br') .or.        &
     &          (nameat(n) == 'Br' .and. nameat(m) == 'Cu')) then
             alfa(nm2) = 23.1_REAL64
             dist(nm2) = 2.464_REAL64
             beta(n,mm2) = 0.057_REAL64*alfa(nm2)
           elseif((nameat(n) == 'Cu' .and. (nameat(m) == 'I ' .or.        &
     &           nameat(m) == ' I')) .or. ((nameat(n) == 'I ' .or.        &
     &           nameat(n) == ' I') .and. nameat(m) == 'Cu')) then
             alfa(nm2) = 22.5_REAL64
             dist(nm2) = 2.617_REAL64
             beta(n,mm2) = 0.091_REAL64*alfa(nm2)

           else
             alfa(nm2) = 100000.0
             dist(nm2) = -1.0
             beta(n,mm2) = UM*alfa(nm2)
           endif

           beta(m,nn2) = beta(n,mm2)
           beta(n,nm2) = sqrt(beta(n,mm2)*beta(n,nn2))
           beta(m,nm2) = sqrt(beta(m,nn2)*beta(m,mm2))

         enddo
         enddo

         if(ntype > 2) then

!          triplets of species beta

           do n=1,ntype
             do m=1,ntype-1
               do k=m+1,ntype
                 if(m /= n .and. k /= n) then
                   mk2 = (k*(k-1))/2 + m
                   mm2 = (m*(m+1))/2
                   kk2 = (k*(k+1))/2
                   beta(n,mk2) = sqrt(beta(n,mm2)*beta(n,kk2))
                 endif
               enddo
             enddo
           enddo

         endif

       endif

       if(iprint /= 0) then

         write(iowrite,*)
         write(iowrite,*)
         write(iowrite,*) '    alpha  constants (N/m) '
         write(iowrite,*)
         write(iowrite,'(6x,10(10x,a2))') (nameat(m),m=1,ntype)
         do n=1,ntype
           write(iowrite,'(2x,a2,4x,10f12.2)') nameat(n),                &
     &                     (alfa((n*(n-1))/2+k),k=1,n)
         enddo

         write(iowrite,*)
         write(iowrite,*) '    distances (Angstrom)'
         write(iowrite,*)
         write(iowrite,'(6x,10(10x,a2))') (nameat(m),m=1,ntype)
         do n=1,ntype
           write(iowrite,'(2x,a2,4x,10f12.3)') nameat(n),                &
     &                     (dist((n*(n-1))/2+k),k=1,n)
         enddo     

         write(iowrite,*)
         write(iowrite,*) '    beta constants (N/m)'
         write(iowrite,*)
         do n=1,ntype
           write(iowrite,*)
           write(iowrite,'(" corner atom:  ",a2)') nameat(n)
           write(iowrite,*)
           write(iowrite,'(6x,10(10x,a2))') (nameat(m),m=1,ntype)
           do m=1,ntype
             write(iowrite,'(2x,a2,4x,10f12.2)') nameat(m),              &
     &                     (beta(n,(m*(m-1))/2+k),k=1,m)
           enddo
         enddo

       endif

!      convert to atomic units

       do n=1,(ntype*(ntype+1))/2
         alfa(n) = alfa(n) * BOHR*BOHR / HARTREE
         dist(n) = dist(n) * ANG / BOHR
         do m=1,ntype
           beta(m,n) = beta(m,n) * BOHR*BOHR / HARTREE
         enddo
       enddo

       return
       end subroutine vff_constants




!>     writes the next cpw.in file

       subroutine write_md_cpwout(filename, meta_pwdat, flgcal,          &
     & adot, ntype, natom, nameat, rat, atmass, alatt,                   &
     & emax, nbandin, nx, ny, nz, sx, sy, sz,                            &
     & lkeat, ltbl,                                                      &
     & mxdtyp, mxdatm)

!      Adapted June 2017. JLM
!      Bug squashed (metadata not from rede) September 2017.
!      Adapted for md
!      copyright  J.L.Martins, INESC-MN.

!      version 1.5 of md

       implicit none

       integer, parameter  :: REAL64 = selected_real_kind(12)

!      input:

       integer, intent(in)                ::  mxdtyp                     !<  array dimension of types of atoms
       integer, intent(in)                ::  mxdatm                     !<  array dimension of types of atoms

       character(len=20), intent(in)      ::  filename                   !<  name of output file
       character(len=250), intent(in)     ::  meta_pwdat                 !<  metadata from cpw_in or PW.DAT

       character(len=6), intent(in)       ::  flgcal                     !<  type of md calculation

       real(REAL64), intent(in)           ::  adot(3,3)                  !<  metric in direct space
       integer, intent(in)                ::  ntype                      !<  number of types of atoms
       integer, intent(in)                ::  natom(mxdtyp)              !<  number of atoms of type i
       character(len=2), intent(in)       ::  nameat(mxdtyp)             !<  chemical symbol for the type i
       real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)       !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i
       real(REAL64), intent(in)           ::  atmass(mxdtyp)             !<  atomic mass (in a.u.) of atom of type i 

       real(REAL64), intent(in)           ::  alatt                      !<  lattice constant

       real(REAL64), intent(in)           ::  emax                       !<  kinetic energy cutoff of plane wave expansion (Hartree).
       integer, intent(in)                ::  nbandin                    !<  target for number of bands      

       integer, intent(in)                ::  nx, ny, nz                 !<  size of the integration mesh in k-space (nx*ny*nz)
       real(REAL64), intent(in)           ::  sx, sy, sz                 !<  offset of the integration mesh (usually 0.5)

       logical, intent(in)                ::  lkeat                      !<  sets the keating option
       logical, intent(in)                ::  ltbl                       !<  toggles XC between TBL and CA

!      local:
       
       real(REAL64)      ::  avec(3,3),bvec(3,3)

       integer            ::  nlatpl                                     ! number of lattice planes
       character(len=3)   ::  vers                                       ! program version
       character(len=50)  ::  title                                      ! title of the calculation
       character(len=140) ::  metadata                                   ! metadata of the calculation

       integer  ::  isl1(3),isl2(3),isl3(3)
       integer  ::  nat
       integer  ::  istart, iend, ioerr

       integer  ::   io                         !  tape number

!      counters

       integer       ::  nt, n, i, j, nc


!      open file

       io = 10
       open(unit = io, file = filename, status='UNKNOWN',                &
     &                  form='FORMATTED')
     
       write(io,'(72a1)') ("#",j=1,72)
       write(io,'("#",70x,"#")')
       write(io,'("#",6x,"cpw.in input file generated by a previous",    &
     &       " mdtest run ",10x,"#")')
       write(io,'("#",70x,"#")')
       write(io,'(72a1)') ("#",j=1,72)

       read(meta_pwdat,'(3(3i4,2x),8x,i3,3x,a3,1x,a50,a140)',            &
     &      iostat=ioerr)                                                &
     &            (isl1(i),i=1,3),(isl2(i),i=1,3),(isl3(i),i=1,3),       &
     &            nlatpl,vers,title,metadata

       if(ioerr == 0) then

         write(io,*)
         write(io,'("#------------------------------------------------")')
         write(io,'("# Rede metadata")')
         write(io,'("#------------------------------------------------")')
         write(io,*)

         write(io,'("%block Rede.Superlattice")')
         write(io,'(12x,3i6)') (isl1(i),i=1,3)
         write(io,'(12x,3i6)') (isl2(i),i=1,3)
         write(io,'(12x,3i6)') (isl3(i),i=1,3)
         write(io,'("%endblock Rede.Superlattice")')
         write(io,*)

         write(io,'("Rede.NumberOfLatticePlanes",3x,i6)') nlatpl
         write(io,*)

         write(io,'("Rede.Version",18x,a3)') vers
         write(io,*)

         write(io,'("Rede.Title",10x,a50)') title
         write(io,*)

         istart = 1
         do i=1,135
           if(metadata(i:i+4) == '#NAME') then
             istart = i+6

             exit

           endif
         enddo

         iend = 140
         do i=1,125
           if(metadata(i:i+4) == '#DATE') then
             iend = i-1
             write(io,'("Rede.Date",21x,a9)') metadata(i+6:i+15)
             write(io,*)

             exit

           endif
         enddo
    
         do i=1,126
           if(metadata(i:i+4) == '#TIME') then
             write(io,'("Rede.Time",21x,a8)') metadata(i+6:i+14)
             write(io,*)

             exit

           endif
         enddo

         write(io,'("Rede.Name",21x,140a1)')                             &
     &                   (metadata(i:i),i=istart,iend)
         write(io,*)

         write(io,'("SystemLabel",19x,140a1)')                           &
     &                   (metadata(i:i),i=istart,iend)
         write(io,*)

       else

         write(io,*)
         write(io,'("SystemLabel",19x,20a1)')                            &
     &                   (meta_pwdat(i:i),i=1,20)
         write(io,*)
       endif

       write(io,*)
       write(io,'("#------------------------------------------------")')
       write(io,'("# Crystal structure")')
       write(io,'("#------------------------------------------------")')
       write(io,*)
       write(io,'("LatticeConstant",10x,f16.8,5x,"bohr")') alatt
       write(io,*)
       
       call adot_to_avec_sym(adot,avec,bvec)

       write(io,'("%block LatticeVectors")')
       write(io,'(3(3x,f16.8))') avec(1,1)/alatt,avec(2,1)/alatt,        &
     &                           avec(3,1)/alatt
       write(io,'(3(3x,f16.8))') avec(1,2)/alatt,avec(2,2)/alatt,        &
     &                           avec(3,2)/alatt
       write(io,'(3(3x,f16.8))') avec(1,3)/alatt,avec(2,3)/alatt,         &
     &                           avec(3,3)/alatt
       write(io,'("%endblock LatticeVectors")')
       write(io,*)

!      Finds total number of atoms

       nat = 0
       do i=1,ntype
         nat = nat + natom(i)
       enddo

       write(io,'("NumberOfSpecies",9x,i8)') ntype
       write(io,*)
       
       write(io,'("NumberOfAtoms",9x,i8)') nat
       write(io,*)
       
       write(io,'("%block Chemical_Species_Label")')
       do i = 1,ntype
         call p_tbl_charge(nameat(i),n)
         write(io,'(i6,3x,i5,3x,a2)') i,n,nameat(i)
       enddo
       write(io,'("%endblock Chemical_Species_Label")')
       write(io,*)
       
       write(io,'("AtomicCoordinatesFormat",5x,"Fractional")')
       write(io,*)
       
       
       write(io,'("%block AtomicCoordinatesAndAtomicSpecies")')
       
       do nt = 1,ntype
       do i = 1,natom(nt)
          write(io,'(3x,3f16.8,4x,i5,5x,"#  ",a2,3x,i5)')                 &
     &          (rat(j,i,nt),j=1,3),nt,nameat(nt),i
       enddo
       enddo
       write(io,'("%endblock AtomicCoordinatesAndAtomicSpecies")')
       write(io,*)

       write(io,'("StructureSource",15x,"mdtest")')
       write(io,*)

       write(io,'("#------------------------------------------------")')
       write(io,'("# Energy cutoff, bands,  and Brillouin mesh")')
       write(io,'("#------------------------------------------------")')
       write(io,*)

       write(io,'("PWEnergyCutoff",13x,f12.4,6x,"hartree")')  emax
       write(io,*)

       write(io,'("NumberOfEigenStates",12x,i8)') nbandin
       write(io,*)

       write(io,'("%block kgrid_Monkhorst_Pack")')
       write(io,'(8x,3(2x,i6),3x,f12.6)')  nx,0,0,sx
       write(io,'(8x,3(2x,i6),3x,f12.6)')  0,ny,0,sy
       write(io,'(8x,3(2x,i6),3x,f12.6)')  0,0,nz,sz
       write(io,'("%endblock kgrid_Monkhorst_Pack")')
       write(io,*)
      

       write(io,*)
       write(io,'("#------------------------------------------------")')
       write(io,'("# Active options")')
       write(io,'("#------------------------------------------------")')
       write(io,*)

       write(io,'("MD.TypeOfRun",18x,a6,8x,"# ONE,EPILBF,MICRO,",        &
     &      "LANG,LBFSYM,VCSLNG,VCSLBF,RSTRT,EPILNG")') flgcal
       write(io,*)

       write(io,'("MD.PotentialModel             KEATNG        ",        &
     &      "# KEATNG (do not use LENJON,LJCLST) ")')
       write(io,*)

       write(io,'("UseSymmetry                   .true.        ",        &
     &      "# .true. , .false. ")')
       write(io,*)

       if(lkeat) then
         write(io,'("MD.UseKeatingCorrections      .true.        ",      &
     &      "# .true. , .false. ")')
       else
         write(io,'("MD.UseKeatingCorrections      .false.       ",      &
     &      "# .true. , .false. ")')
       endif

       write(io,*)

       write(io,'("MD.CG.UseFixedkplusG          .true.        ",        &
     &      "# .true. , .false. ")')
       write(io,*)

       write(io,'("TypeOfScfDiag                 PW            ",        &
     &      "# PW,AO,AOJC,AOJCPW")')
       write(io,*)

       write(io,'("DualApproximation             .true.        ",        &
     &      "#  .true. , .false.")')
       write(io,*)

       if(ltbl) then
         write(io,'("XC.Authors                    TBL           ",      &
     &      "# CA, PBE, TBL")')
       else
         write(io,'("XC.Authors                    CA            ",      &
     &      "# CA, PBE, TBL")')
       endif
       write(io,*)

       write(io,'("Xc.TBL.C                      1.04          ",        &
     &      "# sets Tran-Blaha constant (if negative use calculated)")')
       write(io,*)

       write(io,'("PrintingLevel                 1             ",        &
     &      "# 1, 2, 3")')
       write(io,*)

       write(io,*)
       write(io,'("#------------------------------------------------")')
       write(io,'("# MD Inactive options")')
       write(io,'("#------------------------------------------------")')
       write(io,*)

       write(io,'("#MD.InitialTemperature        300 K           #")')
       write(io,'("#MD.TargetTemperature         300 K           #")')
       write(io,'("#MD.TargetPressure            0 GPa           #")')
       write(io,*)
       write(io,'("#MD.NumberOfSteps             10              #")')
       write(io,'("#MD.LengthTimeStep            2.4 fs          #")')
       write(io,'("#MD.FrictionFracInvTimeStep   20.0            #")')
       write(io,*)
       write(io,'("#MD.CG.Tolerance         0.0001 ''har/bohr''    #")')
       write(io,'("#MD.CG.StepMax                0.01 bohr       #")')
       write(io,'("#MD.CG.FixedkplusGTol    0.01 ''har/bohr''      #")')
       write(io,*)
       write(io,'("#%block MD.TargetStress                       #")')
       write(io,'("   0.0 0.0 0.0                                #")')
       write(io,'("   0.0 0.0 0.0                                #")')
       write(io,'("   0.0 0.0 0.0                                #")')
       write(io,'("#%endblock   MD.TargetStress                  #")')
       write(io,*)
       write(io,'("#MD.CellMass                  10.0            #")')
       write(io,'("#MD.Seed                      76978           #")')

       write(io,*)
       write(io,'("#------------------------------------------------")')
       write(io,'("# Electronic Structure Inactive options")')
       write(io,'("#------------------------------------------------")')
       write(io,*)

       write(io,'("#MaxSCFIterations             20              #")')
       write(io,*)
       write(io,'("#MaxSCFIterations             20              #")')
       write(io,'("#TypeOfPseudoMixing           BROYD1          #",     &
     &      " BROYD1, BFGS#")')
       write(io,*)
       write(io,'("#ElectronicTemperature        1000 K          #")')
       write(io,'("#TypeOfPseudopotential        PSEUKB          #",     &
     &      " PSEUKB")')
       write(io,*)
       write(io,'("#ScfTolerance                 0.00005         #")')
       write(io,'("#DiagTolerance                0.0001          #")')
       write(io,'("#SymmTolerance                1.0E-5          #")')

       close(unit = io)

       return

       end subroutine write_md_cpwout

