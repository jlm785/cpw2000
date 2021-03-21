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

!>     reads the parameters for molecular dynamics calculations
!>     from a file opened by a previous edsf_init

       subroutine read_esdf_md(ipr,                                      &
     & tempk,tempinik,nstep,tstep,flgkeat,                               &
     & beta,iseed,pgtol,dxmax,press,strext,celmas,flgkplusg,epskplusg)

!      Written June 17, 2017, from previous code. jlm
!      Modified, documentation, kplusg,August 10 2019. JLM
!      copyright inesc-mn/Jose Luis Martins/Carlos Loia Reis


!      version 4.94
!      version 1.51 of md

       use esdf

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  ipr                        !<  print control: ipr < 3 no printing, except warnings.


!      output
       
       real(REAL64), intent(out)          ::  tempk                      !<  ionic temperature (in Kelvin)
       real(REAL64), intent(out)          ::  tempinik                   !<  initial ionic temperature (in Kelvin)
       integer, intent(out)               ::  nstep                      !<  number of MD or minimization steps
       real(REAL64), intent(out)          ::  tstep                      !<  time step for MD (in atomic units)

       character(len=6), intent(out)      ::  flgkeat                    !<  adds keating force field correction to energy

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

       real(REAL64)               ::  betafrac
       integer                    ::  nlines
       integer                    ::  ioerr

!      parameters

       real(REAL64), parameter ::  ZERO = 0.0_REAL64
       real(REAL64), parameter ::  AUTOGPA = 29421.58_REAL64
       real(REAL64), parameter ::  AUTOFS = 0.02418884_REAL64

!      counters

       integer   ::  i, j

 
!      ionic temperature (in Kelvin) for langevin calculations

       tempk = 300.0
       tempk = esdf_physical('MD.TargetTemperature',tempk,'k')              

       if(tempk < ZERO) then
         tempk = ZERO
         write(6,*)
         write(6,*) '   WARNING'
         write(6,*)
         write(6,*) '   Negative temperature of ions.  Set to 0'
       endif

       if(tempk > 3000.0) then
         write(6,*)
         write(6,*) '   WARNING'
         write(6,*)
         write(6,*) '   Temperature for ions looks very high'
         write(6,'("    The value of tempk is: ",e12.3)') tempk
       endif


!      tempinik   initial temperature (in kelvin)


       tempinik = 1000.0
       tempinik = esdf_physical('MD.InitialTemperature',tempinik,'k')                 

       if(tempinik < ZERO) then
         tempinik = ZERO
         write(6,*)
         write(6,*) '   WARNING'
         write(6,*)
         write(6,*) '   Negative temperature of ions.  Set to 0'
       endif

       if(tempinik > 3000.0) then
         write(6,*)
         write(6,*) '   WARNING'
         write(6,*)
         write(6,*) '   Temperature for ions looks very high'
         write(6,'("    The value of tempinik is: ",e12.3)') tempinik
       endif


!      number of steps

       nstep = 100
       nstep = esdf_integer('MD.NumberOfSteps',nstep)
 
       if(nstep < 0) then
         write(6,*)
         write(6,*) '   WARNING'
         write(6,*)
         write(6,*) '   Number of MD steps is negative'
         write(6,*) '   The program will probably finish after the'
         write(6,*) '   initial calculations'
         write(6,'("    The value of nstep is: ",i8)') nstep
       endif

      
!      time step (in atomic units, 1 au = 2.4 10**-17 s)


!       tstep = 80.0
       tstep = 100.0
       tstep = esdf_physical('MD.LengthTimeStep',tstep*AUTOFS,'fs')
       tstep = tstep/AUTOFS

       if(tstep < 10.0) then
         write(6,*)
         write(6,*) '   WARNING'
         write(6,*)
         write(6,*) '   Time step looks very small'
         write(6,'("   The value of tstep is: ",2e12.3," fs")') tstep,   &
     &        tstep*AUTOFS
       endif

       if(tstep > 500.0) then
         write(6,*)
         write(6,*) '   WARNING'
         write(6,*)
         write(6,*) '   Time step looks very large'
         write(6,'("   The value of tstep is: ",2e12.3," fs")') tstep,   &
     &        tstep*AUTOFS
       endif


!      Friction coefficient.  It takes betafrac timesteps to dissipate
!      a reasonable energy.

       betafrac=20.0
       betafrac=esdf_double('MD.FrictionFracInvTimeStep',betafrac)

       if(betafrac < 5.0) then
         write(6,*)
         write(6,*) '   WARNING'
         write(6,*)
         write(6,*) '   Time for dissipation looks very small'
         write(6,'("   The value of betafrac is: ",e12.3)') betafrac
       endif

       if(betafrac > 200.0) then
         write(6,*)
         write(6,*) '   WARNING'
         write(6,*)
         write(6,*) '   Time for dissipation looks very large'
         write(6,'("   The value of betafrac is: ",e12.3)') betafrac
       endif

       beta = 1.0/(betafrac*tstep)
       


!      parameters for l-bfgs (kind of conjugate gradient) optimization


!       pgtol = 0.001
!       dxmax = 0.05

       pgtol = 0.0001
       dxmax = 0.01

       pgtol = esdf_physical('MD.CG.Tolerance',pgtol,'har/bohr')

       dxmax = esdf_physical('MD.CG.StepMax',dxmax,'bohr')

       if(pgtol < 1.0e-6) then
         write(6,*)
         write(6,*) '   WARNING'
         write(6,*)
         write(6,*) '   Tolerance for minimization looks very small'
         write(6,'("   The value of pgtol is: ",e12.3)') pgtol
       endif

       if(pgtol > 0.01) then
         write(6,*)
         write(6,*) '   WARNING'
         write(6,*)
         write(6,*) '   Tolerance for minimization looks very large'
         write(6,'("   The value of pgtol is: ",e12.3)') pgtol
       endif

       if(dxmax < 0.0001) then
         write(6,*)
         write(6,*) '   WARNING'
         write(6,*)
         write(6,*) '   Max step in minimization looks very small'
         write(6,'("   The value of dxmax is: ",e12.3)') dxmax
       endif

       if(dxmax > 0.5) then
         write(6,*)
         write(6,*) '   WARNING'
         write(6,*)
         write(6,*) '   Max step in minimization looks very large'
         write(6,'("   The value of dxmax is: ",e12.3)') dxmax
       endif


!      Add an ad-hoc Keating correction to tetrahedrally coordinated semicontuctors

       if(esdf_boolean('MD.UseKeatingCorrections',.false.)) then
        flgkeat='KEATNG'
       else
        flgkeat='      '       
       endif


!      initial seed for random numbers

       iseed = 87697
       iseed = esdf_integer('MD.Seed',iseed)


!      external stress and pressure (29421.58 converts from GPa to au)


       press = ZERO
       press=esdf_physical('MD.TargetPressure',press,'gpa')

       if(abs(press) > 1000.0) then
         write(6,*)
         write(6,*) '   WARNING'
         write(6,*)
         write(6,*) '   The module of pressure looks very large'
         write(6,'("   The value is: ",e12.3," GPa")') press
       endif

!      converts to atomic units

       press = press / AUTOGPA


!      the units of stress in contravariant coordinates
!      are more complicated. from the equation
!         strcont(i,j) = stress(i,j) + strkin(i,j) - strext(i,j) -
!                      press*vcell*adotm1(i,j)
!      we should be able figure out how to fill strext
!      for a cubic crystal, and for a shear component
!      we should multiply by the lattice constant and
!      divide by 29421.58 the stress in GPa

       do i=1,3
       do j=1,3
         strext(i,j) = ZERO
       enddo
       enddo

       if(esdf_block('MD.TargetStress',nlines)) then
         if(nlines /= 3) then
            write(6,*)
            write(6,*) '  Stopped in read_esdf_md  '
            write(6,*) '  Wrong Number of Lines in Stress Block!'

            stop

         endif

         ioerr = 0
         do i=1,nlines
           read(block_data(i),*,iostat=ioerr)                            &
     &              strext(1,i), strext(2,i),strext(3,i)
         
           if(ioerr /= 0) then
             write(6,*)
             write(6,*) '    Stopped in read_esdf_md'
             write(6,*) '    error reading MD.TargetStress'

             stop

           endif

         enddo
         do i = 1,3
         do j = 1,3
           strext(i,j) = strext(i,j) / AUTOGPA
         enddo
         enddo
      endif


!      fictitious cell mass for vcs md


       celmas = 10.0
       celmas=esdf_double('MD.CellMass',celmas)

       if(celmas < ZERO) then
         celmas = 10.0
         write(6,*)
         write(6,*) '   WARNING'
         write(6,*)
         write(6,*) '   Negative fictitious cell mass.  Set to 10'
       endif

       flgkplusg = .FALSE.
       flgkplusg=esdf_boolean('MD.CG.UseFixedkplusG',flgkplusg)

       epskplusg = 0.01
       epskplusg=esdf_physical('MD.CG.FixedkplusGTol',epskplusg,'har/bohr')
       if(epskplusg <= ZERO) then
         epskplusg = ZERO
         flgkplusg = .FALSE.
       endif

!      prints the results

       if(ipr > 2) then
         write(6,*)
         write(6,*)
         write(6,*)  '  The values set by read_esdf_md are:'
         write(6,*)
         write(6,'("   The value of tempk is: ",e12.3)') tempk
         write(6,'("   The value of tempinik is: ",e12.3)') tempinik
         write(6,'("   The value of nstep is: ",i8)') nstep
         write(6,'("   The value of tstep is: ",e12.3)') tstep
         write(6,'("   The value of beta is: ",2e12.3)') beta,betafrac
         write(6,'("   The value of pgtol is: ",e12.3)') pgtol
         write(6,'("   The value of dxmax is: ",e12.3)') dxmax
         write(6,*)  '  The value of flgkeat is: ',flgkeat
         write(6,'("   The value of iseed is: ",i12)') iseed
         write(6,'("   The value of press is: ",e12.3)') press
         write(6,'("   The value of strext is: ",20x,"au",33x,"GPa")')
         write(6,'(27x,3e12.3,3x,3f10.2,/,27x,3e12.3,3x,3f10.2,/,        &
     &    27x,3e12.3,3x,3f10.2)')                                        &
     &      ((strext(i,j),i=1,3),(strext(i,j)*AUTOGPA,i=1,3),j=1,3)
         write(6,'("   The value of celmas is: ",e12.3)') celmas
         write(6,'("   The value of flgkplusg is: ",l1)') flgkplusg
         write(6,'("   The value of epsklusg is: ",e12.3)') epskplusg
         write(6,*)
         write(6,*)
       endif

       return
       end subroutine read_esdf_md
