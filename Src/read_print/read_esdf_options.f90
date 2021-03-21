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

!>     reads the computational options from a file
!>     opened by a previous edsf_init

       subroutine read_esdf_options(flgcal,ipr,                          &
     & flgscf,flgdal,flgmix,teleck,                                      &
     & itmax,epscv,epscvao,epspsi,symkip,symtol)

!      Written June 17, 2017, from previous code. jlm
!      Modified, documentation, August 10 2019. JLM
!      Modified, symkip for EPILBF, October 13, 2019. JLM
!      copyright inesc-mn/Jose Luis Martins/Carlos Loia Reis


!      version 4.94

       use esdf

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)


!      input

       character(len=6), intent(in)       ::  flgcal                     !<  type of calculation

       integer, intent(in)                ::  ipr                        !<  print control: ipr < 3 no printing, except warnings.

!      output

       character(len=6), intent(out)      ::  flgscf                     !<  type of self consistent field and diagonalization
       character(len=4), intent(out)      ::  flgdal                     !<  whether the dual approximation is used
       character(len=6), intent(out)      ::  flgmix                     !<  choice of potential mixing
       
       real(REAL64), intent(out)          ::  teleck                     !<  electronic temperature (in Kelvin)

       integer, intent(out)               ::  itmax                      !<  maximum number of self consistency cycles
       real(REAL64), intent(out)          ::  epscv                      !<  convergence criteria for potential self consistency
       real(REAL64), intent(out)          ::  epscvao                    !<  convergence criteria for potential self consistency in atomic orbitals
       real(REAL64), intent(out)          ::  epspsi                     !<  convergence criteria for iterative diagonalization
       
       logical, intent(out)               ::  symkip                     !<  whether symmetry should be conserved
       real(REAL64), intent(out)          ::  symtol                     !<  tolerance for symmetry recognition subroutines

!      parameters

       real(REAL64), parameter ::  ZERO = 0.0_REAL64


!      type of self-consistent field calculation

!      PW     : old style calculation using only plane wave basis. It has therefore "full" accuracy defined by emax
!      AO     : it uses an atomic basis orbital (although integrals are calculated with plane-waves). limited accuracy, very fast.
!      AOJC   : the results of the atomic orbital basis are improved with a single jacobian relaxation. reasonable accuracy, fast.
!      AOJCPW : atomic orbital basis are used to obtain a better starting point for normal calculation. "full" accuracy, speed gain.

       flgscf = '    PW'
!       flgscf = 'AO    '
!       flgscf = 'AOJC  '
!       flgscf = 'AOJCPW'

       flgscf=esdf_string('TypeOfScfDiag',flgscf)
       call chrcap(flgscf,6)
       if(flgscf == 'PW    ') then
            flgscf ='    PW' 
       endif


!      decides if the "dual" aproximation is used. as usual the price for a faster
!      calculation is a smaller accuracy. pressure is a result that is sensitive
!      to this approximation. if pressure is accurate enough, all other
!      results should be fine.

!       flgdal = '    '
       flgdal = 'DUAL'

       if(esdf_boolean('DualApproximation',.true.)) then
            flgdal='DUAL'
       else 
            flgdal = '    '
       endif


!      chooses between the type of mixing. which is best may depend
!      on system, but bfgs is usually better, and is the only currently implemented.

!       flgmix = 'BROYD1'
       flgmix = 'BFGS  '

       flgmix=esdf_string('TypeOfPseudoMixing',flgmix)
       call chrcap(flgmix,6)


!      maximum number of self consistency cycles. Increase to whatever needed

       itmax = 30
       itmax = esdf_integer('MaxSCFIterations',itmax)              

       if(itmax < 1) then
         itmax = 1
         write(6,*)
         write(6,*) '   WARNING'
         write(6,*)
         write(6,*) '   Number of SCF iterations set to 1'
       endif

       if(itmax > 100) then
         write(6,*)
         write(6,*) '   WARNING'
         write(6,*)
         write(6,'("    The value of itmax seems too large: ",i6)')   &
     &            itmax
       endif

       if(itmax < 5) then
         write(6,*)
         write(6,*) '   WARNING'
         write(6,*)
         write(6,'("    The value of itmax seems too small: ",i6)')   &
     &            itmax
       endif


!      convergence criteria for self consistency. the value of 0.00001 is
!      a safe bet. using 0.00005 saves some time, but you should check
!      if it is accurate enough.

       epscv = 0.00001
!       epscv = 0.0005
       epscv= esdf_double('ScfTolerance',epscv)

       if(epscv > 0.001) then
         write(6,*)
         write(6,*) '   WARNING'
         write(6,*)
         write(6,'("    The value of epscv seems too large: ",e12.3)')   &
     &            epscv
       endif

       if(epscv < 1.0e-7) then
         write(6,*)
         write(6,*) '   WARNING'
         write(6,*)
         write(6,'("    The value of epscv seems too small: ",e12.3)')   &
     &            epscv
       endif

!      convergence criteria for self consistency with atomic orbitals

       if(flgscf == 'AOJCPW') then
         epscvao = 20.0*epscv
       else
         epscvao = epscv
       endif


!      Convergence criteria for diagonalization.  The value of 0.00001 is
!      a safe bet.  Using 0.0001 saves some time, but you should check
!      if it is accurate enough.


       epspsi = 0.00001
!       epspsi = 0.001

       epspsi= esdf_double('DiagTolerance',epspsi)       

       if(epspsi > 0.01) then
         write(6,*)
         write(6,*) '   WARNING'
         write(6,*)
         write(6,'("    The value of epspsi seems too large: ",e12.3)')  &
     &            epspsi
       endif

       if(epspsi < 1.0e-8) then
         write(6,*)
         write(6,*) '   WARNING'
         write(6,*)
         write(6,'("    The value of epspsi seems too small: ",e12.3)')  &
     &            epspsi
       endif


!      electronic temperature (in Kelvin). Use zero, or ionic temperature in most cases.


       teleck = ZERO
!       teleck = tempk
       teleck = esdf_physical('ElectronicTemperature',teleck,'k')              

       if(teleck < ZERO) then
         teleck = ZERO
         write(6,*)
         write(6,*) '   WARNING'
         write(6,*)
         write(6,*) '   Negative temperature of electrons.  Set to 0'
       endif

       if(teleck > 5000.0) then
         write(6,*)
         write(6,*) '   WARNING'
         write(6,*)
         write(6,*) '   Temperature for electrons looks very high'
         write(6,'("    The value of teleck is: ",e12.3)') teleck
       endif


!      conservation of symmetry, and tolerance of symmetry recognition

       symtol = 1.0D-5
       symtol = esdf_double('SymmTolerance',symtol)

       if(symtol > 0.01) then
         write(6,*)
         write(6,*) '   WARNING'
         write(6,*)
         write(6,'("    The value of symtol seems too large: ",e12.3)')  &
     &            symtol
       endif

       if(symtol < 1.0e-8) then
         write(6,*)
         write(6,*) '   WARNING'
         write(6,*)
         write(6,'("    The value of symtol seems too small: ",e12.3)')  &
     &            symtol
       endif

       symkip = .FALSE.
       if(flgcal == 'ONE   ') symkip = .TRUE.
       if(flgcal == 'ONEVRD') symkip = .TRUE.
       if(flgcal == 'LBFSYM') symkip = .TRUE.
       if(flgcal == 'VCSLBF') symkip = .TRUE.
       if(flgcal == 'EPILBF') symkip = .TRUE.

       symkip = esdf_boolean('UseSymmetry',symkip)


!      prints the results

       if(ipr > 2) then
         write(6,*)
         write(6,*)
         write(6,*)  '  The values set by read_esdf_options are:'
         write(6,*)
         write(6,*)  '  The value of flgscf is: ',flgscf
         write(6,*)  '  The value of flgdal is: ',flgdal
         write(6,*)  '  The value of flgmix is: ',flgmix
         write(6,'("   The value of itmax is: ",i6)') itmax
         write(6,'("   The value of epscv is: ",e12.3)') epscv
         write(6,'("   The value of epscvao is: ",e12.3)') epscvao
         write(6,'("   The value of epspsi is: ",e12.3)') epspsi
         write(6,'("   The value of teleck is: ",e12.3)') teleck
         write(6,'("   The value of symtol is: ",e12.3)') symtol
         write(6,*)  '  The value of symkip is: ',symkip
         write(6,*)
         write(6,*)
       endif

       return
       end subroutine read_esdf_options
