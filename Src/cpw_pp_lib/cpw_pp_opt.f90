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

!>  Driver for the calculation of the dielectric matrix
!>
!>  \author       Carlos Loia Reis, Jose Luis Martins
!>  \version      5.08
!>  \date         20 January 2022.
!>  \copyright    GNU Public License v2

subroutine cpw_pp_opt(ioreplay,                                          &
           dims_, flags_, crys_, recip_, spaceg_, pseudo_, atorb_,       &
           pwexp_, strfac_,  vcomp_,                                     &
           efermi, meta_cpw2000, title, subtitle,                        &
           epspsi, icmax)

! Breakup of cpw_pp_band_dos_opt. 20 Janeiro 2022. JLM
! But tracked from Pedro Borlido report.  23 October 2023. JLM

  use cpw_variables

  implicit none

!  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, intent(in)                     ::  ioreplay                   !<  tape number for reproducing calculations

  type(dims_t)                            ::  dims_                      !<  array dimensions
  type(flags_t)                           ::  flags_                     !<  computational flags
  type(crys_t)                            ::  crys_                      !<  crystal structure
  type(recip_t)                           ::  recip_                     !<  reciprocal space information
  type(spaceg_t)                          ::  spaceg_                    !<  space group information
  type(pseudo_t)                          ::  pseudo_                    !<  pseudo-potential (Kleinman-Bylander)
  type(atorb_t)                           ::  atorb_                     !<  atomic orbitals in G-space
  type(pwexp_t)                           ::  pwexp_                     !<  plane-wave expansion choices
  type(strfac_t)                          ::  strfac_                    !<  structure factors
  type(vcomp_t)                           ::  vcomp_                     !<  local potential contributions

 ! other input


  real(REAL64), intent(in)                ::  efermi                     !<  eigenvalue of highest occupied state (T=0) or fermi energy (T/=0), Hartree

  character(len=250), intent(in)          ::  meta_cpw2000               !<  metadata from cpw2000
  character(len=50), intent(in)           ::  title                      !<  title for plots
  character(len=140), intent(in)          ::  subtitle                   !<  subtitle for plots

  real(real64), intent(in)                :: epspsi                      !<  accuracy of eigenvalues
  integer, intent(in)                     :: icmax                       !<  maximum number of iterations for diagonalization

! other variables

  integer                 ::  ios
  integer                 ::  imeth, idiag


  integer                 :: iguess                                      !  kept for compatibility

  character(len=1)        ::  yesno

  character(len=4)        ::  diag_type                                  !  selects diagonalization, 'pw  ','ao  ','aojc'
  logical                 ::  lworkers                                   !  use workers in calculation

  real(real64)            ::  xsvd, csvd

  character(len=12)       ::  file_dos_mesh                              !  SHOULD BE INPUT AND PROPAGATED

! constants

  real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter  :: HARTREE = 27.21138386_REAL64

  file_dos_mesh = 'DOS_MESH.DAT'

   write(6,*)
   write(6,*)  '  The mesh will be read from file  ', file_dos_mesh
   write(6,*)
   write(6,*)  '  Which method you want to use?'
   write(6,*)  '  1:  full plane wave basis diagonalization'
   write(6,*)  '  2:  diagonalization in a Luttinger-Kohn basis'

   read(5,*,iostat=ios) imeth
   write(ioreplay,*) imeth,'   method'

   if(ios /= 0) then
     imeth = 0
     write(6,*)
     write(6,*) '  error processing default input '
     write(6,*)

     return

   endif

   if(imeth == 1) then

     write(6,*)
     write(6,*)  '  Are you a human (y/n)? '
     write(6,*)  '  You need to be very careful if you answer no!'

     read(5,*) yesno
     write(ioreplay,*) yesno,'   Blade Runner test'
     if(yesno == 'N' .or. yesno == 'n') then
       lworkers = .TRUE.
     else
       lworkers = .FALSE.
     endif

     diag_type = 'pw  '

     if(atorb_%latorb) then
       write(6,*)
       write(6,*)  '  which diagonalization do you want to use'
       write(6,*)  '  1)  full pw diagonalization'
       write(6,*)  '  2)  diagonalization in atomic orbitals followed by jacobian relaxation'
       write(6,*)  '  3)  diagonalization in atomic orbitals'

       read(5,*,iostat=ios) idiag
       write(ioreplay,*) idiag,'   diagonal. method in DOS'

       if(idiag == 2) then
         diag_type = 'aojc'
         write(6,*) '  using aojc diagonalization'
       elseif(idiag == 3) then
         diag_type = 'ao  '
         write(6,*) '  using ao diagonalization'
       else
         write(6,*) '  using pw diagonalization'
       endif
     endif

     call out_opt(diag_type, lworkers,                                   &
     title, subtitle,                                                    &
     pwexp_%emax, flags_%flgdal, flags_%flgpsd,                          &
     iguess, epspsi, icmax, pseudo_%ztot,                                &
     crys_%adot, crys_%ntype, crys_%natom, crys_%rat,                    &
     spaceg_%ntrans, spaceg_%mtrx,                                       &
     recip_%ng, recip_%kgv,recip_%phase,recip_%conj,                     &
     recip_%ns, recip_%inds, recip_%kmax,                                &
     recip_%indv, recip_%ek,                                             &
     strfac_%sfact, strfac_%icmplx,                                      &
     vcomp_%veff,                                                        &
     pseudo_%nq, pseudo_%delq, pseudo_%vkb, pseudo_%nkb,                 &
     atorb_%latorb, atorb_%norbat,atorb_%nqwf,atorb_%delqwf,             &
     atorb_%wvfao,atorb_%lorb,                                           &
     dims_%mxdtyp, dims_%mxdatm, dims_%mxdgve, dims_%mxdnst,             &
     dims_%mxdlqp, dims_%mxdcub, dims_%mxdlao)

  elseif(imeth == 2) then

     write(6,*)
     write(6,*)  '  Using Generalized Luttinger-Kohn interpolation'
     write(6,*)

     write(6,*)
     write(6,*)  '  Are you a human (y/n)? '
     write(6,*)  '  You need to be very careful if you answer no!'

     read(5,*) yesno
     write(ioreplay,*) yesno,'   Blade Runner test'
     if(yesno == 'N' .or. yesno == 'n') then
       lworkers = .TRUE.
     else
       lworkers = .FALSE.
     endif

     write(6,*)
     write(6,*) '  Do you want to use the default GLK values'
     write(6,*)

     read(5,*) yesno
     write(ioreplay,*) yesno,'   default GLK values'

     if(yesno == 'N' .or. yesno == 'n') then

       diag_type = 'pw  '

       if(atorb_%latorb) then
         write(6,*)
         write(6,*)  '  which diagonalization do you want to use'
         write(6,*)  '  1)  full pw diagonalization'
         write(6,*)  '  2)  diagonalization in atomic orbitals followed by jacobian relaxation'
         write(6,*)  '  3)  diagonalization in atomic orbitals'

         read(5,*,iostat=ios) idiag
         write(ioreplay,*) idiag,'   diagonal. method in LK DOS'

         if(idiag == 2) then
           diag_type = 'aojc'
           write(6,*) '  using aojc diagonalization'
         elseif(idiag == 3) then
           diag_type = 'ao  '
           write(6,*) '  using ao diagonalization'
         else
           write(6,*) '  using pw diagonalization'
         endif
       endif


       write(6,*)
       write(6,*) '  Maximum number of bands to keep after SVD?'
       write(6,*) '  Multiple of nband, answer between 1 and 4.'
       write(6,*)

       read(5,*,iostat=ios) csvd
       write(ioreplay,*) csvd,'   csvd in GLK'
       if(ios /= 0 .or. csvd < UM) then
         csvd = UM
         write(6,*)
         write(6,'("  csvd = ",f14.6)') csvd
         write(6,*)
       endif
       if(csvd > 4*UM) then
         csvd = 4*UM
         write(6,*)
         write(6,'("  csvd = ",f14.6)') csvd
         write(6,*)
       endif

       write(6,*)
       write(6,*) '  Singular value threshold to keep bands after SVD?'
       write(6,*) '  Answer between 0 and 0.5.'
       write(6,*)

       read(5,*,iostat=ios) xsvd
       write(ioreplay,*) xsvd,'   xsvd in GLK'
       if(ios /= 0 .or. xsvd < 0.00001) then
         xsvd = 0.00001
         write(6,*)
         write(6,'("  xsvd = ",f14.6)') xsvd
         write(6,*)
       endif
       if(xsvd > 0.5) then
         xsvd = 0.5
         write(6,*)
         write(6,'("  xsvd = ",f14.6)') xsvd
         write(6,*)
       endif

     else

       csvd = 1.5
       xsvd = 0.00001

       diag_type = 'pw  '

     endif


     call out_opt_glk(diag_type, lworkers, xsvd, csvd,                   &
     title, subtitle,                                                    &
     pwexp_%emax, flags_%flgdal, flags_%flgpsd,                          &
     iguess, epspsi, icmax, pseudo_%ztot,                                &
     crys_%adot, crys_%ntype, crys_%natom, crys_%rat,                    &
     spaceg_%ntrans, spaceg_%mtrx,                                       &
     recip_%ng, recip_%kgv,recip_%phase,recip_%conj,                     &
     recip_%ns, recip_%inds, recip_%kmax,                                &
     recip_%indv, recip_%ek,                                             &
     strfac_%sfact, strfac_%icmplx,                                      &
     vcomp_%veff,                                                        &
     pseudo_%nq, pseudo_%delq, pseudo_%vkb, pseudo_%nkb,                 &
     atorb_%latorb, atorb_%norbat, atorb_%nqwf, atorb_%delqwf,           &
     atorb_%wvfao,atorb_%lorb,                                           &
     dims_%mxdtyp, dims_%mxdatm, dims_%mxdgve, dims_%mxdnst,             &
     dims_%mxdlqp, dims_%mxdcub, dims_%mxdlao)


  endif

  return
end subroutine cpw_pp_opt

