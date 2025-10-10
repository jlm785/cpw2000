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

!>  Reads data for the calculation.
!>  Interfaces between variables in modules and real ones.
!>
!>  \author       Jose Luis Martins and Carlos Loia Reis
!>  \version      5.12
!>  \date         November 2019, 10 October 2025.
!>  \copyright    GNU Public License v2

subroutine cpw_read_data(fcpwin, fpwdat, iopw, vdriv,                    &
     iprglob, meta_cpw2000, meta_pwdat, symkip, symtol,                  &
     flags_, crys_, pwexp_, kpoint_, xc_, acc_, moldyn_, vcsdyn_,        &
     total_, ewald_, filename_, dims_)

! Adapted November 2019 from cpw.f90. JLM
! Modified, metadata bug, May 2020. JLM, CLR
! Modified, close unit, 11 June 2020. JLM
! Modified, icdiagmax, indentation, 14 May 2025. JLM
! Modified, filenames for pseudos and psi to disk. 10 October 2025. JLM

  use cpw_variables

  implicit none

  type(flags_t)                      ::  flags_                          !<  computational flags
  type(crys_t)                       ::  crys_                           !<  crystal structure
  type(pwexp_t)                      ::  pwexp_                          !<  plane-wave expansion choices
  type(kpoint_t)                     ::  kpoint_                         !<  k-point data
  type(xc_t)                         ::  xc_                             !<  exchange and correlation choice
  type(acc_t)                        ::  acc_                            !<  accuracy parameters
  type(moldyn_t)                     ::  moldyn_                         !<  molecular dynamics variables
  type(vcsdyn_t)                     ::  vcsdyn_                         !<  variational cell shape molecular dynamics variables
  type(dims_t)                       ::  dims_                           !<  array dimensions
  type(enfrst_t)                     ::  total_                          !<  Total energy force stress
  type(enfrst_t)                     ::  ewald_                          !<  Ewald energy force stress
  type(filename_t)                   ::  filename_                       !<  Information about used files

  character(len=*), intent(in)       ::  fcpwin                          !<  filename with input data (new style)
  character(len=*), intent(in)       ::  fpwdat                          !<  filename with input data (old style)

  character(len=4), intent(in)       ::  vdriv                           !<  version of the calling program
  integer, intent(in)                ::  iopw                            !<  tape number

  integer, intent(out)               ::  iprglob                         !<  controls the amount of printing by subroutines

  character(len=250), intent(out)    ::  meta_cpw2000                    !<  metadata coded in cpw2000
  character(len=250), intent(out)    ::  meta_pwdat                      !<  metadata from pw.dat

  logical, intent(out)               ::  symkip                          !<  whether symmetry should be conserved
  real(REAL64), intent(out)          ::  symtol                          !<  tolerance for symmetry recognition subroutines

  integer      ::  ioerr
  logical               ::  lgeom                                        !  indicates if geometry was successfully read.
  logical               ::  lbz                                          !  indicates if Brillouin Zone data was successfully read.
  integer      ::  ipr

  integer    ::    i

  ipr = 3
  call size_mxdtyp_mxdatm_esdf(ipr, fcpwin, lgeom,                       &
          dims_%mxdtyp, dims_%mxdatm)

  do i = 1,250
    meta_pwdat(i:i) = ' '
    meta_cpw2000(i:i) = ' '
  enddo

  if(.NOT. lgeom) then
    call size_mxdtyp_mxdatm(fpwdat,iopw,dims_%mxdtyp,dims_%mxdatm)
  endif

  allocate(crys_%natom(dims_%mxdtyp))
  allocate(crys_%rat(3,dims_%mxdatm,dims_%mxdtyp))
  allocate(crys_%atmass(dims_%mxdtyp))
  allocate(crys_%nameat(dims_%mxdtyp))

  allocate(moldyn_%vat(3,dims_%mxdatm,dims_%mxdtyp))
  allocate(moldyn_%rat1(3,dims_%mxdatm,dims_%mxdtyp))
  allocate(moldyn_%frc1(3,dims_%mxdatm,dims_%mxdtyp))

  allocate(total_%force(3,dims_%mxdatm,dims_%mxdtyp))
  allocate(ewald_%force(3,dims_%mxdatm,dims_%mxdtyp))


  call read_esdf(fcpwin, vdriv,                                          &
      flags_%flgcal, flags_%flgkeat, flags_%flgpsd, flags_%flgscf,       &
          flags_%flgdal, flags_%flgmix,                                  &
      crys_%adot, crys_%ntype, crys_%natom, crys_%nameat, crys_%rat,     &
          crys_%atmass, crys_%alatt,                                     &
      lgeom,                                                             &
      pwexp_%emax, pwexp_%nbandin,                                       &
      kpoint_%nx, kpoint_%ny, kpoint_%nz, kpoint_%sx, kpoint_%sy,        &
          kpoint_%sz,                                                    &
          lbz,                                                           &
      meta_cpw2000,                                                      &
      symkip, symtol,                                                    &
      xc_%author, xc_%tblaha,                                            &
      iprglob,                                                           &
      acc_%itmax, acc_%icdiagmax, acc_%epscv, acc_%epscvao, acc_%epspsi, &
      moldyn_%tempk,                                                     &
      pwexp_%teleck,                                                     &
      moldyn_%tempinik, moldyn_%nstep, moldyn_%tstep, moldyn_%beta,      &
      moldyn_%iseed,                                                     &
      moldyn_%pgtol, moldyn_%dxmax,                                      &
      vcsdyn_%press, vcsdyn_%strext, vcsdyn_%celmas,                     &
      pwexp_%lkplusg, pwexp_%epskplusg,                                  &
      filename_%pseudo_path, filename_%pseudo_suffix,                    &
      filename_%itape_pseudo,                                            &
      filename_%save_psi_path, filename_%itape_save_psi,                 &
      dims_%mxdtyp, dims_%mxdatm)


  meta_pwdat = meta_cpw2000


! prints type of calculation near the top of output file

  call tpage_calc(flags_%flgcal)

  if(.NOT. lgeom) then

!   compatibility with some old versions

    write(6,*)
    write(6,*) ' WARNING '
    write(6,*) ' Using ',fpwdat,' file'
    write(6,*)

    open(UNIT=iopw, FILE=fpwdat, STATUS='OLD', FORM='FORMATTED', IOSTAT=ioerr)

    if(ioerr /= 0) then
      write(6,*)
      write(6,*)
      write(6,*) '  cpw:   Unable to open PW.DAT file ',fpwdat

      stop

    endif

    read(iopw,'(a250)', iostat=ioerr) meta_pwdat
    if(ioerr == 0) then
      write(6,*) meta_pwdat(1:60)
      write(6,*) meta_pwdat(61:110)
      write(6,*) meta_pwdat(111:250)
    else
      write(6,*) '  cpw:   Unable to read titles'
    endif
    write(6,*)
    write(6,*)

!   gets crystal data

    call read_data(crys_%adot, crys_%ntype, crys_%natom,                 &
       crys_%nameat, crys_%rat, crys_%atmass, crys_%alatt,               &
       dims_%mxdtyp, dims_%mxdatm)

    if(.NOT. lbz) then
      read(iopw,'(1x,f9.4)') pwexp_%emax

      read(iopw, '(i5, 15x, 3i5, 5x, 3f5.2)') pwexp_%nbandin,            &
        kpoint_%nx, kpoint_%ny, kpoint_%nz,                              &
        kpoint_%sx, kpoint_%sy, kpoint_%sz

    endif

    close(UNIT=iopw)

!   copies meta_pwdat to meta_cpw2000

    meta_cpw2000 = meta_pwdat


  endif

  return

end subroutine cpw_read_data
