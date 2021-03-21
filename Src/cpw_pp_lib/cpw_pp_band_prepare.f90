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

!>  Calculates several quantities related to the effective potential
!>  that are neede to calculate the bands at a given k-point.

subroutine cpw_pp_band_prepare(ioreplay,                                 &
    dims_, crys_, spaceg_, recip_, recip_in_, pwexp_, strfac_,           &
    vcomp_, vcompin_,                                                    &
    emaxin)

! written February 2, 2020 from previous code. JLM
! Modified, consistent space group, mstar bug. 17 January 2021. JLM
! Modified, bug in initialization of chd. 12 February 2021. JLM
! copyright  Jose Luis Martins/INESC-MN

! version 4.99

  use cpw_variables

  implicit none

!  integer, parameter          :: REAL64 = selected_real_kind(12)

  type(dims_t)                       ::  dims_                           !<  array dimensions
  type(crys_t)                       ::  crys_                           !<  crystal structure
  type(vcomp_t)                      ::  vcompin_                        !<  local potential contributions
  type(recip_t)                      ::  recip_in_                       !<  reciprocal space information
  type(recip_t)                      ::  recip_                          !<  reciprocal space information
  type(spaceg_t)                     ::  spaceg_                         !<  space group information
  type(pwexp_t)                      ::  pwexp_                          !<  plane-wave expansion choices
  type(strfac_t)                     ::  strfac_                         !<  structure factors
  type(vcomp_t)                      ::  vcomp_                          !<  local potential contributions

! input

  integer, intent(in)                :: ioreplay                         !<  tape number for reproducing calculations

  real(REAL64), intent(in)           ::  emaxin                          !<  kinetic energy cutoff of plane wave expansion (hartree).

! other variables

  integer           ::  ipr
  integer           ::  isym

  integer           ::  ngmax
  integer           ::  istatus                                          !  istatus = 0, successful; 1 not closed; 2 no inverse; 3 inconsistent with atomic positions

! allocatable local arrays       

  complex(REAL64), allocatable       ::  chd(:,:,:)                      !  effective potential on the grid

! counters
  integer    ::  i,j,k

! constants

  real(REAL64), parameter    :: ZERO = 0.0_REAL64
  real(REAL64), parameter    :: TOL = 1.0E-07_REAL64
  complex(REAL64), parameter :: C_ZERO = cmplx(ZERO,ZERO,REAL64)


  
  isym = 1
  ipr = 1
  if(spaceg_%ntrans == 0) then
    call sym_identify(isym, ipr, TOL,                                    &
       spaceg_%ntrans, spaceg_%mtrx, spaceg_%tnp,                        &
       crys_%adot, crys_%ntype, crys_%natom, crys_%rat,                  &
       dims_%mxdtyp, dims_%mxdatm)
  else
    ipr = 0
    call sym_test(ipr, TOL, istatus,                                     &
       spaceg_%ntrans, spaceg_%mtrx, spaceg_%tnp,                        &
       crys_%ntype, crys_%natom, crys_%rat, crys_%adot,                  &
       dims_%mxdtyp, dims_%mxdatm)
    if(istatus /= 0) then
      write(6,*) '  STOPPED in cpw_band_prepare'
      write(6,*) '  inconsistent space group, istatus = ',istatus

      stop

    endif
  endif

  write(6,*)
  write(6,'("  The original calculation used a maximum energy",          &
     " PW cutoff of",f10.3," Hartree")') emaxin
  write(6,*) '  Enter maximum energy in Hartree '
  read(5,*) pwexp_%emax
  write(ioreplay,*) pwexp_%emax,'   emax'

  call size_g_space(pwexp_%emax, crys_%adot,                             &
     spaceg_%ntrans, spaceg_%mtrx,                                       &
     dims_%mxdgve, dims_%mxdnst, dims_%mxdcub)

  allocate(recip_%kgv(3,dims_%mxdgve))
  allocate(recip_%phase(dims_%mxdgve))
  allocate(recip_%conj(dims_%mxdgve))
  allocate(recip_%inds(dims_%mxdgve))
  allocate(recip_%indv(dims_%mxdcub))
  allocate(recip_%mstar(dims_%mxdnst))
  allocate(recip_%ek(dims_%mxdnst))

  allocate(recip_%izstar(dims_%mxdnst))
  
  ipr = 1
  call g_space(ipr, pwexp_%emax,                                         &
      crys_%adot, spaceg_%ntrans, spaceg_%mtrx, spaceg_%tnp,             &
      recip_%ng, recip_%kgv, recip_%phase, recip_%conj,                  &
      recip_%inds, recip_%kmax, recip_%indv, recip_%ns, recip_%mstar,    &
      recip_%ek, recip_%izstar,                                          &
      dims_%mxdgve, dims_%mxdnst, dims_%mxdcub)
  
  deallocate(recip_%izstar)

  allocate(strfac_%sfact(dims_%mxdtyp,dims_%mxdnst))

  call structure_factor(strfac_%sfact, strfac_%icmplx,                   &
      recip_%ng, recip_%kgv, recip_%phase, recip_%conj, recip_%inds,     &
      crys_%ntype, crys_%natom, crys_%rat,                               &
      dims_%mxdtyp, dims_%mxdatm, dims_%mxdgve, dims_%mxdnst)

  allocate(chd(-recip_%kmax(1):recip_%kmax(1),                           &
               -recip_%kmax(2):recip_%kmax(2),                           &
               -recip_%kmax(3):recip_%kmax(3)))
  allocate(vcomp_%veff(dims_%mxdnst))
  
  do k = -recip_%kmax(3),recip_%kmax(3)
  do j = -recip_%kmax(2),recip_%kmax(2)
  do i = -recip_%kmax(1),recip_%kmax(1)
    chd(i,j,k) = C_ZERO
  enddo
  enddo
  enddo

  ngmax = 0
  do i=1,min(recip_%ns,recip_in_%ns)
    do j=ngmax+1,ngmax+recip_in_%mstar(i)
      if(abs(recip_in_%kgv(1,j)) <= recip_%kmax(1) .and.                 &
         abs(recip_in_%kgv(2,j)) <= recip_%kmax(2) .and.                 &
         abs(recip_in_%kgv(3,j)) <= recip_%kmax(3) ) then

        if(recip_in_%conj(j) > ZERO) then
          chd(recip_in_%kgv(1,j),recip_in_%kgv(2,j),recip_in_%kgv(3,j)) =      &
                        vcompin_%veff(i)*conjg(recip_in_%phase(j))
        else
          chd(recip_in_%kgv(1,j),recip_in_%kgv(2,j),recip_in_%kgv(3,j)) =      &
                        conjg(vcompin_%veff(i))*recip_in_%phase(j)
        endif

      endif
    enddo
    ngmax = ngmax + recip_in_%mstar(i)
    
  enddo
  

! collects v_effective in stars

  call cube_to_star(vcomp_%veff, recip_%kmax,chd,                        &
      recip_%ng, recip_%kgv, recip_%phase, recip_%conj, recip_%ns,       &
      recip_%mstar,                                                      &
      dims_%mxdgve, dims_%mxdnst)

  return
end subroutine cpw_pp_band_prepare

