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

!>  Computes v hartree and adds v exchange-correlation.

subroutine v_hartree_xc(ipr, author, tblaha, adot, exc, strxc, rhovxc,   &
    vhar, vxc, den, denc, rholap, tau,                                      &
    ng, kgv, phase, conj, ns, inds, kmax, mstar, ek,                           &
    mxdgve, mxdnst)

! Adapted from Sverre Froyen plane wave program
! written june 8 1987. jlm
! modified february 18 1990. jlm
! modified 4.0 16 october 1993. jlm
! modified (conj) 20-22 march 1999. jlm
! modified (chd) 24 july 2002. jlm
! meta-GGA added. CLR
! modified, f90, meta-GGA, etc. 30 September 2015.  JLM
! modified, rholap, tau, 7 October 2015. JLM
! Modified documentation, January 2020. JLM
! Modified ipr, icheck, 13 February 2021. JLM
! copyright INESC-MN/Jose Luis Martins/Carlos Loia Reis

! Version 4.99


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdnst                          !<  array dimension for g-space stars

  integer, intent(in)                ::  ipr                             !<  print switch
  character(len=*), intent(in)       ::  author                          !<  type of xc wanted (ca=pz , pbe, tbl)
  real(REAL64), intent(in)           ::  tblaha                          !<  Tran-Blaha constant
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space

  complex(REAL64), intent(in)        ::  den(mxdnst)                     !<  density for the prototype G-vector
  complex(REAL64), intent(in)        ::  denc(mxdnst)                    !<  core density for the prototype G-vector
  complex(REAL64), intent(in)        ::  rholap(mxdnst)
  complex(REAL64), intent(in)        ::  tau(mxdnst)

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
  complex(REAL64), intent(in)        ::  phase(mxdgve)                   !<  phase factor of G-vector n
  real(REAL64), intent(in)           ::  conj(mxdgve)                    !<  is -1 if one must take the complex conjugate of x*phase
  integer, intent(in)                ::  ns                              !<  number os stars with length less than gmax
  integer, intent(in)                ::  inds(mxdgve)                    !<  star to which g-vector n belongs
  integer, intent(in)                ::  kmax(3)                         !<  max value of |kgv(i,n)|
  integer, intent(in)                ::  mstar(mxdnst)                   !<  number of g-vectors in the j-th star
  real(REAL64), intent(in)           ::  ek(mxdnst)                      !<  kinetic energy (hartree) of g-vectors in star j


! output

  real(REAL64), intent(out)          ::  exc                             !<  total XC energy (int rho * e_XC in Hartree).
  real(REAL64), intent(out)          ::  rhovxc                          !<  correction for XC energy (int rho * v_XC in Hartree).
  real(REAL64), intent(out)          ::  strxc(3,3)                      !<  contribution of xc to the stress tensor (contravariant,Hartree)
  complex(REAL64), intent(out)       ::  vhar(mxdnst)                    !<  Hartree potential (in Hartree) for the prototype G-vector
  complex(REAL64), intent(out)       ::  vxc(mxdnst)                     !<  XC potential (in Hartree) for the prototype G-vector

! local allocatable arrays

  real(REAL64), allocatable          ::  vxcmsh(:)
  real(REAL64), allocatable          ::  rhomsh(:)
  real(REAL64), allocatable          ::  rholapmsh(:)
  real(REAL64), allocatable          ::  taumsh(:)
  real(REAL64), allocatable          ::  wrkfft(:)
  complex(REAL64), allocatable       ::  chd(:)
  complex(REAL64), allocatable       ::  dentot(:)
  complex(REAL64), allocatable       ::  vxcg(:)

! local variables

  real(REAL64)    ::  vcell, bdot(3,3)
  integer         ::  mxdwrk, mxdfft
  integer         ::  nsfft(3), n1,n2,n3, id, ntot
  real(REAL64)    ::  fac
  integer         ::  ncheck(4)

! counters

  integer       ::  i

! parameters

  real(REAL64), parameter     ::  PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter     ::  ZERO = 0.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)


  call adot_to_bdot(adot, vcell, bdot)

! printout charge density

  if (ipr > 1) then
    write(6,*)
    write(6,*) 'charge density'
    write(6,*)
    do i = 1,min(ns,20)
      write(6,'(i8,4f14.6)') i,den(i),denc(i)
    enddo
  endif

  call size_fft(kmax, nsfft, mxdfft, mxdwrk)

  n1 = nsfft(1)
  n2 = nsfft(2)
  n3 = nsfft(3)
  id = nsfft(1) + 1
  ntot = id * n2 * n3

  allocate(dentot(mxdnst))
  allocate(rhomsh(mxdfft))
  allocate(vxcmsh(mxdfft))
  allocate(chd(mxdfft))
  allocate(wrkfft(mxdwrk))

  do i = 1,ns
    dentot(i) = den(i) + denc(i)
  enddo

  call mesh_set(ipr, 'v_hartree_xc', adot, dentot, rhomsh, ncheck,       &
      ng, kgv, phase, conj, inds, kmax,                                  &
      mxdgve, mxdnst, mxdfft)

  if(ncheck(1) /= n1 .or. ncheck(2) /= n2 .or. ncheck(3) /= n3           &
      .or. ncheck(4) /= id) then
    write(6,*)
    write(6,*)  "  STOPPED in v_hartree_xc:  inconsistency in mesh_set:"
    write(6,'("  in v_hxc ",4i6,"  in mesh_set ",4i6)') n1,n2,n3,id,     &
               (ncheck(i),i=1,4)
    write(6,*)

    stop

  endif

!-----------------------clr----------------------------------------

! compute exchange and correlation

  if(author=="TBL") then

    allocate(rholapmsh(mxdfft))
    allocate(taumsh(mxdfft))

    call mesh_set(ipr, "rholap", adot, rholap, rholapmsh, ncheck,        &
        ng, kgv, phase, conj, inds, kmax,                                &
        mxdgve, mxdnst, mxdfft)

  if(ncheck(1) /= n1 .or. ncheck(2) /= n2 .or. ncheck(3) /= n3           &
      .or. ncheck(4) /= id) then
    write(6,*)
    write(6,*)  "  STOPPED in v_hartree_xc:  inconsistency in mesh_set:"
    write(6,'("  in v_hxc ",4i6,"  in mesh_set ",4i6)') n1,n2,n3, id,    &
               (ncheck(i),i=1,4)
    write(6,*)

    stop

  endif

    call mesh_set(ipr, "tau", adot, tau, taumsh, ncheck,                 &
        ng, kgv, phase, conj, inds, kmax,                                &
        mxdgve, mxdnst, mxdfft)

  if(ncheck(1) /= n1 .or. ncheck(2) /= n2 .or. ncheck(3) /= n3           &
      .or. ncheck(4) /= id) then
    write(6,*)
    write(6,*)  "  STOPPED in v_hartree_xc:  inconsistency in mesh_set:"
    write(6,'("  in v_hxc ",4i6,"  in mesh_set ",4i6)') n1,n2,n3, id,    &
               (ncheck(i),i=1,4)
    write(6,*)

    stop

  endif

    call xc_cell2(author, tblaha, id, n2, n1,n2,n3,                      &
        rhomsh, taumsh, rholapmsh, vxcmsh, adot, exc, rhovxc, strxc )


  deallocate(rholapmsh)
  deallocate(taumsh)

  else

    call xc_cell( author, id, n2, n1, n2, n3,                            &
        rhomsh, vxcmsh, adot, exc, rhovxc, strxc )

  endif
!-----------------------clr----------------------------------------

! fills chd

  do i = 1,ntot
    chd(i) = cmplx(vxcmsh(i),ZERO,REAL64)
  enddo


! fourier transform exchange and correlation

  call cfft_c16(chd, id, n1,n2,n3, 1, wrkfft, mxdwrk)

! compute v hartree


  vhar(1) = C_ZERO
  fac = 2*PI/vcell
  do i=2,ns
    vhar(i) = fac*den(i)/ek(i)
  enddo

! add the exchange correlation, with symmetrization
! if the fft mesh is reduced noise will appear

  allocate(vxcg(mxdgve))

  call mesh_fold(vxcg, chd, id, n1,n2,n3,                                &
      ng, kgv,                                                           &
      mxdgve, mxdfft)

  call star_of_g_fold(vxc, vxcg, .FALSE.,                                &
      ng, phase, conj, ns, inds, mstar,                                  &
      mxdgve, mxdnst)

  deallocate(dentot)
  deallocate(rhomsh)
  deallocate(vxcmsh)
  deallocate(chd)
  deallocate(wrkfft)
  deallocate(vxcg)

  return
end subroutine v_hartree_xc
