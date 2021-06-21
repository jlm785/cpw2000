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

!>  Calculates G-vectors for reciprocal lattice points inside
!>  a sphere of radius gmax.
!>
!>  \author       Jose Luis Martins
!>  \version      4.99
!>  \date         1989-2021
!>  \copyright    GNU Public License v2

subroutine g_space(ipr, emax,                                            &
    adot, ntrans, mtrx, tnp,                                             &
    ng, kgv, phase, conj,                                                &
    inds, kmax, indv, ns, mstar, ek, izstar,                             &
    mxdgve, mxdnst, mxdcub)

! written September 16 1989. jlm
! modified February 18 1990. jlm
! modified 14 October 93. jlm
! modified (aminv) 19 January 1996. jlm
! modified (conj,Hartree) 20-22 March 1999 and 12 April. jlm
! modified for f90 and cleaned, 8 January 2014. JLM
! Modified, documentation, January 2020. JLM
! Modifiied ng = sum mstar check, delta. 9 March 2021. JLM
! copyright INESC-MN/Jose Luis Martins

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdnst                          !<  array dimension for g-space stars
  integer, intent(in)                ::  mxdcub                          !<  array dimension for 3-index g-space

  integer, intent(in)                ::  ipr                             !<  print switch
  real(REAL64), intent(in)           ::  emax                            !<  kinetic energy cutoff of plane wave expansion (hartree).
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  integer, intent(in)                ::  ntrans                          !<  number of symmetry operations in the factor group
  integer, intent(in)                ::  mtrx(3,3,48)                    !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group
  real(REAL64), intent(in)           ::  tnp(3,48)                       !<  2*pi* i-th component (in lattice coordinates) of the fractional translation vector associated with the k-th symmetry operation of the factor group

! output

  integer, intent(out)               ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(out)               ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
  complex(REAL64), intent(out)       ::  phase(mxdgve)                   !<  phase factor of G-vector n
  real(REAL64), intent(out)          ::  conj(mxdgve)                    !<  is -1 if one must take the complex conjugate of x*phase

  integer, intent(out)               ::  ns                              !<  number os stars with length less than gmax
  integer, intent(out)               ::  inds(mxdgve)                    !<  star to which g-vector n belongs
  integer, intent(out)               ::  kmax(3)                         !<  max value of |kgv(i,n)|
  integer, intent(out)               ::  indv(mxdcub)                    !<  kgv(i,indv(jadd)) is the g-vector associated with jadd. jadd is defined by the g-vector components and kmax
  real(REAL64), intent(out)          ::  ek(mxdnst)                      !<  kinetic energy (hartree) of g-vectors in star j
  integer, intent(out)               ::  mstar(mxdnst)                   !<  number of g-vectors in the j-th star
  integer, intent(out)               ::  izstar(mxdnst)                  !<  is 0 if the phase=0

! local allocatable arrays

  real(REAL64), allocatable     ::  gm(:)
  real(REAL64), allocatable     ::  gmtmp(:)
  integer,allocatable           ::  indx(:)
  integer,allocatable           ::  kgvtmp(:,:)

! local variables
  real(REAL64)  ::  gmax                                                 !  gmax is maximum length of the generated g-vectors
  real(REAL64)  ::  vcell, bdot(3,3), gmax2, gmod, gmr, diff
  integer       ::  ktran(3,48), kk(3)
  real(REAL64)  ::  aminv(3)
  integer       ::  ncount, ktot, iadd, istar, istop
  integer       ::  iprot, itest, iavail                                 !  prototype vector in star, vector being tested, first vector not yet indexed
  complex(REAL64)  ::  expifi(48)
  real(REAL64)  ::  rgm
  real(REAL64)  ::  fi
  complex(REAL64)  ::  phsum
  integer       ::  ncmax, nckeep

! counters
  integer       ::  i, j, k, nc

! parameters
  real(REAL64), parameter ::  PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter ::  DELTA = 0.0001_REAL64
  real(REAL64), parameter ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64


! calculates vcell,bdot

  call adot_to_bdot(adot,vcell,bdot)

! gmax is the maximum length of the generated g-vectors

  gmax = 2*sqrt(2*emax)

! estimate limits for kgv(i)
! these will be reset at the end of the routine

  aminv(1) = 2*PI / sqrt(adot(1,1))
  aminv(2) = 2*PI / sqrt(adot(2,2))
  aminv(3) = 2*PI / sqrt(adot(3,3))
  kmax(1) = int(gmax/aminv(1)) + 1
  kmax(2) = int(gmax/aminv(2)) + 1
  kmax(3) = int(gmax/aminv(3)) + 1
  gmax2 = gmax*gmax

! G = k1*b1+k2*b2+k3*b3

  allocate(gm(mxdgve))
  allocate(gmtmp(mxdgve))
  allocate(kgvtmp(3,mxdgve))

  ng = 0
  do i=-kmax(1),kmax(1)
  do j=-kmax(2),kmax(2)
  do k=-kmax(3),kmax(3)

    gmod = i*(bdot(1,1)*i+bdot(1,2)*j+bdot(1,3)*k) +                     &
           j*(bdot(2,1)*i+bdot(2,2)*j+bdot(2,3)*k) +                     &
           k*(bdot(3,1)*i+bdot(3,2)*j+bdot(3,3)*k)

    if (gmod <= gmax2) then

      ng = ng + 1

      if(ng > mxdgve) then

        write(6,*)
        write(6,'("    STOPPED in g_space with gmod = ",e10.3,           &
                & "  and  mxdgve = ",i9)') gmod, mxdgve
        write(6,'("    try  mxdgve of the order of ",i9)')               &
                    int(0.017*vcell*gmax*gmax*gmax)

        stop

      endif

      gmtmp(ng) = gmod
      kgvtmp(1,ng) = i
      kgvtmp(2,ng) = j
      kgvtmp(3,ng) = k

    endif
  enddo
  enddo
  enddo

  if(ng <= 1) then

    write(6,*)
    write(6,'("  STOPPED in g_space  gmax = ",e10.3,                      &
            & "  no. of g-vec.= ",i6)') gmax,ng

    stop

  endif


! sort the g-vectors according to length and collect
! heapsort algorithm from numerical recipes

  allocate(indx(ng))

  call sort(ng,gmtmp,indx)

  do i = 1,ng
    gm(i) = gmtmp(indx(i))
    kgv(1,i) = kgvtmp(1,indx(i))
    kgv(2,i) = kgvtmp(2,indx(i))
    kgv(3,i) = kgvtmp(3,indx(i))
  enddo

  deallocate(gmtmp)
  deallocate(kgvtmp)
  deallocate(indx)


! sort g vectors of same length by the
! kgv(1,i), then kgv(2,i), then the
! kgv(3,i) value.   this makes the kgv
! sort machine independent.

  istar = 2

  do i = 2,ng

    if( i == ng) then
      diff = 1.0
    else
      diff = gm(i+1)-gm(istar)
    endif

    if(diff > delta) then
      istop = i
      do j = istar,istop-1
        do k = j+1,istop
          if (kgv(1,k) < kgv(1,j)) then

            rgm = gm(j)
            kk(1) = kgv(1,j)
            kk(2) = kgv(2,j)
            kk(3) = kgv(3,j)
            gm(j) = gm(k)
            kgv(1,j) = kgv(1,k)
            kgv(2,j) = kgv(2,k)
            kgv(3,j) = kgv(3,k)
            gm(k) = rgm
            kgv(1,k) = kk(1)
            kgv(2,k) = kk(2)
            kgv(3,k) = kk(3)

          elseif (kgv(1,k) == kgv(1,j)) then

            if (kgv(2,k) < kgv(2,j)) then

              rgm = gm(j)
              kk(1) = kgv(1,j)
              kk(2) = kgv(2,j)
              kk(3) = kgv(3,j)
              gm(j) = gm(k)
              kgv(1,j) = kgv(1,k)
              kgv(2,j) = kgv(2,k)
              kgv(3,j) = kgv(3,k)
              gm(k) = rgm
              kgv(1,k) = kk(1)
              kgv(2,k) = kk(2)
              kgv(3,k) = kk(3)

            elseif (kgv(2,k) == kgv(2,j)) then

              if (kgv(3,k) < kgv(3,j)) then

                rgm = gm(j)
                kk(1) = kgv(1,j)
                kk(2) = kgv(2,j)
                kk(3) = kgv(3,j)
                gm(j) = gm(k)
                kgv(1,j) = kgv(1,k)
                kgv(2,j) = kgv(2,k)
                kgv(3,j) = kgv(3,k)
                gm(k) = rgm
                kgv(1,k) = kk(1)
                kgv(2,k) = kk(2)
                kgv(3,k) = kk(3)

              endif

            endif

          endif

        enddo

      enddo

      istar = istop + 1

    endif

  enddo

! collects g-vectors into stars by symmetry
! first star is (0 0 0)

  ns = 1
  ek(1) = ZERO
  phase(1) = cmplx(UM,ZERO,REAL64)
  conj(1) = UM
  inds(1) = 1
  mstar(1) = 1

  iprot = 1
  iavail = 2
  itest = 2

! this loop will be run at most ng*ng times, in practice it will be ~ng times

  ncmax = huge(ng)
  if(nint(sqrt(real(ncmax))) > ng) ncmax =  ng*(ng-1)/2

  do nc = 1,ncmax

    if(itest > ng) then
      diff = 1.0
    else
      diff = abs(gm(itest)-gm(iprot))
    endif

    if(diff > DELTA) then

!     new star

      iprot = iavail
      ns = ns + 1

      if(ns > mxdnst) then

        write(6,*)
        write(6,'("   STOPPED in g_space.  gmax**2 = ",e10.3,            &
           & "   mxdnst = ",i8)') gm(iprot),mxdnst
        write(6,'("   try mxdnst of order ",i8)') ng/ntrans

        stop

      endif

      ek(ns) = gm(iprot) / 2
      do i=1,ntrans
        ktran(1,i) = mtrx(1,1,i)*kgv(1,iprot) +                          &
                     mtrx(1,2,i)*kgv(2,iprot) +                          &
                     mtrx(1,3,i)*kgv(3,iprot)
        ktran(2,i) = mtrx(2,1,i)*kgv(1,iprot) +                          &
                     mtrx(2,2,i)*kgv(2,iprot) +                          &
                   mtrx(2,3,i)*kgv(3,iprot)
        ktran(3,i) = mtrx(3,1,i)*kgv(1,iprot) +                          &
                     mtrx(3,2,i)*kgv(2,iprot) +                          &
                     mtrx(3,3,i)*kgv(3,iprot)
      enddo

      do i = 1,ntrans
        fi = tnp(1,i)*ktran(1,i) + tnp(2,i)*ktran(2,i) +                 &
             tnp(3,i)*ktran(3,i)
        expifi(i) = cmplx(cos(fi),sin(fi),REAL64)
      enddo

      phsum = cmplx(ZERO,ZERO,REAL64)
      ncount = 0
      do i=1,ntrans
        if ((ktran(1,i)-kgv(1,iprot) == 0) .and.                         &
            (ktran(2,i)-kgv(2,iprot) == 0) .and.                         &
            (ktran(3,i)-kgv(3,iprot) == 0)) then
          phsum = phsum + expifi(i)
          ncount = ncount + 1
        endif
      enddo
      phase(iprot) = phsum/ncount
      conj(iprot) = UM
      mstar(ns) = 1
      inds(iprot) = ns

!     searches for the inverse of star member

      do i=iprot+1,ng

        if ((kgv(1,i)+kgv(1,iprot) == 0) .and.                           &
            (kgv(2,i)+kgv(2,iprot) == 0) .and.                           &
            (kgv(3,i)+kgv(3,iprot) == 0)) then

!         exchanges i with iprot

          if(i /= iprot+1) then
            kk(1) = kgv(1,iprot+1)
            kk(2) = kgv(2,iprot+1)
            kk(3) = kgv(3,iprot+1)
            kgv(1,iprot+1) = kgv(1,i)
            kgv(2,iprot+1) = kgv(2,i)
            kgv(3,iprot+1) = kgv(3,i)
            kgv(1,i) = kk(1)
            kgv(2,i) = kk(2)
            kgv(3,i) = kk(3)
            gmr = gm(iprot+1)
            gm(iprot+1) = gm(i)
            gm(i) = gmr
          endif

!         verifies if the inverse is related by symmetry

          ncount = 0
          phsum = cmplx(ZERO,ZERO,REAL64)
          do j=1,ntrans
            if ((ktran(1,j)-kgv(1,iprot+1) == 0) .and.                   &
                (ktran(2,j)-kgv(2,iprot+1) == 0) .and.                   &
                (ktran(3,j)-kgv(3,iprot+1) == 0)) then
              phsum = phsum + expifi(j)
              ncount = ncount + 1
            endif
          enddo

          if(ncount /= 0) then
            phase(iprot+1) = phsum/ncount
            conj(iprot+1) = UM
          else
            phase(iprot+1) = phase(iprot)
            conj(iprot+1) = -UM
          endif

          mstar(ns) = mstar(ns) + 1
          inds(iprot+1) = ns

          exit

        endif

      enddo

!     updates iavail and itest

      iavail = iprot + 2
      itest = iprot + 2

    else

!     we are not starting a new star
!     first check if itest is related by symmetry to iprot

      ncount = 0
      phsum = cmplx(ZERO,ZERO,REAL64)
      do i=1,ntrans
        if ((ktran(1,i)-kgv(1,itest) == 0) .and.                         &
            (ktran(2,i)-kgv(2,itest) == 0) .and.                         &
            (ktran(3,i)-kgv(3,itest) == 0)) then
          ncount = ncount + 1
           phsum = phsum + expifi(i)
        endif
      enddo

      if(ncount /= 0) then

!       it is related by symmetry

!       exchanges itest with iavail

        if(itest /= iavail) then

          kk(1) = kgv(1,iavail)
          kk(2) = kgv(2,iavail)
          kk(3) = kgv(3,iavail)
          kgv(1,iavail) = kgv(1,itest)
          kgv(2,iavail) = kgv(2,itest)
          kgv(3,iavail) = kgv(3,itest)
          kgv(1,itest) = kk(1)
          kgv(2,itest) = kk(2)
          kgv(3,itest) = kk(3)
          gmr = gm(iavail)
          gm(iavail) = gm(itest)
          gm(itest) = gmr

        endif

        phase(iavail) = phsum/ncount
        conj(iavail) = UM
        mstar(ns) = mstar(ns) + 1
        inds(iavail) = ns

!       searches for the inverse of star member

        do i=iavail+1,ng

          if ((kgv(1,i)+kgv(1,iavail) == 0) .and.                        &
              (kgv(2,i)+kgv(2,iavail) == 0) .and.                        &
              (kgv(3,i)+kgv(3,iavail) == 0)) then

!           exchanges i with iavail+1

            if(i /= iavail+1) then

              kk(1) = kgv(1,iavail+1)
              kk(2) = kgv(2,iavail+1)
              kk(3) = kgv(3,iavail+1)
              kgv(1,iavail+1) = kgv(1,i)
              kgv(2,iavail+1) = kgv(2,i)
              kgv(3,iavail+1) = kgv(3,i)
              kgv(1,i) = kk(1)
              kgv(2,i) = kk(2)
              kgv(3,i) = kk(3)
              gmr = gm(iavail+1)
              gm(iavail+1) = gm(i)
              gm(i) = gmr

            endif

!           verifies if the inverse is related by symmetry

           ncount = 0
           phsum = cmplx(ZERO,ZERO,REAL64)
           do j=1,ntrans
              if ((ktran(1,j)-kgv(1,iavail+1) == 0) .and.                &
                  (ktran(2,j)-kgv(2,iavail+1) == 0) .and.                &
                  (ktran(3,j)-kgv(3,iavail+1) == 0)) then
                phsum = phsum + expifi(j)
                ncount = ncount + 1
              endif
            enddo

            if(ncount /= 0) then
              phase(iavail+1) = phsum/ncount
              conj(iavail+1) = UM
            else
              phase(iavail+1) = phase(iavail)
              conj(iavail+1) = -UM
            endif

            mstar(ns) = mstar(ns) + 1
            inds(iavail+1) = ns

            exit

          endif

        enddo

!     updates iavail and itest

        if(itest == iavail) then
          itest = itest +2
        else
          itest = itest +1
        endif
        iavail = iavail +2

      else

!       itest it is not related by symmetry
!       increases itest counter for next vector

        itest = itest + 1
      endif

    endif

    nckeep = nc

    if(iavail > ng) exit

  enddo

  deallocate(gm)

  if(nckeep == ncmax) then

    write(6,*)
    write(6,*) '  STOPPED in g_space.  '
    write(6,*) '  Congratulations you are tackling a huge problem'
    write(6,*) '  Change the kind of nc, ncmax, nckeep'

    stop

  endif

! reset kmax and check size of indv

  kmax(1) = 0
  kmax(2) = 0
  kmax(3) = 0
  do i = 1,ng
    if (kgv(1,i) > kmax(1)) kmax(1)=kgv(1,i)
    if (kgv(2,i) > kmax(2)) kmax(2)=kgv(2,i)
    if (kgv(3,i) > kmax(3)) kmax(3)=kgv(3,i)
  enddo
  ktot = (2*kmax(1)+1)*(2*kmax(2)+1)*(2*kmax(3)+1)

  if(ktot > mxdcub) then

    write(6,*)
    write(6,'("   STOPPED in gspace  change mxdcub",                     &
         & " to at least  ",i9)') ktot

    stop

  endif

! compute index array indv

  do i = 1,ktot
    indv(i) = 0
  enddo
  do i = 1,ng
    iadd = kgv(1,i)+kmax(1)
    iadd = iadd*(2*kmax(2)+1) + kgv(2,i)+kmax(2)
    iadd = iadd*(2*kmax(3)+1) + kgv(3,i)+kmax(3)+1
    indv(iadd) = i
  enddo

! checks stars with zero phase factors

  ncount = 1
  izstar(1) = 1
  do i = 2,ns
    if(abs(phase(ncount+1) - cmplx(UM,ZERO,REAL64)) > DELTA) then
      do j=ncount+1,ncount+mstar(i)
        phase(j) = cmplx(ZERO,ZERO,REAL64)
      enddo
      izstar(i) = 0
    else
      izstar(i) = 1
    endif
    ncount = ncount + mstar(i)
  enddo

! paranoid check

  if(ncount /= ng) then

    write(6,*)
    write(6,'("  STOPPED in g_space. Number of G-vectors in stars ",i9,  &
            & "  inconsistent with total number ",i9)') ncount,ng

        stop
  endif

  if (ipr > 0) then

    write(6,*)
    write(6,'(1x,i8," G-vectors are set up in ",i7," stars,",            &
              & "  kmax = ",3i4)') ng,ns,(kmax(i),i=1,3)

    if (ipr > 1) then

!     print g-vectors

      write(6,*)
      write(6,'("  G-vectors ")')
      write(6,*)
      istop = 0
      do i=1,ns
        write(6,*)
        write(6,'("  Star no",i8,"  with",i3," members.  Kinetic",       &
            & " energy = ",f12.4,4x,i3)') i,mstar(i),ek(i),izstar(i)
        write(6,*)
        write(6,*)
        write(6,'("       i      inds     kx   ky   kz        phase")')
        istar = istop+1
        istop = istar+mstar(i)-1
        do j=istar,istop
          write(6,'(i8,2x,i8,2x,3i5,f10.5,3x,f10.5,5x,f5.2)')            &
               j,inds(j),(kgv(k,j),k=1,3),phase(j),conj(j)
        enddo
      enddo

    endif
  endif

  return

end subroutine g_space
