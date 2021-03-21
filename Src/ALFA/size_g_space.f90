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

!>  Calculates the size of arrays in g-space
!>
!>  \author       Jose Luis Martins
!>  \version      4.99
!>  \date         2013-2021
!>  \copyright    GNU Public License v2

subroutine size_g_space(emax, adot, ntrans, mtrx,                        &
    mxdgve, mxdnst, mxdcub)

! Written December 18, 2013 from the g_space.f subroutine. jlm
! Modified, documentation, January 2020. JLM
! Modified, same ordering as g_space, delta. 11 March 2021. JLM
! copyright  Jose Luis Martins/INESC-MN

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)           ::  emax                       !<  kinetic energy cutoff of plane wave expansion (hartree).
  real(REAL64), intent(in)           ::  adot(3,3)                  !<  metric in direct space
  integer, intent(in)                ::  ntrans                     !<  number of symmetry operations in the factor group
  integer, intent(in)                ::  mtrx(3,3,48)               !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group

! output

  integer, intent(out)               ::  mxdgve                     !<  array dimension for g-space vectors
  integer, intent(out)               ::  mxdnst                     !<  array dimension for g-space stars
  integer, intent(out)               ::  mxdcub                     !<  array dimension for 3-index g-space

! local allocatable arrays

  real(REAL64), allocatable     ::  gmtmp(:), gm(:)
  integer,allocatable           ::  kgvtmp(:,:), kgv(:,:)
  integer,allocatable           ::  indx(:)

! local variables
  real(REAL64)  ::  gmax                                            !  gmax is maximum length of the generated g-vectors
  real(REAL64)  ::  vcell, bdot(3,3), gmax2, gmod, gmr, diff
  integer       ::  kmax(3), kcub, km(3), ktran(3,48), kk(3)
  real(REAL64)  ::  aminv(3)
  integer       ::  ng, ns, ncount, ncmax, nckeep
  integer       ::  iprot, itest, iavail                            !  prototype vector in star, vector being tested, first vector not yet indexed
  integer       ::  istar, istop
  real(REAL64)  ::  rgm

! counters

  integer       ::  i, j, k, nc

! parameters

  real(REAL64), parameter :: PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter :: DELTA = 0.0001_REAL64


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

! g = k1*b1+k2*b2+k3*b3

  kcub = (2*kmax(1)+1) * (2*kmax(2)+1) * (2*kmax(3)+1)
  allocate(gmtmp(kcub))
  allocate(kgvtmp(3,kcub))

  ng = 0
  km(1) = 0
  km(2) = 0
  km(3) = 0
  do i=-kmax(1),kmax(1)
  do j=-kmax(2),kmax(2)
  do k=-kmax(3),kmax(3)

    gmod = i*(bdot(1,1)*i+bdot(1,2)*j+bdot(1,3)*k) +                &
           j*(bdot(2,1)*i+bdot(2,2)*j+bdot(2,3)*k) +                &
           k*(bdot(3,1)*i+bdot(3,2)*j+bdot(3,3)*k)

    if (gmod <= gmax2) then
      ng = ng + 1

      gmtmp(ng) = gmod
      kgvtmp(1,ng) = i
      kgvtmp(2,ng) = j
      kgvtmp(3,ng) = k
      if(km(1) < abs(i)) km(1) = abs(i)
      if(km(2) < abs(j)) km(2) = abs(j)
      if(km(3) < abs(k)) km(3) = abs(k)
    endif
  enddo
  enddo
  enddo

  if(ng <= 1) then
    write(6,*)
    write(6,'("  stopped in size_g_space  gmax = ",e10.3,           &
              "  no. of g-vec.= ",i6)') gmax,ng

    stop

  endif

  mxdgve = ng
  mxdcub = (2*km(1)+1) * (2*km(2)+1) * (2*km(3)+1)

! sort the g-vectors according to length and collect
! heapsort algorithm from numerical recipes

  allocate(indx(ng))

  call sort(ng,gmtmp,indx)

  allocate(gm(ng))
  allocate(kgv(3,ng))

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
      do i=1,ntrans
        ktran(1,i) = mtrx(1,1,i)*kgv(1,iprot) +                     &
                     mtrx(1,2,i)*kgv(2,iprot) +                     &
                     mtrx(1,3,i)*kgv(3,iprot)
        ktran(2,i) = mtrx(2,1,i)*kgv(1,iprot) +                     &
                     mtrx(2,2,i)*kgv(2,iprot) +                     &
                   mtrx(2,3,i)*kgv(3,iprot)
        ktran(3,i) = mtrx(3,1,i)*kgv(1,iprot) +                     &
                     mtrx(3,2,i)*kgv(2,iprot) +                     &
                     mtrx(3,3,i)*kgv(3,iprot)
      enddo

!     searches for the inverse of star member

      do i=iprot+1,ng

        if ((kgv(1,i)+kgv(1,iprot) == 0) .and.                      &
            (kgv(2,i)+kgv(2,iprot) == 0) .and.                      &
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
      do i=1,ntrans
        if ((ktran(1,i)-kgv(1,itest) == 0) .and.                    &
            (ktran(2,i)-kgv(2,itest) == 0) .and.                    &
            (ktran(3,i)-kgv(3,itest) == 0)) then
          ncount = ncount + 1
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

!       searches for the inverse of star member

        do i=iavail+1,ng

          if ((kgv(1,i)+kgv(1,iavail) == 0) .and.                   &
              (kgv(2,i)+kgv(2,iavail) == 0) .and.                   &
              (kgv(3,i)+kgv(3,iavail) == 0)) then

!           exchanges i with iavail+1

            if(i .ne. iavail+1) then
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

  if(nckeep == ncmax) then

    write(6,*)
    write(6,*) '  STOPPED in size_g_space.  '
    write(6,*) '  Congratulations you are tackling a huge problem'
    write(6,*) '  Change the kind of nc, ncmax, nckeep'

    stop

  endif

  mxdnst = ns

  deallocate(gm)
  deallocate(kgv)

  return

  end subroutine size_g_space
