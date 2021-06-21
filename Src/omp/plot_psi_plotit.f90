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

!>  Plots or prepares the  1D, 2D, or 3D plots of the wave-functions
!>  for a given k-point.

subroutine plot_psi_plotit(ioreplay, nc,                                 &
              adot, ntype, natom, nameat, rat, zv,                       &
              ng, kgv, kmax,                                             &
              mtxd, neig, isort, psi, ei,                                &
              mxddim, mxdbnd, mxdtyp, mxdatm, mxdgve)


! Written 20 February 2018 based on cpw_analysis_sub, out_band_onek,
! and rho_v_plot_sub.
! Documentation, merge of with(out) spin-orbit 3 February 2021. JLM

! copyright  Jose Luis Martins/INESC-MN

! version 4.99

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  ioreplay                        !<  tape number for reproducing calculations

  integer, intent(in)                ::  nc                              !<  = 1, no spin-orbit, = 2 with spin-orbit

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  character(len=2), intent(in)       ::  nameat(mxdtyp)                  !<  chemical symbol for the type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i
  real(REAL64), intent(in)           ::  zv(mxdtyp)                      !<  valence of atom with type i

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
  integer, intent(in)                ::  kmax(3)                         !<  max value of |kgv(i,n)|

  integer, intent(in)                ::  mtxd                            !<  dimension of the hamiltonian
  integer, intent(in)                ::  neig                            !<  number of eigenvectors
  integer, intent(in)                ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian
  complex(REAL64), intent(in)        ::  psi(nc*mxddim,nc*mxdbnd)        !<  |psi> 
  real(REAL64), intent(in)           ::  ei(nc*mxdbnd)                   !<  eigenvalue (hartree)       

! local allocatable arrays

  complex(REAL64), allocatable       ::  denu(:)                         !  unsymmetrized density (in G)
  real(REAL64), allocatable          ::  rhomsh(:)
  real(REAL64), allocatable          ::  wrkfft(:)
  complex(REAL64), allocatable       ::  chd(:)
  complex(REAL64), allocatable       ::  chd1(:)
  complex(REAL64), allocatable       ::  chd2(:)
  integer, allocatable               ::  ipoint(:)

! local variables

  integer         ::  iorb                      !  chosen eigenvector
  integer         ::  iplot                     !  plot dimensionality
  integer         ::  ipsi                      !  1: |psi|**2, 2: Re(psi), 3: Im(psi)

  integer         ::  mxdfft,mxdwrk
  integer         ::  nsfft(3)
  integer         ::  n1, n2, n3, id
  integer         ::  nn1, nn2, nn3
  integer         ::  kd1, kd2, kd3
  integer         ::  k1, k2, k3
  integer         ::  ntot, it
 
! parameters

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter  ::  C_I = cmplx(ZERO,UM,REAL64)
  real(REAL64), parameter     ::  EV = 27.211385_REAL64

! counters

  integer    ::  i, j, k, m



! initial calculations


  call size_fft(kmax, nsfft, mxdfft, mxdwrk)

  if(nc == 1) then
    allocate(chd(mxdfft))
  elseif(nc == 2) then
    allocate(chd(mxdfft))
    allocate(chd1(mxdfft))
    allocate(chd2(mxdfft))
  else
    write(6,'("   STOPPED in plot_psi_plotit:   nc = ",i8)') nc

    stop

  endif

  allocate(wrkfft(mxdwrk))

  n1 = nsfft(1)
  n2 = nsfft(2)
  n3 = nsfft(3)
  id = n1
  ntot = id * n2 * n3
  nn1 = (n1-1) / 2
  nn2 = (n2-1) / 2
  nn3 = (n3-1) / 2
  kd1 = 0
  kd2 = 0
  kd3 = 0

! fills array ipoint (gather-scatter index)

  allocate(ipoint(mtxd))

  do i = 1,mtxd
    it = isort(i)
    k1 = kgv(1,it)
    if (iabs(k1) > kd1) kd1 = iabs(k1)
    if (k1 < 0) k1 = n1 + k1
    k2 = kgv(2,it)
    if (iabs(k2) > kd2) kd2 = iabs(k2)
    if (k2 < 0) k2 = n2 + k2
    k3 = kgv(3,it)
    if (iabs(k3) > kd3) kd3 = iabs(k3)
    if (k3 < 0) k3 = n3 + k3
    ipoint(i) = (k3*n2 + k2)*id + k1 + 1
  enddo
  if (kd1 > nn1 .or. kd2 > nn2 .or. kd3 > nn3) then
    write(6,*)
    write(6,'("     STOPPED in psi_plot_plotit      ",                   &
       &  "size of matrix ",i7," fft mesh ",3i5)') mtxd,n1,n2,n3

    stop

  endif

  allocate(rhomsh(mxdfft))
  allocate(denu(ng))



! loop over possible plots

  do k = 1,100

    write(6,*)
    write(6,*) '  Which orbital do you want to plot?'
    write(6,*) '  1 to ',neig
    write(6,*)
    write(6,*) '  0 ends ploting for this k-point'
    write(6,*)

    read(5,*) iorb
    write(ioreplay,'(2x,i8,"   orbital choice")') iorb

    if(iorb < 0 .or. iorb > neig) then

      iorb = 0
      write(6,*)
      write(6,*) '  Wrong answer '
      write(6,*)

    endif

    if(iorb == 0) then

      write(6,*)
      write(6,*)  '  ending plots for this k-point'
      write(6,*)

      exit

    endif

    write(6,'("  The energy of orbital ",i5," is:",f12.3,"eV")')         &
                           iorb,ei(iorb)*EV
 
!$omp parallel do default(shared) private(i)
    do i = 1,ntot
      rhomsh(i) = ZERO
    enddo
!$omp end parallel do

!$omp parallel do default(shared) private(i)
    do i = 1,ntot
      chd(i) = C_ZERO
    enddo
!$omp end parallel do

    if(nc == 1) then
!$omp parallel do default(shared) private(i)
      do i = 1,mtxd
        chd(ipoint(i)) = psi(i,iorb)
      enddo
!$omp end parallel do

!     fourier transform to real space

      call cfft_wf_c16(chd, id, n1,n2,n3, kd1,kd2,kd3,-1, wrkfft, mxdwrk)

    else

!$omp parallel do default(shared) private(i)
      do i = 1,ntot
        chd1(i) = C_ZERO
      enddo
!$omp end parallel do

!$omp parallel do default(shared) private(i)
      do i = 1,ntot
        chd2(i) = C_ZERO
      enddo
!$omp end parallel do

!$omp parallel do default(shared) private(i)
      do i = 1,mtxd
        chd1(ipoint(i)) = psi(2*i-1,iorb)
        chd2(ipoint(i)) = psi(2*i  ,iorb)
      enddo
!$omp end parallel do

!     fourier transform to real space

      call cfft_wf_c16(chd1, id, n1,n2,n3, kd1,kd2,kd3, -1, wrkfft, mxdwrk)
    
      call cfft_wf_c16(chd2, id, n1,n2,n3, kd1,kd2,kd3, -1, wrkfft, mxdwrk)

    endif

!   loop over components of psi and type of plot

    do j = 1,100

      if(nc == 1) then

        write(6,*)
        write(6,*) '  What do you want to plot?'
        write(6,*)
        write(6,*) '  0)  choose another orbital'
        write(6,*) '  1)  square of absolute psi (default) '
        write(6,*) '  2)  real part of periodic part of psi'
        write(6,*) '  3)  imaginary part of periodic part of psi'
        write(6,*)
        
        read(5,*) ipsi
        write(ioreplay,'(2x,i8,"   choice of component of psi")') ipsi
        
        if(ipsi < 0 .or. ipsi > 3) then

          ipsi = 1
          write(6,*)
          write(6,*) '  Wrong answer, using default'
          write(6,*)

        endif

        if(ipsi == 0) then

          write(6,*)
          write(6,*)  '  ending plots for this orbital'
          write(6,*)

          exit

        elseif(ipsi == 3) then

!$omp parallel do default(shared) private(i)
          do i=1,ntot
            rhomsh(i) = -real(chd(i)*C_I,REAL64)
          enddo
!$omp end parallel do

        elseif(ipsi == 2) then

!$omp parallel do default(shared) private(i)
          do i=1,ntot
            rhomsh(i) = real(chd(i),REAL64)
          enddo
!$omp end parallel do

        else

!$omp parallel do default(shared) private(i)
          do i=1,ntot
            rhomsh(i) = real(chd(i)*conjg(chd(i)),REAL64)
          enddo
!$omp end parallel do

        endif

      else

        write(6,*)
        write(6,*) '  What do you want to plot?'
        write(6,*)
        write(6,*) '  0)  choose another orbital'
        write(6,*) '  1)  square of absolute psi (default)'
        write(6,*) '  2)  real part of periodic part of psi up'
        write(6,*) '  3)  imaginary part of periodic part of psi up'
        write(6,*) '  4)  real part of periodic part of psi down'
        write(6,*) '  5)  imaginary part of periodic part of psi down'
        write(6,*)


        read(5,*) ipsi
        write(ioreplay,'(2x,i8,"   choice of component of psi")') ipsi

        if(ipsi < 0 .or. ipsi > 5) then

          ipsi = 1
          write(6,*)
          write(6,*) '  Wrong answer, using default'
          write(6,*)

        endif

        if(ipsi == 0) then

          write(6,*)
          write(6,*)  '  ending plots for this orbital'
          write(6,*)

          exit

        elseif(ipsi == 5) then

!$omp parallel do default(shared) private(i)
          do i=1,ntot
            rhomsh(i) = -real(chd2(i)*C_I,REAL64)
          enddo
!$omp end parallel do

        elseif(ipsi == 4) then

!$omp parallel do default(shared) private(i)
          do i=1,ntot
            rhomsh(i) = real(chd2(i),REAL64)
          enddo
!$omp end parallel do

        elseif(ipsi == 3) then

!$omp parallel do default(shared) private(i)
          do i=1,ntot
            rhomsh(i) = -real(chd1(i)*C_I,REAL64)
          enddo
!$omp end parallel do

        elseif(ipsi == 2) then

!$omp parallel do default(shared) private(i)
          do i=1,ntot
            rhomsh(i) = real(chd1(i),REAL64)
          enddo
!$omp end parallel do

      else

!$omp parallel do default(shared) private(i)
          do i=1,ntot
            rhomsh(i) = real(chd1(i)*conjg(chd1(i)),REAL64) +              &
                        real(chd2(i)*conjg(chd2(i)),REAL64)
          enddo
!$omp end parallel do

        endif

      endif

!     fourier transform to momentum space


!$omp parallel do default(shared) private(i)
      do i=1,ntot
        chd(i) = cmplx(rhomsh(i),ZERO,REAL64)
      enddo
!$omp end parallel do
           

      call rfft_c16(chd, id, n1,n2,n3, 1, wrkfft, mxdwrk)

!     initialize denu

      call mesh_fold(denu,chd,id,n1,n2,n3,                               &
      ng,kgv,                                                            &
      mxdgve,mxdfft)


      do m = 1,100

        write(6,*)
        write(6,*) '  What kind of plot do you want?'
        write(6,*)
        write(6,*) '  0)  Exit to next wave-function component'
        write(6,*) '  1)  One-dimensional planar and z-average '
        write(6,*) '  2)  Two dimensional (contour) plot'
        write(6,*) '  3)  Three dimensional plot'
        write(6,*)
        
        read(5,*) iplot
        write(ioreplay,'(2x,i8,"   plot dimensionality")') iplot
        
        if(iplot < 0 .or. iplot > 3) then

          iplot = 0
          write(6,*)
          write(6,*) '  Wrong answer '
          write(6,*)

        endif

        if(iplot == 0) then

          write(6,*)
          write(6,*)  '  ending plots for this component'
          write(6,*)

          exit

        elseif(iplot == 1) then

          call plot_average_simple(ioreplay, denu,                       &
              adot, ntype, natom, nameat, rat, zv,                       &
              ng, kgv,                                                   &
              mxdtyp, mxdatm, mxdgve)

        elseif(iplot == 2) then

          call plot_contour(ioreplay, denu, adot, ng, kgv, mxdgve)

        elseif(iplot == 3) then

          call plot_contour3D(ioreplay, denu,                            &
              adot, ntype, natom, nameat, rat,                           &
              ng, kgv,                                                   &
              mxdtyp, mxdatm, mxdgve)

        endif

      enddo                                !  loop over possible dimensionality

    enddo                                  !  loop over possible components

  enddo                                    !  loop over possible orbitals 

  deallocate(chd)
  if(nc == 2) then
    deallocate(chd1)
    deallocate(chd2)
  endif
  deallocate(wrkfft)
  deallocate(ipoint)
  deallocate(rhomsh)
  deallocate(denu)

  return
end subroutine plot_psi_plotit

