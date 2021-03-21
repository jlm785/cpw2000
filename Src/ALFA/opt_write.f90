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

!>  writes to files the several optical functions  
!>  derived from the complex dielectric function


subroutine opt_write(ix, iy, iotape, title, subtitle,                    &
                     nhist, ehist, e_re, e_im, vcell, ztot, ispin)

! Written by Jose Luis Martins extracting previous code
! by Carlos Loia Reis. 20 October 2020
! copyright  Carlos Loia Reis/INESC-MN

! version 4.98

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)              ::  ix, iy                            !<  tensor components

  integer, intent(in)              ::  iotape                            !<  tape number
  character(len=50), intent(in)    ::  title                             !<  title for plots
  character(len=140), intent(in)   ::  subtitle                          !<  subtitle for plots

  integer, intent(in)              ::  nhist                             !<  number of histogram points
  real(REAL64), intent(in)         ::  ehist(nhist)                      !<  energies    
  real(REAL64), intent(in)         ::  e_re(nhist), e_im(nhist)          !<  real and imaginary parts of the dielectric function

  real(REAL64), intent(in)         ::  vcell                             !<  volume of primitive cell
  real(REAL64), intent(in)         ::  ztot                              !<  total number of electrons
  integer, intent(in)              ::  ispin                             !<  ispin = 1 with spin-orbit, ispin = 2 no spin-orbit
   
! local allocatable arrays

  real(REAL64), allocatable        ::  fhist(:)                         !  quantity to be plotted  
  real(REAL64), allocatable        ::  n_re(:), n_im(:)                 !  real and imaginary part of refractive index

! local variables

  real(REAL64)       ::  e_mod
  complex(REAL64)    ::  ratio

! filenames

  character(len=6)   ::  fsuffix(0:3,0:3)
  character(len=4)   ::  flso
  character(len=8)   ::  fprefix

  character(len=40)  ::  filename_out
  character(len=50)  ::  func                                            !  name of function

! constants

  real(REAL64), parameter    ::  UM = 1.0_REAL64, ZERO = 0.0_REAL64
  complex(REAL64), parameter ::  C_I = cmplx(ZERO,UM,REAL64)

  real(REAL64), parameter    ::  TENM6 = 1.0E-6_REAL64
  real(REAL64), parameter    ::  TENM2 = 0.01_REAL64

  real(REAL64), parameter    ::  HARTREE = 27.211386246_REAL64
  real(REAL64), parameter    ::  EPSILON_0 = 8.854187813E-12_REAL64
  real(REAL64), parameter    ::  AUT = 2.4188843266E-17_REAL64
  real(REAL64), parameter    ::  ALPHA = UM / 137.036
  real(REAL64), parameter    ::  BOHR = 0.5291772109E-10_REAL64 

! counters

  integer          :: i

  fprefix = 'optical_'

  fsuffix(0,0) = ".dat"

  fsuffix(1,1) = "xx.dat"
  fsuffix(2,2) = "yy.dat"
  fsuffix(3,3) = "zz.dat"

  fsuffix(1,2) = "xy.dat"
  fsuffix(2,3) = "yz.dat"
  fsuffix(1,3) = "xz.dat"

  fsuffix(2,1) = "yx.dat"
  fsuffix(3,2) = "zy.dat"
  fsuffix(3,1) = "zx.dat"

  if (ispin == 1) then
     if(ix == 0 .and. iy == 0) then
        flso = "_so"
     else
        flso = "_so_"
     endif
  else
     if(ix == 0 .and. iy == 0) then
        flso = " "
     else
        flso = "_"
     endif
  endif

  allocate(fhist(nhist))
  allocate(n_re(nhist))
  allocate(n_im(nhist))

  do i = 1,nhist
    e_mod = sqrt(e_re(i)*e_re(i) + e_im(i)*e_im(i))
    n_re(i) = sqrt( (e_mod + e_re(i)) / 2 )
    n_im(i) = sqrt( (e_mod - e_re(i)) / 2 )
  enddo



  filename_out = fprefix//"epsilon_im"//trim(flso)//fsuffix(ix,iy)

  func = 'Imaginary Part of Dielectric Function             '

  call opt_write_file(trim(filename_out), iotape, title, subtitle, func, &
                      nhist, ehist, e_im, vcell, ztot)


  filename_out = fprefix//"epsilon_re"//trim(flso)//fsuffix(ix,iy)

  func = 'Real Part of Dielectric Function                  '

  call opt_write_file(trim(filename_out), iotape, title, subtitle, func, &
                      nhist, ehist, e_re, vcell, ztot)



  do i = 1,nhist
    fhist(i) = TENM6 * ehist(i)*e_im(i) * EPSILON_0/AUT 
  enddo

  filename_out = fprefix//"sigma_re"//trim(flso)//fsuffix(ix,iy)

  func = 'Real Optical conductivity (MSiemen/meter)         '

  call opt_write_file(trim(filename_out), iotape, title, subtitle, func, &
                      nhist, ehist, fhist, vcell, ztot)


  do i = 1,nhist
    fhist(i) = -TENM6 * ehist(i)*(e_re(i)-1) * EPSILON_0/AUT
  enddo

  filename_out = fprefix//"sigma_im"//trim(flso)//fsuffix(ix,iy)

  func = 'Imaginary Optical conductivity (MSiemen/meter)    '

  call opt_write_file(trim(filename_out), iotape, title, subtitle, func, &
                      nhist, ehist, fhist, vcell, ztot)



  filename_out = fprefix//"refractive_n_re"//trim(flso)//fsuffix(ix,iy)

  func = 'Refractive index n (real part of n)               '

  call opt_write_file(trim(filename_out), iotape, title, subtitle, func, &
                      nhist, ehist, n_re, vcell, ztot)


  filename_out = fprefix//"refractive_k_im"//trim(flso)//fsuffix(ix,iy)

  func = 'Refractive index k (imaginary part of n)          '

  call opt_write_file(trim(filename_out), iotape, title, subtitle, func, &
                      nhist, ehist, n_im, vcell, ztot)



  do i = 1,nhist
    fhist(i) = TENM2 * 2*ehist(i)*n_im(i) * ALPHA/BOHR
  enddo

  filename_out = fprefix//"absorption_coef"//trim(flso)//fsuffix(ix,iy)

  func = 'Absorption coeficient in cm^-1                    '

  call opt_write_file(trim(filename_out), iotape, title, subtitle, func, &
                      nhist, ehist, fhist, vcell, ztot)



  do i = 1,nhist
    ratio = (n_re(i) + C_I*n_im(i) - 1) / (n_re(i) + C_I*n_im(i) + 1)
    fhist(i) = real(ratio*conjg(ratio),REAL64)
  enddo

  filename_out = fprefix//"reflectivity"//trim(flso)//fsuffix(ix,iy)

  func = 'Reflectivity                                      '

  call opt_write_file(trim(filename_out), iotape, title, subtitle, func, &
                      nhist, ehist, fhist, vcell, ztot)

 

  do i = 1,nhist
    fhist(i) = real(C_I/(e_re(i) + C_I*e_im(i)),REAL64)
  enddo

  filename_out = fprefix//"loss_function"//trim(flso)//fsuffix(ix,iy)

  func = 'Loss function                                     '

  call opt_write_file(trim(filename_out), iotape, title, subtitle, func, &
                      nhist, ehist, fhist, vcell, ztot)

  deallocate(fhist)
  deallocate(n_re)
  deallocate(n_im)

  return
end subroutine opt_write


