!>    Fits an equation of state (Birch or Murnaghan)
!>    to a few E(V) data points for several structures.
!>    Different ways of ploting the results are available.
!>    The Murnaghan fit was written by D. Wood
!>    The Birch fit was added by jlm. Sep 89
!>    Changed to hartree jlm. jun 04
!>    Changed to f90 and gnuplot. June 2022. JLM

!>    Murnaghan equation of state:
!>    F. D. Murnaghan, Proc. Natl. Acad. Sci. 30, 244 (1944).
!>    Given a set of ordered pairs (x(i),y(i)),i = 1,npt
!>    where npt is the total number of points provided to the fit.
!>    Fit is to the form
!>    y = total energy per primitive cell in Hartree
!>      =  a + b*x + c*x**d
!>    where x  =  volume of primitive cell in cubic atomic units

!>    Birch equation of state:
!>    F. Birch J. Geophys. Res. 83, 1257 (1978).
!>    The funcyion y is different...

      program eqst

      implicit none

      integer, parameter  ::  REAL64 = selected_real_kind(12)
      integer, parameter  ::  mxdstr = 10                                !<  maximum number of structures
      integer, parameter  ::  mxdnpt = 40                                !<  maximum number of energy points (per structure)

      character(len=5)          ::  ftype                                !<  type of function
      integer                   ::  nstr                                 !<  number of structures
      real(REAL64)              ::  volfac                               !<  volume is volfac*alatt**3
      integer                   ::  npt                                  !<  number of calculated energies

      real(REAL64)              ::  xar, yar, xscl, yscl, fn
      common/params/ xar(mxdnpt),yar(mxdnpt),xscl(mxdnpt),                &
     &          yscl(mxdnpt),fn,npt

      real(REAL64)              ::  alat(mxdnpt)
      real(REAL64)              ::  wvzero(mxdstr),wezero(mxdstr)
      real(REAL64)              ::  wbzero(mxdstr),wbprim(mxdstr)
      real(REAL64)              ::  xarz(mxdnpt,mxdstr)
      real(REAL64)              ::  yarz(mxdnpt,mxdstr)
      integer                   ::  nptz(mxdstr)
      character(len=10)         ::  label(mxdstr)

      real(REAL64)              ::  a, b, c, dd
      common/rmurna/ a,b,c,dd

      real(REAL64), external    ::  fmurna,fbirch, apmurn
      real(REAL64), external    ::  zbrent

      real(REAL64)              ::  y, yav, y2av
      real(REAL64)              ::  vary, sqvry
      real(REAL64)              ::  bpmin, bpmax                         !<  bounds for B'
      real(REAL64)              ::  below, above, eps
      real(REAL64)              ::  emin
      integer                   ::  nmin
      real(REAL64)              ::  sq, qq
      real(REAL64)              ::  dy, percnt

      real(REAL64)              ::  goodns, qual                         !<  quality of fit
      real(REAL64)              ::  aeq                                  !<  equilibrium lattice constant
      real(REAL64)              ::  bprim                                !<  B'
      real(REAL64)              ::  bzero                                !<  B(0)
      real(REAL64)              ::  bulkmd                               !<  Bulk modulus
      real(REAL64)              ::  ezero                                !<  E(aeq)
      real(REAL64)              ::  vmin
      real(REAL64)              ::  vzero

      integer                   ::  ioerr
      character(len=100)        ::  fline

      real(REAL64)              ::  pmax
      integer                   ::  npts



!     constants

      real(REAL64), parameter   ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64

!     counters

      integer                   ::  ns, mm, ii, jj, n



!     reads the number of structures and type of fit

      read(5,*) ftype
      if(ftype /= 'MURNA' .and. ftype /= 'murna' .and.                   &
     &   ftype /= 'BIRCH' .and. ftype /= 'birch') stop 'ftype'

      read(5,*) nstr
      if(nstr > mxdstr) stop 'nstr'

!     loop over structures

      do ns = 1,nstr

!       read number of points

        read(5,*) npt
        if(npt > mxdnpt .or. npt < 4) stop 'npt'
        fn = UM*npt

!       volfac is the volume of primitive cell in (lattice constant)**3,
!       e.g. for fcc bravais lattice, volfac = .25
!       if volfac is negative volume instead of lattice constants are read

        do ii = 1,100
          fline(ii:ii) = ' '
        enddo

        read(5,'(a100)') fline
        read(fline,*,iostat=ioerr) volfac, label(ns)

        if(ioerr /= 0) then
          read(fline,*) volfac
          label(ns) = '          '
        endif

!       reads in the data points. the lattice constant is
!       in atomic units, the energy in Hartree and converted to Rydberg.
!       For angstroms exchange commented lines

        do mm = 1,npt
          read (5,*) alat(mm),yar(mm)
          if(volfac < ZERO) then
            xar(mm) = alat(mm)
!           xar(mm) = alat(mm)/0.529177**3
          else
            xar(mm) = volfac*(alat(mm))**3
!           xar(mm) = volfac*(alat(mm)/0.529177)**3
          endif
!         back to Rydberg :-(
          yar(mm)  =  2*yar(mm)
        enddo

!       evaluate ybar,y2bar for use in the variance
!       and average of each variable, in order to scale the array y
!       to make the least squares algorithm more stable

        yav = ZERO
        y2av =ZERO
        do ii = 1,npt
          y = yar(ii)
          yav = yav+y
          y2av = y2av+y*y
        enddo

        vary = (y2av-yav*yav/fn)/fn
        yav = yav/fn
        sqvry = sqrt(vary)
        write(6,901) ftype,npt,yav,sqvry
  901  format(/,2x,a5,' fit from ',i3,' (x,y) data pairs read in',//,    &
     &   5x,'mean of y = ',f10.5,2x,'variance = ',f10.5,5x,              &
     &   '---y array are re-scaled',//)

        do ii = 1,npt
          xscl(ii) = xar(ii)
          yscl(ii) = (yar(ii)-yav)/sqvry
        enddo

!       fit of the equation of state

        if(ftype == 'murna' .or. ftype == 'MURNA') then

!         dd is evaluated from a non linear fit
!         limits for the search of optimum bprim
!         fmurna(below)*fmurna(above) < 0

          bpmin = 2*UM
          bpmax = 7*UM
          below  =  UM-bpmax
          above  =  UM-bpmin
          eps = 0.0001

          dd = zbrent(fmurna,below,above,eps)

        else

!         estimates the equilibrium volume

          eps = 0.000001
          nmin = 1
          emin = yscl(1)
          do n = 1,npt
            if(yscl(n) < emin) then
              nmin = n
              emin = yscl(n)
            endif
          enddo
          vmin = xar(nmin)
          below = 0.75*vmin
          above = 1.5*vmin
          do n = 1,npt
            if(xar(n) .gt. below .and. xar(n) .lt. vmin-eps)             &
     &          below = xar(n)
            if(xar(n) .lt. above .and. xar(n) .gt. vmin+eps)             &
     &          above = xar(n)
          enddo

          call brent(fbirch,below,vmin,above,ezero,vzero,bzero,          &
     &               bprim,eps)
        endif

        write(6,902)
 902  format(2x,'cols are x(in),y(in),alat,y(in,scaled),',               &
     &      'y(fit,scaled), and % diff',/)

!       use sq for sum of (y-yscl)**2

        sq = ZERO
        do n = 1,npt
          if(ftype == 'murna' .or. ftype == 'MURNA') then
            qq = apmurn(xscl(n))
          else
            call eofv(qq,xscl(n),ezero,vzero,bzero,bprim)
          endif
          dy = yscl(n)-qq
          percnt = 100.*dy
          sq = sq+dy*dy
          write(6,903) xar(n),yar(n),alat(n),yscl(n),qq,percnt
 903  format(2x,f10.5,4x,2(f10.5,2x),2x,2(f10.5,2x),4x,f10.5,3x,f8.3)
        enddo
        if(npt > 4) then
          goodns = sqrt(sq/float(npt-4))
          write(6,904) goodns
 904  format(/,2x,'goodness of fit = scaled variance = ',                &
     &     '(mn sqr)/(npt-4) = ',/,3x,f12.8,/)
          qual = sqvry*goodns/sqrt(float(npt))
          write(6,905) qual
  905  format(/3x,'typical error in y',3x,e14.6,/)
       else
          write(6,906)
 906  format(/,2x,'number of data pairs read in  =  4',                  &
     &      '--perfect fit',/)
        endif

        if(ftype == 'murna' .or. ftype == 'MURNA') then

!         now find the coefficients a,b,c,d for the original, unscaled y variabl

          a = yav+a*sqvry
          b = sqvry*b
          c = sqvry*c
          write(6,907) a,b,c,dd
 907  format(5x,'for unscaled quantities, y = a+bx+c*x**d where ',/,     &
     &     4x,'a,b,c,d= ',1pe18.10,3e18.10/)
          bprim = UM-dd
          bzero = b*bprim
          vzero = (-c*dd/b)**(UM/bprim)
          ezero = a-b*bprim*vzero/dd
        else
          ezero = yav+ezero*sqvry
          bzero = bzero*sqvry
        endif

!       now calculate the actual bulk modulus in gpa

        bulkmd = 14710.8*bzero
        if(volfac > ZERO) then
          aeq = (vzero/volfac)**(1./3.)
        else
          aeq = vzero**(1./3.)
        endif
        write(6,908) aeq,vzero,ezero,bulkmd,bprim
 908  format(/,'   aeq(a.u.)=',f8.4,'   Veq(a.u.)=',f9.3,            &
     &    '   E0(Ryd)=',f10.5,/,'   B0(GPa)=',f6.2,'   B0PRIM= ',f6.2,/)

!       stores for plot subroutine

        wvzero(ns)  =  vzero
        wezero(ns)  =  ezero
        wbzero(ns)  =  bzero
        wbprim(ns)  =  bprim
        nptz(ns)  =  npt
        do jj = 1,npt
          xarz(jj,ns) = xar(jj)
          yarz(jj,ns) = yar(jj)
        enddo
      enddo


       read(5,*) pmax,npts

      call gpeqst(wezero,wvzero,wbzero,wbprim,nstr,label,ftype,pmax,npts)
!      call tdeqst(wezero,wvzero,wbzero,wbprim,nstr,ftype)
!      call tdeofv(xarz,yarz,nptz,wezero,wvzero,wbzero,                   &
!     &            wbprim,nstr,ftype)

      call gpeofv(xarz,yarz,nptz,wezero,wvzero,wbzero,                   &
     &            wbprim,nstr,label,ftype)

      stop
      end program eqst



      function fmurna(d)

!     Murnaghan equation of state

      implicit none
      integer, parameter  :: REAL64 = selected_real_kind(12)
      integer, parameter  ::  mxdnpt = 40                                !<  maximum number of energy points (per structure)

!     input

      real(REAL64), intent(in)    ::  d

!     value

      real(REAL64)                ::  fmurna

!     commons

      integer                   ::  npt                                  !<  number of calculated energies
      real(REAL64)              ::  xar, yar, xscl, yscl, fn
      common/params/ xar(mxdnpt),yar(mxdnpt),xscl(mxdnpt),               &
     &                  yscl(mxdnpt),fn,npt

      real(REAL64)              ::  a, b, c, dd
      common/rmurna/ a,b,c,dd

!     local

      real(REAL64)              ::  sx, sy,sxd, sxy, sx2, sx1pd
      real(REAL64)              ::  sxdy, sx2d, sxdp1l
      real(REAL64)              ::  sxdlxy, sxdlx, sx2dlx
      real(REAL64)              ::  x, alx, xd, y
      real(REAL64)              ::  alpha, beta, gamma
      real(REAL64)              ::  u, w

!     parameters

      real(REAL64), parameter   ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64

!     counters

      integer    ::  i


      dd = d
      sx = ZERO
      sy = ZERO
      sxd = ZERO
      sxy = ZERO
      sx2 = ZERO
      sx1pd = ZERO
      sxdy = ZERO
      sx2d = ZERO
      sxdlxy = ZERO
      sxdlx = ZERO
      sx2dlx = ZERO
      sxdp1l = ZERO
      do i = 1,npt
        x = xscl(i)
        alx = log(x)
        xd = x**d
        y = yscl(i)
        sx = sx+x
        sy = sy+y
        sxd = sxd+xd
        sxy = sxy+x*y
        sx2 = sx2+x*x
        sx1pd = sx1pd+x*xd
        sxdy = sxdy+y*xd
        sx2d = sx2d+xd*xd
        sxdlxy = sxdlxy+alx*y*xd
        sxdlx = sxdlx+alx*xd
        sx2dlx = sx2dlx+xd*xd*alx
        sxdp1l = sxdp1l+x*alx*xd
      enddo
      alpha = -sxd/sx2d
      beta = -sx1pd/sx2d
      gamma = sxdy/sx2d
      u = -(sx+alpha*sx1pd)/(sx2+beta*sx1pd)
      w = (sxy-gamma*sx1pd)/(sx2+beta*sx1pd)
      a = -(sxd*(gamma+beta*w)+sx*w-sy)/(fn+u*sx+(alpha+beta*u)*sxd)
      b = u*a+w
      c = (alpha+beta*u)*a+beta*w+gamma
! these are the fit parameters
      fmurna = sxdlxy-a*sxdlx-b*sxdp1l-c*sx2dlx
!      write(6,777) a,b,c,dd,f
! 777 format(2x,'a,b,c,dd,f =  ',5e14.5)

      return
      end function fmurna


      function apmurn(x)

!     Murnaghan equation of state

      implicit none
      integer, parameter  :: REAL64 = selected_real_kind(12)

!     input

      real(REAL64), intent(in)    ::  x

!     value

      real(REAL64)                ::  apmurn
      real(REAL64)                ::  a, b, c, dd
      common/rmurna/ a,b,c,dd

      apmurn = a + b*x + c*x**dd

      return
      end function apmurn


      function zbrent(func,x1,x2,tol)

      implicit none
      integer, parameter  :: REAL64 = selected_real_kind(12)

!     input

      real(REAL64), intent(in)    ::  x1, x2
      real(REAL64), intent(in)    ::  tol

      real(REAL64), external      ::  func

!     value

      real(REAL64)                ::  zbrent

!     local

      real(REAL64)      ::  a, b, c, d, e
      real(REAL64)      ::  fa, fb, fc
      real(REAL64)      ::  p, q , r , s
      real(REAL64)      ::  xm
      real(REAL64)      ::  tol1, tol2

      logical           ::  lfail

!     parameters

      integer, parameter        ::  ITMAX = 100
      real(REAL64), parameter   ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
      real(REAL64), parameter   ::  EPS = 3.0e-8_REAL64

!     counters

      integer    ::  iter

      a = x1
      b = x2
      fa = func(a)
      fb = func(b)

      if(fb*fa > 0.) stop  'root must be bracketed for zbrent'

      fc = fb
      lfail = .TRUE.
      do iter = 1,ITMAX
        if(fb*fc > ZERO) then
          c = a
          fc = fa
          d = b-a
          e = d
        endif
        if(abs(fc) < abs(fb)) then
          a = b
          b = c
          c = a
          fa = fb
          fb = fc
          fc = fa
        endif
        tol1 = 2*EPS*abs(b)+tol/2
        xm = (c-b)/2
        if(abs(xm) <= tol1 .or. fb == 0.)then
          lfail = .FALSE.

          exit

        endif
        if(abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
          s = fb/fa
          if(a == c) then
            p = 2*xm*s
            q = UM-s
          else
            q = fa/fc
            r = fb/fc
            p = s*(2*xm*q*(q-r)-(b-a)*(r-UM))
            q = (q-UM)*(r-UM)*(s-UM)
          endif
          if(p > ZERO) q = -q
          p = abs(p)
          if(2*p < min(3*xm*q-abs(tol1*q),abs(e*q))) then
            e = d
            d = p/q
          else
            d = xm
            e = d
          endif
        else
          d = xm
          e = d
        endif
        a = b
        fa = fb
        if(abs(d) > tol1) then
          b = b+d
        else
          b = b+sign(tol1,xm)
        endif
        fb = func(b)
      enddo

      if(lfail) stop 'zbrent exceeded maximum iterations'

      zbrent = b

      return
      end function zbrent

      subroutine gpeofv(xarz,yarz,nptz,ez,vz,bz,bp,nstr,label,ftype)

!     writes tape for td plotting software

      implicit none
      integer, parameter  ::  REAL64 = selected_real_kind(12)
      integer, parameter  ::  mxdstr = 10                                !<  maximum number of structures
      integer, parameter  ::  mxdnpt = 40                                !<  maximum number of energy points (per structure)

!     input

      real(REAL64), intent(in)      ::  xarz(mxdnpt,mxdstr)              !<
      real(REAL64), intent(in)      ::  yarz(mxdnpt,mxdstr)

      integer, intent(in)           ::  nptz(mxdstr)
      real(REAL64), intent(in)      ::  ez(mxdstr)
      real(REAL64), intent(in)      ::  vz(mxdstr)
      real(REAL64), intent(in)      ::  bz(mxdstr)
      real(REAL64), intent(in)      ::  bp(mxdstr)
      integer, intent(in)           ::  nstr                             !<  number of structures
      character(len=10)             ::  label(mxdstr)                    !<  label of the structure
      character(len=5), intent(in)  ::  ftype                            !<  type of function

!     local

      real(REAL64)          ::  xminz(mxdstr), xmaxz(mxdstr)
      real(REAL64)          ::  xmin, xmax, emin, emax
      integer               ::  npt
      real(REAL64)          ::  a, b, c, dd
      real(REAL64)          ::  x, y
      real(REAL64)          ::  step

      integer               ::  iodat, iocom

!     constants

      real(REAL64), parameter   ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
      real(REAL64), parameter   ::  EV = 13.606

!     counters

      integer                   ::  i, j

       iodat = 9
       iocom = 10
       open(unit=iodat,file='eofv_dat.gp')
       open(unit=iocom,file='eofv_com.gp')

       xmin = xarz(1,1)
       xmax = xarz(1,1)
       emin = yarz(1,1)
       emax = yarz(1,1)
       do j = 1,nstr
         xminz(j) = xarz(1,j)
         xmaxz(j) = xarz(1,j)
         npt = nptz(j)
         do i = 1,npt
           if(xarz(i,j) < xminz(j)) xminz(j) = xarz(i,j)
           if(xarz(i,j) > xmaxz(j)) xmaxz(j) = xarz(i,j)
           if(yarz(i,j) < emin) emin = yarz(i,j)
           if(yarz(i,j) > emax) emax = yarz(i,j)
         enddo
         if(xminz(j) < xmin) xmin = xminz(j)
         if(xmaxz(j) > xmax) xmax = xmaxz(j)
       enddo
       emin = emin-0.05*(emax-emin)
       emin = emin*EV
       emax = emax*EV

       write(iocom,'("set term wxt persist")')
       write(iocom,'("# set output ''eofv.pdf''")')
       write(iocom,'("# set term pdf color font ''Helvetica,12''")')
       write(iocom,'("set title ''Equation of state E(V)''")')

       write(iocom,'("set format x ''%.0f''")')
       write(iocom,'("set format y ''%.1f''")')

       write(iocom,'("set yrange [",f10.2,":",f10.2,"]")')  emin, emax
       write(iocom,'("set xrange [",f8.1,":",f8.1,"]")')  xmin, xmax

       write(iocom,'("set xlabel ''V (a.u.)''")')
       write(iocom,'("set ylabel ''E (eV)'' ")')
       write(iocom,*)

       write(iocom,'("plot \")')

       do j = 1,nstr
         npt = nptz(j)
         do i = 1,npt
           write(iodat,'(2f15.8)') xarz(i,j),yarz(i,j)*EV
         enddo
         write(iodat,*)
         write(iodat,*)
         write(iocom,'("''eofv_dat.gp'' using 1:2 index ",i3," w p ls ",    &
     &        i3," ti ''",a10," calculated'', \")')   &
     &            j-1, j,label(j)
       enddo

       if(ftype == 'murna' .or. ftype == 'MURNA') then
         do j = 1,nstr
           dd = 1.0-bp(j)
           b = bz(j)/bp(j)
           c = -(vz(j)**bp(j))*b/dd
           a = ez(j)+b*bp(j)*vz(j)/dd
           step = 0.01*(xmaxz(j)-xminz(j))
           do i = 1,121
             x = xminz(j) + dble(i-11)*step
             y = a + b*x + c*x**dd
             write(iodat,'(2f15.8)') x,y*EV
           enddo
           write(iocom,'("''eofv_dat.gp'' using 1:2 index ",i3,          &
     &            " w li ls ",i3," ti ''",a10," interpolated'', \")')    &
     &            nstr+j-1, j,label(j)
           write(iodat,*)
           write(iodat,*)
         enddo
       else
         do j = 1,nstr
           step = 0.01*(xmaxz(j)-xminz(j))
           do i = 1,121
             x = xminz(j) + dble(i-11)*step
             call eofv(y,x,ez(j),vz(j),bz(j),bp(j))
             write(iodat,'(2f15.8)') x,y*EV
           enddo
           write(iocom,'("''eofv_dat.gp'' using 1:2 index ",i3," w li ls ",i3,", \")')   &
     &            nstr+j-1, j
           write(iodat,*)
           write(iodat,*)
         enddo
       endif
       write(iocom,*)

       close(unit=iocom)
       close(unit=iodat)

       return
       end subroutine gpeofv


      subroutine gpeqst(ez,vz,bz,bp,nstr,label,ftype,pmax,npts)

!     writes files for gnuplot plot of the equation of state

      implicit none
      integer, parameter  ::  REAL64 = selected_real_kind(12)
      integer, parameter  ::  mxdstr = 10                                !<  maximum number of structures

!     input

      real(REAL64), intent(in)      ::  ez(mxdstr)
      real(REAL64), intent(in)      ::  vz(mxdstr)
      real(REAL64), intent(in)      ::  bz(mxdstr)
      real(REAL64), intent(in)      ::  bp(mxdstr)
      integer, intent(in)           ::  nstr                             !<  number of structures
      character(len=10)             ::  label(mxdstr)                    !<  label of the structure
      character(len=5), intent(in)  ::  ftype                            !<  type of function

      real(REAL64), intent(in)      ::  pmax                             !<  maximum pressure for plot
      integer, intent(in)           ::  npts                             !<  number of points in plot

!     local

      real(REAL64)                  ::  h(mxdstr)

      real(REAL64)          ::  emin, vref, hmin
      integer               ::  jmin, jmp
      real(REAL64)          ::  p                                        !  pressure
      real(REAL64)          ::  xx, vr, vmin

      integer               ::  iodat, iocom

!     constants

      real(REAL64), parameter   ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
      real(REAL64), parameter   ::  GPA = 14710.8

!     counters

      integer                   ::  i, j

       iodat = 9
       iocom = 10
       open(unit=iodat,file='eqst_dat.gp')
       open(unit=iocom,file='eqst_com.gp')

       write(iocom,'("set term wxt persist")')
       write(iocom,'("# set output ''eqst.pdf''")')
       write(iocom,'("# set term pdf color font ''Helvetica,12''")')
       write(iocom,'("set title ''Equation of state p(V)''")')

       write(iocom,'("set format x ''%.0f''")')
       write(iocom,'("set format y ''%.2f''")')

!      gets volumes.  Uses Birch for estimate

       emin = ez(1)
       jmin = 1
       do i = 1,nstr
         if(ez(i) < emin) then
           jmin = i
           emin = ez(i)
         endif
       enddo
       vref = vz(jmin)
       vmin = UM
       do i = 1,nstr
         call vofp(vr,pmax/GPA,vz(i),bz(i),bp(i))
         if(vr/vref < vmin) vmin = vr/vref
       enddo
       vmin = 0.95*vmin

       write(iocom,'("set yrange [",f8.2,":",f8.2,"]")')  vmin, UM
       write(iocom,'("set xrange [",f8.1,":",f8.1,"]")')  ZERO, pmax

       write(iocom,'("set xlabel ''p (GPa)''")')
       write(iocom,'("set ylabel ''V / V_0'' ")')
       write(iocom,*)

       write(iocom,'("plot \")')

       do j = 1,npts
         p = (j-1)*pmax/(GPA*(npts-1))
         jmp = jmin

         if(ftype == 'murna' .or. ftype == 'MURNA') then
           do i = 1,nstr
             xx = UM + bp(i)*p/bz(i)
             xx = xx ** ((bp(i)-UM)/bp(i))
             h(i) = ez(i) + bz(i)*vz(i)*(xx-UM)/(bp(i)-UM)
           enddo
         else
           do i = 1,nstr
             call hofp(h(i),p,ez(i),vz(i),bz(i),bp(i))
           enddo
         endif

         hmin = h(1)
         jmin = 1
         do i = 1,nstr
           if(h(i) < hmin) then
             jmin = i
             hmin = h(i)
           endif
         enddo

         if(jmin /= jmp) then
           write(6,'(" TRANSITION PRESSURE P= ",F8.3,"  BETWEEN",2I4)')  &
     &            p*GPA,jmp,jmin
! 200   FORMAT(' TRANSITION PRESSURE P= ',F8.3,'  BETWEEN',2I4)
           write(iocom,'("''eqst_dat.gp'' using 1:2 index ",i3,          &
     &            " w li ls ",i3," ti ''",a10,"'', \")')                 &
     &            jmp-1, jmp,label(jmp)
           write(iodat,*)
           write(iodat,*)
        endif

         if(ftype == 'murna' .or. ftype == 'MURNA') then
           xx = UM + bp(jmin)*p/bz(jmin)
           xx = UM / xx ** (UM/bp(jmin))
           vr = xx*vz(jmin)/vref
         else
           call vofp(vr,p,vz(jmin),bz(jmin),bp(jmin))
           vr=vr/vref
         endif

         write(iodat,'(2f15.8)') p*GPA, vr
! 100   FORMAT(2F15.8)

       enddo
       write(iocom,'("''eqst_dat.gp'' using 1:2 index ",i3," w li ls ",  &
     &      i3," ti ''",a10,"''")')                                     &
     &      jmp-1, jmp,label(jmp)

       close(unit=iocom)
       close(unit=iodat)
       return
       end subroutine gpeqst




      subroutine pofv(p, v, vzero, bzero, bprim)

!     Birch equation of state p(V)

      implicit none
      integer, parameter  :: REAL64 = selected_real_kind(12)

!     input

      real(REAL64), intent(in)    ::  v, vzero, bzero, bprim

!     output

      real(REAL64)                ::  p

!     other variables

      real(REAL64)                ::  x, csi, f

!     constants

      real(REAL64), parameter    ::  UM = 1.0_REAL64

      x = vzero/v
      csi = 3 - 3*bprim/(4.0*UM)
      f = (x**((2*UM)/(3*UM))-UM) / 2
      p = 3*bzero*f*((UM+2*f)**((5*UM)/(2*UM)))*(UM-2*csi*f)

      return
      end subroutine pofv



      subroutine bofv(b, v, vzero, bzero, bprim)

!     Birch equation of state B(V)

      implicit none
      integer, parameter  :: REAL64 = selected_real_kind(12)

!     input

      real(REAL64), intent(in)    ::  v, vzero, bzero, bprim

!     output

      real(REAL64)                ::  b

!     other variables

      real(REAL64)                ::  x, csi, x23

!     constants

      real(REAL64), parameter    ::  UM = 1.0_REAL64

      x = vzero/v
      csi = 3 - 3*bprim/(4.0*UM)
      x23 = x**((2*UM)/(3*UM))
      b = -bzero*(x**((5*UM)/(3*UM)))*                                      &
     &              (5-7*x23+csi*(5+x23*(-14+x23*9))) / (2*UM)

      return
      end subroutine bofv



      subroutine eofv(e, v, ezero, vzero, bzero, bprim)

!     Birch equation of state E(V)

      implicit none
      integer, parameter  :: REAL64 = selected_real_kind(12)

!     input

      real(REAL64), intent(in)    ::  v, ezero, vzero, bzero, bprim

!     output

      real(REAL64)                ::  e

!     other variables

      real(REAL64)                ::  x, csi, x23

!     constants

      real(REAL64), parameter    ::  UM = 1.0_REAL64

      x = vzero/v
      csi = 3 - 3*bprim/(4.0*UM)
      x23 = x**((2*UM)/(3*UM))
      e = ezero + ((9*UM)/(8*UM))*bzero*vzero*                             &
     &    (1+2*csi/(3*UM)+x23*(-2-2*csi+x23*(1+2*csi-2*csi*x23/(3*UM))))

      return
      end subroutine eofv



      subroutine vofp(v, p, vzero, bzero, bprim)

!     Birch equation of state V(p)

      implicit none
      integer, parameter  :: REAL64 = selected_real_kind(12)

!     input

      real(REAL64), intent(in)    ::  p, vzero, bzero, bprim

!     output

      real(REAL64)                ::  v

!     other variables

      real(REAL64)                ::  vguess, pguess, bguess, vold

!     constants

      real(REAL64), parameter   ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
      real(REAL64), parameter   ::  EPS = 1.0E-10_REAL64

!     counters

      integer                     ::  ic, ic2

!     guess from murnaghan

       vguess = vzero*(bzero/(bzero+bprim*p))**(UM/bprim)
       call pofv(pguess, vguess, vzero, bzero, bprim)

!      newton step

       do ic = 1,50
         call bofv(bguess, vguess, vzero, bzero, bprim)
         vold = vguess
         vguess = vguess*(UM+(pguess-p)/bguess)

!        tries to avoid divergence of newton step

         if(vguess < ZERO) vguess = abs(vold) / 2

         do ic2 = 1,20
           call pofv(pguess, vguess, vzero, bzero, bprim)
           if(pguess < ZERO .and. vguess < vzero) then
             vguess = min(2*vguess,0.9*vzero)
           else
             exit
           endif
         enddo
         if(pguess < ZERO .and. vguess < vzero) STOP 'vofp2'

         if(abs(p-pguess) < EPS) exit
       enddo
       if(abs(p-pguess) >= EPS) STOP 'vofp'

       v = vguess

       return
       end subroutine vofp



      subroutine hofp(h, p, ezero, vzero, bzero, bprim)

!     Birch equation of state H(p)

      implicit none

      integer, parameter  :: REAL64 = selected_real_kind(12)

!     input

      real(REAL64), intent(in)    ::  p, ezero, vzero, bzero, bprim

!     output

      real(REAL64)                ::  h

!     other variables

      real(REAL64)                ::  v, e

      call vofp(v, p, vzero, bzero, bprim)
      call eofv(e, v, ezero, vzero, bzero, bprim)
      h = e + p*v

      return
      end subroutine hofp



      function f1(v, vzero)

!     birch equation of state (fit)

      implicit none
      integer, parameter  :: REAL64 = selected_real_kind(12)

!     input

      real(REAL64), intent(in)    ::  v, vzero

!     value

      real(REAL64)                ::  f1

!     other variables

      real(REAL64)                ::  x, x23

!     constants

      real(REAL64), parameter    ::  UM = 1.0_REAL64

      x = vzero/v
      x23 = x**((2*UM)/(3*UM))
!      f1=1.125+x23*(-2.25+x23*1.125)
      f1 = 9*(1 + x23*(-2+x23)) / (8*UM)

      return
      end function f1



      function f2(v, vzero)

!     birch equation of state (fit)

      implicit none
      integer, parameter  :: REAL64 = selected_real_kind(12)

!     input

      real(REAL64), intent(in)    ::  v, vzero

!     value

      real(REAL64)                ::  f2

!     other variables

      real(REAL64)                ::  x, x23

!     constants

      real(REAL64), parameter    ::  UM = 1.0_REAL64

      x = vzero/v
      x23 = x**((2*UM)/(3*UM))
!      f2 = 0.75+x23*(-2.25+x23*(2.25-x23*0.75))
      f2 = 3*(1 + x23*(-3 + x23*(3-x23))) / (4*UM)

      return
      end function f2



      subroutine brent(f,ax,bx,cx,ez,xmin,bz,bp,tol)

      implicit none
      integer, parameter  :: REAL64 = selected_real_kind(12)

!     input

      real(REAL64), intent(in)    ::  ax, bx, cx
      real(REAL64), intent(in)    ::  tol

      real(REAL64), external      ::  f

!     output

      real(REAL64), intent(out)   ::  xmin


!     input and output

      real(REAL64), intent(inout)    ::  ez, bz, bp

!     local

      real(REAL64)      ::  a, b
      real(REAL64)      ::  v, w, x
      real(REAL64)      ::  fv, fw, fx
      real(REAL64)      ::  e
      real(REAL64)      ::  xm
      real(REAL64)      ::  p, q, r
      real(REAL64)      ::  u, d
      real(REAL64)      ::  fu, etemp
      real(REAL64)      ::  tol1, tol2

      logical      ::  ljump, lfail

!     parameters

      integer, parameter        ::  ITMAX = 100
      real(REAL64), parameter   ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
      real(REAL64), parameter   ::  R5 = sqrt(5*UM)
      real(REAL64), parameter   ::  CGOLD = (R5-UM)/(R5+UM)
      real(REAL64), parameter   ::  ZEPS = 1.0e-10_REAL64

!     counters

      integer    ::  iter

      a = min(ax,cx)
      b = max(ax,cx)
      v = bx
      w = v
      x = v
      e = ZERO
      fx = f(x,ez,bz,bp)
      fv = fx
      fw = fx
      lfail  =  .true.

      do iter = 1,ITMAX
        xm = (a+b) / 2
        tol1 = TOL*abs(x)+ZEPS
        tol2 = 2*tol1
        if(abs(x-xm) <= (tol2-(b-a)/2)) then
          lfail  =  .false.

          exit

        endif
        ljump  =  .true.
        if(abs(e) > tol1) then
          r = (x-w)*(fx-fv)
          q = (x-v)*(fx-fw)
          p = (x-v)*q-(x-w)*r
          q = 2*(q-r)
          if(q > ZERO) p = -p
          q = abs(q)
          etemp = e
          e = d
          if(abs(p) >= abs(q*etemp/2) .or. p <= q*(a-x) .or.            &
     &          p >= q*(b-x)) then
            ljump  =  .true.
          else
            d = p/q
            u = x+d
            if(u-a < tol2 .or. b-u < tol2) d = sign(tol1,xm-x)
            ljump  =  .false.
          endif
        endif

        if(ljump) then
          if(x >= xm) then
            e = a-x
          else
            e = b-x
          endif
          d = CGOLD*e
        endif

        if(abs(d) >= tol1) then
          u = x+d
        else
          u = x+sign(tol1,d)
        endif
        fu = f(u,ez,bz,bp)
        if(fu <= fx) then
          if(u >= x) then
            a = x
          else
            b = x
          endif
          v = w
          fv = fw
          w = x
          fw = fx
          x = u
          fx = fu
        else
          if(u < x) then
            a = u
          else
            b = u
          endif
          if(fu <= fw .or. w == x) then
            v = w
            fv = fw
            w = u
            fw = fu
          else if(fu <= fv .or. v == x .or. v == w) then
            v = u
            fv = fu
          endif
        endif
      enddo

      if(lfail) stop 'brent exceed maximum iterations.'

      xmin = x
      fx = f(x,ez,bz,bp)

      return
      end subroutine brent



      function fbirch(vzero, ezero, bzero, bprim)

      implicit none
      integer, parameter  :: REAL64 = selected_real_kind(12)

!     input

      real(REAL64), intent(in)    ::  vzero

!     output

      real(REAL64), intent(out)   ::  ezero,bzero,bprim

!     value

      real(REAL64)                ::  fbirch

!     common

      real(REAL64)                ::  xar, yar, xscl, yscl, fn
      integer                     ::  npt
      common/params/ xar(40), yar(40), xscl(40), yscl(40), fn, npt

!     other variables

      real(REAL64)                ::  a(3,3), b(3), c(3)
      real(REAL64)                ::  csi, sumx, yfit

!     external functions

      real(REAL64), external      ::  f1, f2

!     constants

      real(REAL64), parameter    :: ZERO = 0.0_REAL64, UM = 1.0_REAL64

!     counters

      integer                     ::  i, j, k

      do j = 1,3
        b(j) = ZERO
        do k = 1,3
          a(k,j) = ZERO
        enddo
      enddo

      do i = 1,npt
        c(1) = UM
        c(2) = f1(xar(i),vzero)
        c(3) = f2(xar(i),vzero)
        do j = 1,3
          b(j) = b(j) + yscl(i)*c(j)
          do k = 1,3
            a(k,j) = a(k,j) + c(k)*c(j)
          enddo
        enddo
      enddo

      call gaussj(a,3,3,b,1,1)

      ezero = b(1)
      bzero = b(2)/vzero
      csi = b(3)/b(2)
      bprim = 4*(UM-csi/(3*UM))

      sumx = ZERO
      do i = 1,npt
        yfit = b(1) + b(2)*f1(xar(i),vzero) + b(3)*f2(xar(i),vzero)
        sumx = sumx + (yscl(i)-yfit)*(yscl(i)-yfit)
      enddo
      fbirch = sumx

      return
      end function fbirch

      subroutine gaussj(a, n, np, b, m, mp)

      implicit none
      integer, parameter  :: REAL64 = selected_real_kind(12)
      integer, parameter  :: NMAX = 50

!     input

      integer, intent(in)         ::  n
      integer, intent(in)         ::  m
      integer, intent(in)         ::  np
      integer, intent(in)         ::  mp

!     input and output

      real(REAL64), intent(inout) ::  a(np,np)
      real(REAL64), intent(inout) ::  b(np,mp)

!     local

      integer              ::  ipiv(nmax), indxr(nmax), indxc(nmax)
      integer              ::  icol, irow
      real(REAL64)         ::  dum, big, pivinv

!      constants

       real(REAL64), parameter    :: ZERO = 0.0_REAL64, UM = 1.0_REAL64

!     counters

      integer                     ::  i, j, k, l, ll

      do j = 1,n
        ipiv(j) = 0
      enddo

      do i = 1,n
        big = ZERO

        do j = 1,n
          if(ipiv(j) /= 1)then
            do k = 1,n
              if (ipiv(k) == 0) then
                if (abs(a(j,k)) >= big)then
                  big = abs(a(j,k))
                  irow = j
                  icol = k
                endif
              else if (ipiv(k) > 1) then
                stop 'singular matrix'
              endif
            enddo
          endif
        enddo

        ipiv(icol) = ipiv(icol)+1
        if (irow /= icol) then
          do l = 1,n
            dum = a(irow,l)
            a(irow,l) = a(icol,l)
            a(icol,l) = dum
          enddo
          do l = 1,m
            dum = b(irow,l)
            b(irow,l) = b(icol,l)
            b(icol,l) = dum
          enddo
        endif

        indxr(i) = irow
        indxc(i) = icol
        if (a(icol,icol) == ZERO) stop 'singular matrix.'
        pivinv = UM/a(icol,icol)
        a(icol,icol) = UM
        do l = 1,n
          a(icol,l) = a(icol,l)*pivinv
        enddo
        do l = 1,m
          b(icol,l) = b(icol,l)*pivinv
        enddo

        do ll = 1,n
          if(ll /= icol)then
            dum = a(ll,icol)
            a(ll,icol) = ZERO
            do l = 1,n
              a(ll,l) = a(ll,l) - a(icol,l)*dum
            enddo
            do l = 1,m
              b(ll,l) = b(ll,l) - b(icol,l)*dum
            enddo
          endif
        enddo

      enddo

      do l = n,1,-1
        if(indxr(l) /= indxc(l))then
          do k = 1,n
            dum = a(k,indxr(l))
            a(k,indxr(l)) = a(k,indxc(l))
            a(k,indxc(l)) = dum
          enddo
        endif
      enddo

      return
      end subroutine gaussj


