c      program epc
c  vax version
c  electroproduction yields of nucleons and pions
c  written by j.s. oconnell and j.w. lightbody, jr.
c  national bureau of standards
c  april 1988
c  transverse scaling region added



c Modified by Carlos Ayerbe to make it a function to be read by 
c C++ code lines modified or comented has marked with (CA)
      real*8 function epc_func(ebeam, z1, n1, partID, p, thp)



      implicit real*8 (a-h,o-z)


      common/sg/ia
      common/qd/qdf
      common/del/ip
      common/m/z,n

      character part*3,infile*80,outfile*80
      logical fileok

      data am/938.28/,amp/139.6/

c redefining the variables with the code nomenclature (CA)
      integer partID !(CA)
      real*8 ebeam
      integer n1 !number of neutrons (NOT NUCLEONS!!)
      integer z1 !number of protons
      real*8 p
      real*8 thp

c      print *,"beam:", ebeam,"z:",z1, " n:", n1, "partID:", partID  !(CA)
c      print *, "p:", p, "angle: ", thp  !(CA)


	if(partID.eq.1)then
	   part = 'p'
	elseif(partID.eq.-1)then
	   part = 'n'
	elseif(partID.eq.2)then
	   part = 'pi+'
	elseif(partID.eq.-2)then
	   part = 'pi-'
	elseif(partID.eq.0)then
	   part = 'pi0'
	endif

 

      z = z1
      n = n1
      e1 = ebeam

c      print *,"*********************" !(CA)
c      print *,"beam:", e1,"z:",z, " n:", n !(CA)


c p and thp are defined the same (CA)



c I comment these lines which read an external file (CA)
c      fileok = .false.
c      do while(.not.fileok)
cccc         call QQPARM ( infile, i )
cccc         if ( i .EQ. 0 ) then
c         if ( iargc() .GE. 1 ) then
c            call getarg ( 1, infile )
c            write ( 6, '(''input file: '',A)' ) infile
c         else
c            write(6,'('' input file: '',$)')
c            read(5,'(a80)') infile
c         endif
c         open(1,file=infile,access='sequential',status='old',err=2000)
c         fileok = .true.
c2000     continue
c         if (.not.fileok)
c     >     write(6,'('' Error opening input file'')')
c      enddo

c      call filename(infile,outfile,'res',.false.)


c COMMENTED TO REMOVE THE OUTPUT FILE (CA)
c      outfile ='deut_p.out'
c      write(6,'('' results on output file: '',a)') outfile
c      open(2,file=outfile,access='sequential',status='unknown')


      pi = acos(-1.)

c I comment these lines which read an external file (CA)
c this need to be tweaked better since the case 0 stops
c the program, but I sent the value 0 as pi0, so 
c there in no part 0 case

c  input of z,n
c      read(1,*) z,n
c input of particle type [p,n,pi+,pi-,pi0, 0 = stop]
c2     read(1,'(a3)') part

      ia = z+n

c COMMENTED TO REMOVE THE OUTPUT FILE (CA)
c      if (ia.eq.2.or.ia.eq.3.or.ia.eq.4
c     >   .or.ia.eq.12.or.ia.eq.16) then
c         write(2,'('' using a = '',i3,'' spectral function'')') ia
c      else if (ia.eq.1) then
c         write(2,'('' nucleon case'')')
c      else
c         write(2,'('' no specific spectral  for this nucleus'')')
c         write(2,'('' a = 16 spectral  will be used'')')
c      endif
c      write(2,'('' '')')
c      write(2,'('' Z='',I3,'', A='',I3,'',  (e,'',a3,'') cross ''
c     >            ''sections in ub/mev/c-sr'')') NINT(z), ia, part
c      write(2,'('' '')')
c      write(2,'('' Eel    mom   angle   kin.E   d2qd    d2qf    d2del
c     >  d2sc     total     totale'')')
c      write(2,'('' MeV   MeV/c   deg     MeV   <-------------ub/MeV-
c     >sr--------------->  ub/MeV-sr'')')

c  'an' is effective number of nucleons for one pion production
      if (part.eq.'p') then
         an = n/3.+2.*z/3.
         ip = 1
      else if (part.eq.'n') then
         an = z/3.+2.*n/3.
         ip = -1
      else if (part.eq.'pi+') then
         an = z/3.
         ip = 2
      else if (part.eq.'pi-') then
         an = n/3.
         ip = 2
      else if (part.eq.'pi0') then
         an = 2.*(n+z)/3.
         ip = 0
      else
         stop
      endif
c      write(6,*)an,ip
      if (ia.eq.1) then
      else if (abs(ip).eq.1) then
         if (ia.gt.1.and.ia.lt.5) then
            dlf = ia
         else
            dlf = 7.
         endif
         al = dlf
         qdf = al*n*z/float(ia)
      else
      endif

c start of looping over electron energy e [MeV] and momenta,angles of proton
c or pion (p [MeV/c],thp [deg]). Now coming from the wrapper
c 1    read(1,*,err=1000,end=1001) e1,p,thp
c      if (e1.le.0.) goto 2
      if (e1.le.0.) goto 1001
      th = thp*pi/180.
      if (abs(ip).eq.1) then
c        print *, " p: ", p, "am:", am
         e = sqrt(p**2+am**2)
         tp = p**2/(e+am)
c        print *, "tp: ", tp, " p: ", p, "e: ", e, "am:", am
         aj = p/e
         if (ia.eq.1) then
            d2qd = 0.
            d2qf = 0.
         else if (ia.gt.1) then
            call dep(e1,tp,th,ip,d2qd)
c           print *, "--->d2qd ", d2qd, "aj", aj
            d2qd = d2qd*aj
c           print *, "***d2qd ", d2qd, "aj", aj
      
c            write(6,*)d2qd
            call ep(e1,tp,th,d2qf)
            d2qf = d2qf*aj
         endif
      else if (abs(ip).eq.2.or.ip.eq.0) then
         e = sqrt(p**2+amp**2)
         tp = p**2/(e+amp)
         aj = p/e
         d2qd = 0.
         d2qf = 0.
      else
         stop
      endif
      if (p.lt.500.) then
         d2sc = 0.
         call delta(e1,tp,th,d2del)
         d2del = an*d2del*aj
      else
         d2del = 0.
         if (abs(ip).eq.2.or.ip.eq.0) then
c            write(6,*)'calling s2pi'
            call s2pi(2,e1,tp,th,d2sc1)
            call s2pi(-2,e1,tp,th,d2sc2)
            if (part.eq.'pi+') then
               d2sc = z*d2sc1+n*d2sc2
            else if (part.eq.'pi-') then
               d2sc = n*d2sc1+z*d2sc2
            else if (part.eq.'pi0') then
               d2sc = ia*d2sc1
            else
               stop
            endif
         else if (abs(ip).eq.1) then
            call s2pi(1,e1,tp,th,d2sc1)
            d2sc = ia*d2sc1
         endif
      endif


c I added this condition, because sometimes the value d2sc return NaN

      if (d2sc.ne.d2sc) then
         total = d2qd+d2qf+d2del
         totale = total/aj
      else
         total = d2qd+d2qf+d2del+d2sc
         totale = total/aj
      endif

c      if (totale.eq.0)

c      print *,"d2qd", d2qd, "d2qf", d2qf, "d2del", d2del,"d2sc", d2sc ! (CA)


c COMMENTED TO REMOVE THE OUTPUT FILE (CA)
c      write(2,'(f7.1,1x,f7.1,1x,f6.2,1x,f7.2,1x,1pe7.1,1x,e7.1,1x,
c     >                e7.1,1x,e7.1,1x,e9.3,1x,e9.3,1x,e9.3)')
c     >                     e1,p,thp,tp,d2qd,d2qf,d2del,d2sc,total,totale
c     >,aj

c not sure of comment this, but it is sending to a commented line (CA)
c      go to 1

c1000  write(6,*)'error reading file'
1001  continue

c      close(1)
c      close(2)

c     print *,"Kinetic Energy", tp ! (CA)
        
      epc_func = totale

c     print *,"coming back (END) ", totale ! (CA)

      return ! because this is a function (CA)


      end

      subroutine filename(string,file,ext,number)

c... this subroutine creates the string file by glueing 'ext' after the
c... last occurence of '.' in the string 'string'.

      implicit none

      character string*(*),file*(*),ext*(*),tmp*80,detnr*3
      integer*4 lastpoint,nrchar,extlength,i,digpoint
      logical   stop,number

c... determine the position in the string 'string' where the last '.'
c... occures.

      nrchar = len(string)
      do i = 1, nrchar
         tmp(i:i) = string(nrchar-i+1:nrchar-i+1)
      enddo
      lastpoint   = index(tmp,'.')

      if (lastpoint .le. 0) then
         lastpoint = index(string,' ')
      else
         lastpoint = nrchar - lastpoint + 1
         stop      = .false.
         digpoint  = 0
         if (number) then
            do 20 i = lastpoint+1, nrchar
               if (.not.stop) then
                  if (ichar(string(i:i)) .ge. ichar('0') .and.
     +                          ichar(string(i:i)) .le. ichar('9')) then
                     digpoint                 = i - lastpoint
                     detnr(digpoint:digpoint) = string(i:i)
                     stop                     = .false.
                  else
                     stop                     = .true.
                  endif
               endif
20          continue
         endif
      endif

c... create file name

      if (index(ext,' ') .ne. 0) then
         extlength                            = index(ext,' ')-1
      else
         extlength                            = len(ext)
      endif
      file(1:lastpoint-1)                     = string(1:lastpoint-1)
      file(lastpoint:lastpoint)               = '.'
      if (digpoint .ne. 0)
     +    file(lastpoint+1:lastpoint+digpoint) = detnr(1:digpoint)
      file(lastpoint+digpoint+1:lastpoint+digpoint+extlength)
     +                                         = ext(1:extlength)


      return
      end

      subroutine vtp(amt,am1,ei,w0,tp,th,gn)

c  tiator-wright virtual photon spectrum
c  phys. rev. c26,2349(1982) and nuc. phys. a379,407(1982)

      implicit real*8 (a-h,o-z)

      data ame/.511/

      pi = acos(-1.)
      ef0 = ei-w0
      if (ef0.ge.ame) then
         aki = sqrt(ei**2-ame**2)
         akf0 = sqrt(ef0**2-ame**2)
         akp = sqrt(tp**2+2.*am1*tp)
         ep = tp+am1
         ar = ei+amt-ep
         br = ef0*(akp*cos(th)-aki)/akf0
c     brp = (akf0/ef0)**2*br
         a = ame**2-ei*ef0
         b = aki*akf0
         d = -ame**2*br*(ei/ef0-1.)/ar
         ap = a-d
         bp = b+d
         an1 = 1./137./2./pi*w0**2/aki**2
         apb = -ame**2*(aki-akf0)**2/(ame**2+ei*ef0+aki*akf0)
         an1 = an1*b/bp*(ar+br)/(ar-ap/bp*br)
         an2 = 1.-2.*a/w0**2
         an4 = ((ap-bp)*(ar+br)/apb/(ar-br))
         if (an4.gt.0.) then
            an2 = an2*log(an4)
            an3 = -4.*b/w0**2
            ane = an1*(an2+an3)
            d0 = amt+ei-ep+ef0/akf0*(akp*cos(th)-aki)
            r = (amt+w0-ep/akp*w0*cos(th))/d0
            gn = ane*r/w0
            if (gn.lt.0.)gn = 0.
            return
         endif
      endif
      gn = 0.

      return
      end

      subroutine dep(e1,tp,th,ip,d2qd)

c  quasi-deuteron cross section

      implicit real*8 (a-h,o-z)

      common/sg/ia
      common/qd/qdf

      data am/939./,amd/1876./

      if (ia.ne.1) then
         pn = sqrt(tp**2+2.*am*tp)
         call kine(amd,am,am,pn,th,w0,thc)
         if (w0.gt.0. .and. w0.lt.e1) then
            w0g = w0/1000.
            call sigd(w0g,thc,ip,dsqd)
            call part(amd,am,am,pn,th,ajt,ajw)
c     cross section in ub/mev-sr
            call vtp(amd,am,e1,w0,tp,th,phi)
            d2qd = qdf*phi*dsqd*ajw*ajt
c           print *, "*-*-*d2qd ", d2qd
            return
         endif
      endif
      d2qd = 0.

      return
      end

      subroutine sigd(e,th,ip,dsqd)

c  deuteron cross section
c  based on fit of thorlacius & fearing
c  phys. rev. c33,1830(1986)
c  photon energy range 10 - 625 mev
c  e[gev] in lab system
c  th[rad] & dsqd[ub/sr] in center-of-momentum system

      implicit real*8 (a-h,o-z)

      dimension c0(8),c1(4),c2(4),c3(4),c4(4)
      dimension a(0:4),b(4,4)

      data c0/2.61e2,-1.10e2,2.46e1,-1.71e1,5.76e0,
     #    -2.05e0,2.67e-1,1.13e2/
      data c1/1.68e1,-4.66e1,2.56e0,-4.72e0/
      data c2/-2.03e2,-8.12e1,-4.05e0,-5.99e0/
      data c3/-1.77e1,-3.74e1,-5.07e-1,-5.40e0/
      data c4/-2.05e0,-7.05e0,9.40e-1,-2.05e0/

      x = cos(th)
      if (e.le..625) then
c  test for neutron
         x = ip*x
c  coeficients
         a(0) = c0(1)*exp(c0(2)*e)+c0(3)*exp(c0(4)*e)
         a(0) = a(0)+(c0(5)+c0(6)*e)/(1.+c0(8)*(e-c0(7))**2)
         dsqd = a(0)*p(0,x)
         do l = 1,4
            b(1,l) = c1(l)
            b(2,l) = c2(l)
            b(3,l) = c3(l)
            b(4,l) = c4(l)
         enddo
         do l = 1,4
            a(l) = b(l,1)*exp(b(l,2)*e)+ b(l,3)*exp(b(l,4)*e)
            dsqd = dsqd+a(l)*p(l,x)
         enddo
      else if (e.lt..700) then
         dsqd = .3
      else if (e.lt..800) then
         dsqd = .15
      else if (e.lt..900) then
         dsqd = .1
      else
         dsqd = 55./(e-.350)
      endif

      return
      end

      real*8 function p(l,x)

c  legendre polynomials

      implicit real*8 (a-h,o-z)

      if (l.eq.0) then
         p = 1.
      else if (l.eq.1) then
         p = x
      else if (l.eq.2) then
         p = .5*(3.*x**2-1.)
      else if (l.eq.3) then
         p = .5*(5.*x**3-3.*x)
      else if (l.eq.4) then
         p = 1./8.*(35.*x**4-30.*x**2+3.)
      else
         p = 0.
      endif

      return
      end

      subroutine delta(e1,tp,th,d2del)

c  photoproduction of nucleons and pions via delta

      implicit real*8 (a-h,o-z)

      common/del/ip

      data am/939./,amp/139./

      if (abs(ip).eq.1) then
         am1 = am
         am2 = amp
      else
         am1 = amp
         am2 = am
      endif
      ep = tp+am1
      pn = sqrt(ep**2-am1**2)
      call kine(am,am1,am2,pn,th,w,tc)
      if (w.gt.0. .and. w.lt.e1) then
         call part(am,am1,am2,pn,th,ajt,ajw)
         call sigma(w,tc,dsigg)
         call vtp(am,am1,e1,w,tp,th,phi)
         d2del = phi*dsigg*ajt
c cross section in ub/mev-sr
         return
      endif
      d2del = 0.

      return
      end

      subroutine part(amt,am1,am2,pn,tn,ajt,ajw)

c  partial derivatives

      implicit real*8 (a-h,o-z)

      pi = acos(-1.)
      dt = pi/50.
      dp = 10.
c  angle
      tnp = tn+dt
      tnm = tn-dt
      call kine(amt,am1,am2,pn,tnp,w,tcp)
      call kine(amt,am1,am2,pn,tnm,w,tcm)
      ajt = (cos(tcp)-cos(tcm))/(cos(tnp)-cos(tnm))
      ajt = abs(ajt)
c  energy
      pnp = pn+dp
      pnm = pn-dp
      call kine(amt,am1,am2,pnp,tn,wp,tc)
      call kine(amt,am1,am2,pnm,tn,wm,tc)
      ajw = (wp-wm)/(pnp-pnm)
      ajw = abs(ajw)

      return
      end

      subroutine kine(amt,am1,am2,pn,th,w,tc)

c  computes cm variables from lab variables

      implicit real*8 (a-h,o-z)

      ep = sqrt(pn**2+am1**2)
      pnt = pn*sin(th)
      pnl = pn*cos(th)
      anum = pn**2+am2**2-(amt-ep)**2
      den = 2.*(pnl+amt-ep)
      w = anum/den
      if (w.le.0.)w = 0.
c  invariant mass
      ww = sqrt(amt**2+2.*w*amt)
c  cm variables
      pct = pnt
      b = w/(amt+w)
      g = (w+amt)/ww
      pcl = g*(pnl-b*ep)
      pcs = pcl**2+pct**2
      pc = sqrt(pcs)
      cthc = pcl/pc
      tc = acos(cthc)

      return
      end

      subroutine sigma(e,thrcm,sigcm)

c  real photon cross section in delta region
c  microbarns per steradian

      implicit real*8 (a-h,o-z)

      gam = 100.
      pi = acos(-1.)
      if (e.gt.420.) then
        sigcm = (1.+420./e)*90./4./pi
      else
        sigcm = 360.*(5.-3.*cos(thrcm)**2)
        sigcm = sigcm/16./pi/(1.+(e-320)**2/gam**2)
      endif

      return
      end

      subroutine ep(e1,tp,thp,dsep)

c  electro proton production cross sections
      implicit real*8 (a-h,o-z)

      common/pco/ph(10),wph(10)

      data aml/.511/

      pi = acos(-1.d0)
      call gausab(10,ph,wph,0.d0,2.*pi,pi)
      ak = sqrt(e1**2-aml**2)
      call sep(ak,tp,thp,dsep)
      dsep = dsep*1.e4
c  cross section in ub/mev-sr

      end

      real*8 function dot(v,u)

      implicit real*8 (a-h,o-z)

      dimension v(3),u(3)

      dot = 0.
      do i = 1,3
         dot = dot+v(i)*u(i)
      enddo
      
      return
      end

      subroutine cross(v,u,w)

      implicit real*8 (a-h,o-z)

      dimension v(3),u(3),w(3)

      w(1) = v(2)*u(3)-v(3)*u(2)
      w(2) = v(3)*u(1)-v(1)*u(3)
      w(3) = v(1)*u(2)-v(2)*u(1)

      return
      end

      subroutine gausab(n,e,w,a,b,c)

      implicit real*8 (a-h,o-z)

      dimension e(*),w(*)

      data eps/1.d-16/

      if (a.ge.c.or.c.ge.b)stop
c           stops program if a, c, b are out of sequence
      pi = acos(-1.d0)
      al = (c*(a+b)-2*a*b)/(b-a)
      be = (a+b-2*c)/(b-a)
      m = (n+1)/2
      dn = n
      do i = 1,m
         di = i
         x = pi*(4.d0*(dn-di)+3.d0)/(4.d0*dn+2.d0)
         xn = (1.d0-(dn-1.d0)/(8.d0*dn*dn*dn))*cos(x)
         if (i.gt.n/2) xn = 0
         do iter = 1,10
            x = xn
            y1 = 1.d0
            y = x
            if (n.ge.2) then
               do j = 2,n
                  dj = j
                  y2 = y1
                  y1 = y
                  y = ((2.d0*dj-1.d0)*x*y1-(dj-1.d0)*y2)/dj
               enddo
            endif
            ys = dn*(x*y-y1)/(x*x-1.d0)
            h = -y/ys
            xn = x+h
            if (abs(h).lt.eps) exit
         enddo
         e(i) = (c+al*x)/(1.d0-be*x)
         e(n-i+1) = (c-al*x)/(1.d0+be*x)
         gew = 2.d0/((1.d0-x*x)*ys*ys)
         w(i) = gew*(al+be*c)/(1.d0-be*x)**2
         w(n-i+1) = gew*(al+be*c)/(1.d0+be*x)**2
      enddo

      return
      end

      subroutine vect(thp,the,phi,p,ak1,ak2)

c  cartesian components of electron and proton vectors

      implicit real*8 (a-h,o-z)

      common/v/ak1v(3),ak2v(3),qv(3),pv(3),pp(3)

      pv(1) = p*sin(thp)
      pv(2) = 0.
      pv(3) = p*cos(thp)
      ak1v(1) = 0.
      ak1v(2) = 0.
      ak1v(3) = ak1
      ak2v(1) = ak2*sin(the)*cos(phi)
      ak2v(2) = ak2*sin(the)*sin(phi)
      ak2v(3) = ak2*cos(the)
      qv(1) = ak1v(1)-ak2v(1)
      qv(2) = ak1v(2)-ak2v(2)
      qv(3) = ak1v(3)-ak2v(3)
      pp(1) = pv(1)-qv(1)
      pp(2) = pv(2)-qv(2)
      pp(3) = pv(3)-qv(3)

      return
      end

      real*8 function amag(v)

      implicit real*8 (a-h,o-z)

      dimension v(3)

      amag = 0.
      do i = 1,3
         amag = amag+v(i)**2
      enddo
      amag = sqrt(amag)

      return
      end

      subroutine lept(e1,e2,ak1,ak2,aml,qs,qus,the,v)

c  lepton factors for coincidence cross section

      implicit real*8 (a-h,o-z)

      dimension v(5)

      v(1) = (qus/qs)**2*(e1*e2+ak1*ak2*cos(the)+aml**2)
      x = ak1*ak2*sin(the)
      v(2) = x**2/qs+qus/2.
      v(3) = qus/qs*x/sqrt(qs)*(e1+e2)
      v(4) = x**2/qs
      v(5) = 0.

      return
      end

      subroutine d4s(ak1,ak2,the,p,pp,thqp,cphip,dsig)

c  fully differential cross section

      implicit real*8 (a-h,o-z)

      dimension v(5),w(5)

      data am/939./,aml/.511/,a/855./

      qs = ak1**2+ak2**2-2.*ak1*ak2*cos(the)
      e1 = sqrt(ak1**2+aml**2)
      e2 = sqrt(ak2**2+aml**2)
      qus = 2.*(e1*e2-ak1*ak2*cos(the)-aml**2)
      sm = 2.*(1.44)**2/qus**2*ak2/ak1
      ep = sqrt(am**2+p**2)
      ps = ep*p
      fns = 1./(1.+qus/a**2)**4
      call lept(e1,e2,ak1,ak2,aml,qs,qus,the,v)
      call form(qs,p,thqp,cphip,w)
      sum = 0.
      do i = 1,5
         sum = sum+v(i)*w(i)
      enddo
      dsig = sm*ps*fns*sum*sgsl(pp)
      return
      end

      subroutine sthe(d2s)

c  integral over electron polar angle

      implicit real*8 (a-h,o-z)

      common/s/ ak1,ak2,the,p,thp
      common/e/th1(12),wt1(12),th2(12),wt2(12)
      common/e1/th3(24),wt3(24)

      d2s1 = 0.
      do i = 1,12
         the = th1(i)
         call sphi(d3s)

         d2s1 = d2s1+d3s*wt1(i)*sin(the)
      enddo
      d2s2 = 0.
      do i = 1,12
         the = th2(i)
         call sphi(d3s)
         d2s2 = d2s2+d3s*wt2(i)*sin(the)
      enddo
      d2s3 = 0.
      do i = 1,24
         the = th3(i)
         call sphi(d3s)
         d2s3 = d2s3+d3s*wt3(i)*sin(the)
      enddo
      d2s = d2s1+d2s2+d2s3

      return
      end

      subroutine sphi(d3s)

c  integrate over electron azimuthal angle

      implicit real*8 (a-h,o-z)

      common/s/ ak1,ak2,the,p,thp
      common/v/ak1v(3),ak2v(3),qv(3),pv(3),pp(3)
      common/pco/ph(10),wph(10)

      dimension qxp(3),ak1x2(3)

      d3s = 0.
      do i = 1,10
         phi = ph(i)
         call vect(thp,the,phi,p,ak1,ak2)
         call cross(qv,pv,qxp)
         call cross(ak1v,ak2v,ak1x2)
c  proton theta
         cthep = dot(pv,qv)/amag(pv)/amag(qv)
         thqp = acos(cthep)
c  proton phi
         cphip = dot(qxp,ak1x2)
         if (cphip.eq.0.) then
            cphip = 1.
         else
            cphip = cphip/amag(qxp)/amag(ak1x2)
         endif
         ppm = amag(pp)
         call d4s(ak1,ak2,the,p,ppm,thqp,cphip,dsig)

         d3s = d3s+dsig*wph(i)
      enddo

      return
      end

      subroutine form(qs,p,thqp,cphip,w)

c  nuclear form factors

      implicit real*8 (a-h,o-z)

      common/m/zz,nn
      common/del/ip

      dimension w(5)

      data am/939./,up/2.79/,un/-1.91/

      if (ip.eq.1) then
         z = zz
         n = 0.
      else if (ip.eq.-1) then
         z = 0.
         n = nn
      else
         z = 0.
         n = 0.
      endif
      y = p/am*sin(thqp)
      w(1) = z
      w(2) = z*y**2
      w(2) = w(2)+(z*up**2+n*un**2)*qs/2./am**2
      w(3) = -2.*z*y*cphip
      w(4) = z*y**2*(2.*cphip**2-1.)
      w(5) = 0.

      return
      end

      subroutine sep(ak,tp,thpp,d2s)

      implicit real*8 (a-h,o-z)

      common/s/ak1,ak2,the,p,thp
      common/e/th1(12),wt1(12),th2(12),wt2(12)
      common/e1/th3(24),wt3(24)

      data am/939./,aml/.511/,be/16./

      pi = acos(-1.d0)
      thp = thpp
      ak1 = ak
      ak2 = ak1-tp-be
c  gaussian points for the
      themax = aml*(ak1-ak2)/ak1/ak2
      call gausab(12,th1,wt1,0.d0,2.*themax,themax)
      call gausab(12,th2,wt2,2.*themax,100.*themax,10.*themax)
      a3 = 100.*themax
      c3 = a3+(pi-a3)/10.
      call gausab(24,th3,wt3,a3,pi,c3)
      p = sqrt(tp**2+2.*tp*am)
      call sthe(d2s)
      if (ak2.le.0.)d2s = 0.

      return
      end

      real*8 function sgsl(p)

      implicit real*8 (a-h,o-z)

c  p integral over sgsl normalized to 1/4pi
      common/sg/ia
      common/qd/qdf

      if (ia.eq.2) then
c  begin 2-h
         pp = p/197.3
         sgs = 3.697-7.428*pp-2.257*pp**2
         sgs = sgs+3.618*pp**3-1.377*pp**4+.221*pp**5-.013*pp**6
         if (sgs.lt.-293.)go to 1
         sgs = exp(sgs)
         sgs = sgs/.18825/4./3.1416/(197.3)**3
         sgsl = sgs/1.
      else if (ia.eq.3) then
c  begin 3-he
         if (-(p/33)**2.lt.-293.)go to 1
         sgs = 2.4101e-6*exp(-p/33)
         sgs = sgs-1.4461e-6*exp(-(p/33)**2)
         sgs = sgs+1.6871e-10*exp(-(p/493)**2)
         sgsl = sgs/2.
      else if (ia.eq.4) then
c   begin 4-he
         if (-(p/113.24)**2.lt.-293.)go to 1
         sgs = 1.39066e-6*exp(-(p/113.24)**2)
         sgs = sgs+3.96476e-9*exp(-(p/390.75)**2)
         sgsl = sgs/2.
         sgsl = sgsl/2./3.1416
      else if (ia.eq.12) then
c  begin 12-c
         if (-(p/127)**2.lt.-293.)go to 1
         sgs = 1.7052e-7*(1.+(p/127)**2)*exp(-(p/127)**2)
         sgs = sgs+1.7052e-9*exp(-(p/493)**2)
         sgsl = sgs/6.
      else
c  begin 16-o
         if (-(p/120)**2.lt.-293.)go to 1
         sgs = 3.0124e-7*(1.+(p/120)**2)*exp(-(p/120)**2)
         sgs = sgs+1.1296e-9*exp(-(p/493)**2)
         sgsl = sgs/8.
      endif
      return
    1 sgsl = 0.

      return
      end

      subroutine s2pi(ip,e1,tp,th,d2sc)

c  integral over scaling cross section

      implicit real*8 (a-h,o-z)

      data am/939./,amp/139./

      if (abs(ip).eq.1) then
c  one pion thr
         ap = am
         am2 = amp
      else if (ip.eq.2.or.ip.eq.0) then
c  one pion thr
         ap = amp
         am2 = am
      else if (ip.eq.-2) then
c  two pion thr
         ap = amp
         am2 = am+amp
      else
         stop
      endif
      p = sqrt(tp**2+2.*ap*tp)
      e = tp+ap
      call kine(am,ap,am2,p,th,thr,tc)
c      write(6,*)'thr is',thr
      if (e1.gt.thr) then
         dw = (e1-thr)/20.
         sum = 0.
         do i = 1,20
            w = thr+(float(i)-.5)*dw
            call vtp(am,amp,e1,w,tp,th,gn)
            call wiser(w/1.e3,p/1.e3,th,f)
            sum = sum+gn*f*dw
         enddo
c         write(6,*)'sum=',sum
         d2sc = sum*p**2/e*1.e-6
         return
      endif
      d2sc = 0.

      return
      end

      subroutine wiser(w,p,th,f)

c  invariant inclusive cross section
c  units in gev

      implicit real*8 (a-h,o-z)

      common/del/ip

      dimension a(7),b(7),c(7)

      data a/5.66e2,8.29e2,1.79,2.10,-5.49,-1.73,0./
      data b/4.86e2,1.15e2,1.77,2.18,-5.23,-1.82,0./
      data c/1.33e5,5.69e4,1.41,0.72,-6.77,1.90,-1.17e-2/
      data am/.939/,amp/.139/

      if (abs(ip).eq.1) then
         ap = am
      else if (abs(ip).eq.2.or.ip.eq.0) then
         ap = amp
      else
         stop
      endif
      e = sqrt(p**2+ap**2)
c  mandelstam variables
      s = 2.*w*am+am**2
c     t = -2.*w*e+2.*w*p*cos(th) +ap**2
      u = am**2-2.*am*e+ap**2
c  fitting variables
      pt = p*sin(th)
      aml = sqrt(pt**2+ap**2)
      call fxr(w,p,e,th,xr)
c  fitted
      if (ip.eq.2.or.ip.eq.0) then
         x1 = a(1)+a(2)/sqrt(s)
         x2 = (1.-xr+a(3)**2/s)**a(4)
         x3 = exp(a(5)*aml)
         x4 = exp(a(6)*pt**2/e)
         f = x1*x2*x3*x4
      else if (ip.eq.-2) then
         x1 = b(1)+b(2)/sqrt(s)
         x2 = (1.-xr+b(3)**2/s)**b(4)
         x3 = exp(b(5)*aml)
         x4 = exp(b(6)*pt**2/e)
         f = x1*x2*x3*x4
      else if (abs(ip).eq.1) then
         x1 = c(1)+c(2)/sqrt(s)
         x2 = (1.-xr+c(3)**2/s)**c(4)
         x3 = exp(c(5)*aml)
         x4 = 1./(1.+abs(u))**(c(6)+c(7)*s)
         f = x1*x2*x3*x4
      else
         stop
      endif

      return
      end

      subroutine fxr(w,p,e,th,xr)

c  computes ratio of cm particle momentum to photon momentum
c  gev units

      implicit real*8 (a-h,o-z)

      data am/.939/

      pt = p*sin(th)
      pl = p*cos(th)
c  lorentz transformation
      b = w/(w+am)
      d = sqrt(2.*w*am+am**2)
      g = (w+am)/d
      bg = b*g
c cm variables
      wc = g*w-bg*w
      plc = g*pl-bg*e
      pc = sqrt(pt**2+plc**2)
      xr = pc/wc

      return
      end


