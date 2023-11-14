!       **********************************************************************

subroutine f2pi_sub(E, x, theta_e, km1, km2, F2pi_ret)

!	PROGRAM TDIS
!       
!       COMPUTES THE LEADING CONTRIBUTION IN THE PION CLOUD MODEL TO
!       PROTON PRODUCTION FROM SU(2) MESON-BARYON CONFIGURATIONS OF
!       THE NUCLEON; AS SUCH, THE MAIN MECHANISM IS EXCHANGE OF THE
!       PION AND RHO
!       
!       MAINLY, THESE ARE FOR THE REACTION e + n --> p + e' + X, AS
!       MIGHT BE MEASURED AT ``BoNuS''
!       
!       WRITTEN: T. Hobbs (APRIL 21, 2014)
!
!       Modified to be used within g4sbs framework
!       Carlos Ayerbe Gayoso (Aug, 2020)
!
!       Parameters to receive:
!       E: energy beam
!       x: x Bjoerken
!       theta_e: scatering electron angle
!       km1: lower limit of scatter proton momentum
!       km2: upper limit of scatter proton momentum
!
!       Note: maybe those limits should modified or elimated, since we want the 
!       structure function for a given x and hadrom momentum
!
!       Modifications: 
!       - removed or commented, all the print statements
!       - removed or commented the file write output
!       - removed or commented the cTEQ stuff (it is calculated in g4sbs)
!
!
!
!       **********************************************************************
!       VARIABLE DECLARATIONS
	IMPLICIT NONE
	INTEGER ix,nx,iy,ikTk,ikTx,FLAG,typ
	PARAMETER (nx=30)
	EXTERNAL fypiN, f_rhoN, f_RhoDel, fypiD, CTQ6PDF
	REAL*8  fypiN, f_rhoN, f_RhoDel, fypiD
        REAL*8  F2pi_ret, F2N_ret !return values
	REAL*8  x,xmin,xmax,xint,y,ymax,y0,yint,L
	REAL*8  theta_e,alpha1,alpha2,E,mN,Q2,pH,EH,z,nu
	REAL*8  xpt(nx),F2pi0,pi,alpha1r, alpha2r, yT
	REAL*8  F2n(nx),F2pi_k(nx),RAT_k(nx)
	REAL*8  km1,km2,ktmink,ktmaxk,ktintk,fpi_intk,frho_intk
	REAL*8  ktk,kmag,pi_kint,rho_kint,F2pi_GRVT,xVpiT,xSpiT
	REAL*8  fpi_x(nx),kTintx,fpi_intx,frho_x(nx),frho_intx
	REAL*8  kTx,kTminx,kTmaxx,kM(nx),cosph
	REAL    CTQ6PDF,u_pro,d_pro,ubar_pro,dbar_pro,F2neu
	REAL*8  Ld,fpiD_intx,frhoD_intx,fpiD_x(nx),frhoD_x(nx)
	REAL*8  HSF0,F2piK(nx) 
    !***********************************************************************
!---------------------------------------------------------------------
!       WE SET THE RENORMALIZATION CUT-OFF PARAMETER LAMBDA "L" = ... IN UNITS OF GeV
!       L = 1.18D0      !COV DIPOLE: HSS NORM
!       L = 1.71D0      !IMF DIPOLE: HSS NORM

!       L = 1.63D0       !HSS +
	L = 1.56D0		!HSS CENT. VALUE
!       L = 1.48D0      !HSS -

!       L = 1.4D0      !WM 1993 PRD
!...    AND USE AN INDEPENDENT PARAMETER FOR THE DELTA:
!       Ld = 0.63D0     !COV DIPOLE: HSS NORM
!       Ld = 1.24D0     !IMF DIPOLE: HSS NORM

!       Ld = 1.46D0     !HSS +
	Ld = 1.39D0		!HSS CENT. VALUE
!       Ld = 1.32D0     !HSS -
!_____________________________________________________________________________
!       KINNI
	yT = 0.05D0		!THE CHOSEN VALUE OF y FOR SF TAGGING

! the bound limits are giving as arguments for the subroutine
! because we want the f2pi value calculated in a certain momentum range (CA)


!       km1 = 0.15D0              !THE LOWER BOUND OF |\vec{k}| !PER THIA KEPPEL!!
!       km2 = km1 + 0.01D0        !THE UPPER BOUND OF |\vec{k}|
!       km2 = km1 + 0.1D0        !THE UPPER BOUND OF |\vec{k}|   A TEST!!!
!       km1 = 0.060D0              !THE LOWER BOUND OF |\vec{k}| !PER THIA KEPPEL!! -- MARCH 26th (CK1)
!       km2 = 0.160D0        !THE UPPER BOUND OF |\vec{k}|

!	km1 = 0.250D0		!THE LOWER BOUND OF |\vec{k}| !PER THIA KEPPEL!! -- MARCH 26th (CK2!)
!	km2 = 0.400D0		!THE UPPER BOUND OF |\vec{k}|



!       km1 = 0.0D0              !THE LOWER BOUND OF |\vec{k}| !NO CUTS!
!       km2 = 10.D0        !THE UPPER BOUND OF |\vec{k}|
!_____  (WANT: km1 = 0.060, 0.080, 0.100, 0.130, 0.160)_________________________
!***********************************************************************
!       HERE WE PLACE A GLOBAL FLAG FOR THE CHOICE OF THE WAVEFUNCTION SUPPRESSION FACTOR
!       typ = 1 !DIPOLE FORM FACTOR
	typ = 2			!EXPONENTIAL FORM FACTOR
!       typ = 3 !COV. DIPOLE FORM FACTOR
!***********************************************************************
!       THIS CHOOSES AMONG THE VARIOUS POSSIBLE COMBINATIONS OF SPLITTING FUNCTION/PDF
	FLAG = 0
!       THE FLAGS TOGGLE AMONG DISSOCIATION MODES AS: 
!       FLAG = 0  --- THE PION CONTRIBUTION  | J = 0 + 1/2
!       FLAG = 1  --- THE RHO CONTRIBUTION   | J = 1 + 1/2
!---------------------------------------------------------------------
!***********************************************************************
!       EXPERIMENTAL INPUTS AS GIVEN BY THIA KEPPEL (11.12.2013) ------------
	pi = 4.D0*DATAN(1.D0)
!	E = 11.D0               !ELECTRON BEAM ENERGY, GeV
!	theta_e = 35.D0         !ELECTRON SCATTERING ANGLE, DEGREES
	pH = 0.325D0            !PROD. HADRON MOMENTUM, GeV: [.05 - .50]
	alpha1 = 30.D0          !PROD. HADRON ANGLE LOWER BOUND, DEGREES
	alpha2 = 70.D0          !PROD. HADRON ANGLE UPPER BOUND, DEGREES
	alpha1r = (alpha1/180.D0)*pi
	alpha2r = (alpha2/180.D0)*pi
!---------------------------------------------------------------------
!***********************************************************************

!       print *, 'IN f2pi_sub.f90'
!       print *, E, x, theta_e, km1, km2

!       WE NOW COMPUTE CONVOLUTIONS AND DISTRIBUTIONS OVER THE RELEVANT x
	ix = 1
	xmin = 0.05D0
	xmax = 0.33D0 ! normal was 0.16
	xint = (xmax-xmin)/nx
!***************************************
	F2pi0 = 0.D0
	HSF0 = 0.D0
!	DO   x=xmin,xmax+xint/10,xint ! no x loop (CA)

	   xpt(ix) = x
!        print *, 'x = ', xpt(ix) 
!---------------------------------------------------------------------
!       CALCULATE USEFUL KINEMATICAL VARIABLES:
	   mN = 0.93891897D0	!TARGET NEUTRON MASS, GeV

	   Q2 = 2.D0*mN*x*E &	![GeV]**2
        * (1.D0 - 1.D0/( (2.D0*E/(x*mN)) &
        * DSIN((theta_e/180.D0)*pi)**2.D0 + 1.D0 ) )

!          print *, 'Q2: ', Q2

	   nu = Q2 / (2.D0 * mN * x) !DIS MOMENTUM EXCHANGE, GeV

	   EH = DSQRT(pH**2.D0 + mN**2.D0)

	   z = EH / nu		!SIDIS 'z' FRACTION

!---------------------------------------------------------------------
!       DEFINE THE BOUNDS OF THE CONVOLUTION INTEGRAL -- x to 1
	   F2piK(ix) = 0.D0
	   iy = 0
	   y0 = x
       ymax = 1.D0      !'ymax' FIXED BY |k| BOUNDS
	   ymax = 0.329D0	!SHARP y CUTOFF (was 0.15)
	   yint = (ymax-y0)/100.D0
	   DO y=y0,ymax,yint
!_____________________________________________________________________________
!_____________________________________________________________________________
!------ THE TAGGED SF FOR FIXED BINS OF |\vec{k}| AT FIXED kT ----------------
!       INTEGRATION RANGES DEPENDENT UPON THE |\vec{k}| MAGNITUDE PER SF TAGGING
	      ikTk = 0
	      kTmink = 0.D0	!FOR FULLY INTEGRATED COMPARISON
	      kTmaxk = 1.D0
	      kTintk = (kTmaxk - kTmink)/50000.D0

	      fpi_intk = 0.D0
	      frho_intk = 0.D0

	      DO kTk = kTmink, kTmaxk, kTintk
!...... ANGLE CONSTRAINT  (45 deg: cosph < 0.707; 30 deg: cosph < 0.866)
		 kmag = DSQRT( kTk**2.D0 + (1.D0/(4.D0*mN**2*(1.D0-y)**2.D0))&
     	      * (kTk**2.D0 + (1.D0 - (1.D0-y)**2.D0)*mN**2.D0)**2.D0 )

		 cosph =  (1.D0/(2.D0*mN*(1.D0-y))) &
                     * (kTk**2.D0 + (1.D0 - (1.D0-y)**2.D0)*mN**2.D0) &
                     / kmag

!       IF (cosph.GT.0.707D0) THEN  !FOR 45 DEGREE CUT
		 IF (cosph.GT.0.866D0) THEN !FOR 30 DEGREE CUT
		    pi_kint = 0.D0
		    rho_kint = 0.D0
		 ELSE
!		    IF (kmag.LT.km1) THEN !COMMENT OUT FOR FULLY INCLUSIVE COMP.!!
!		       pi_kint = 0.D0
!		       rho_kint = 0.D0
!		    ELSEIF (kmag.GT.km2) THEN
!		       pi_kint = 0.D0
!		       rho_kint = 0.D0
!		    ELSE
		       pi_kint = 2.D0*kTk*fypiN(y,kTk,L,typ,0) ! factor two is due to isospin
		       rho_kint = 2.D0*kTk*f_rhoN(y,kTk,L,typ,0) ! if we want proton (pi0), output divided by 2
!		    ENDIF
		 ENDIF

		 IF (ikTk.EQ.0) THEN
		    fpi_intk = fpi_intk + pi_kint
		    frho_intk = frho_intk + rho_kint
		 ELSE IF (ikTk/2*2.NE.ikTk) THEN
		    fpi_intk = fpi_intk + 4.D0*pi_kint
		    frho_intk = frho_intk + 4.D0*rho_kint
		 ELSE IF (ikTk/2*2.EQ.ikTk) THEN
		    fpi_intk = fpi_intk + 2.D0*pi_kint
		    frho_intk = frho_intk + 2.D0*rho_kint
		 ENDIF
		 
		 ikTk = ikTk + 1
	      ENDDO ! end DO kTk
	      fpi_intk = (kTintk/3.D0)*fpi_intk
	      frho_intk = (kTintk/3.D0)*frho_intk
             
          
	      CALL GRV (x/y,Q2,xVpiT,xSpiT)

	      F2pi_GRVT = (5.D0/9.D0) * (xVpiT + 2.D0 * xSpiT)

       

!       ------------ INTEGRATED, OVER FINITE y ------------------------
	      IF (FLAG.EQ.0) THEN
		 F2pi_k(ix) = fpi_intk * F2pi_GRVT
!PION CONTRIBUTION
	      ELSE IF (FLAG.EQ.1) THEN
		 F2pi_k(ix) = frho_intk * F2pi_GRVT
!RHO CONTRIBUTION
	      ENDIF
!________________________________________________________________________________________
!________________________________________________________________________________________
!       WE ADD THE EVALUATED INTEGRAND TO THE CUMULANT IN A STANDARD RUNGE-KUTTA METHOD
!-----  FOR THE x DISTRIBUTIONS
	      IF (iy.EQ.0) THEN
		 F2piK(ix) = F2piK(ix) + F2pi_k(ix)

	      ELSE IF (iy/2*2.NE.iy) THEN
		 F2piK(ix) = F2piK(ix) + 4.D0*F2pi_k(ix)
	      ELSE IF (iy/2*2.EQ.iy) THEN
		 F2piK(ix) = F2piK(ix) + 2.D0*F2pi_k(ix)
	      ENDIF

	      iy = iy + 1
	   ENDDO !end DO y
	   F2piK(ix) = (yint/3.D0) * F2piK(ix)

           F2pi_ret = F2piK(ix)! the pion SF to return (CA)
!	   print*, x, km1, F2piK(ix),   F2piK99(ix) 
!	   print*, x, km1, F2piK(ix),   F2pi_ret 

!_____________________________________________________________________________
!_____________________________________________________________________________
!------ THE NEUTRON SF F_2 AS GENERATED FROM CTEQ PDFs ASSUMING CS -----------
!       N.B.: MUST BE TYPED AS A REAL FUNCTION, ETC.!
!	   CALL SETCTQ6(1)	! CTEQ 'MS-bar' SCHEME. 
!	   u_pro = CTQ6PDF (1, REAL(x), SQRT(REAL(Q2)))
!	   ubar_pro = CTQ6PDF (-1, REAL(x), SQRT(REAL(Q2)))
!	   d_pro = CTQ6PDF (2, REAL(x), SQRT(REAL(Q2)))
!	   dbar_pro = CTQ6PDF (-2, REAL(x), SQRT(REAL(Q2)))

!	   F2neu = 2.*x * ((4./9.)*(d_pro + dbar_pro) &
!	   + (1./9.)*(u_pro + ubar_pro))
!	   F2n(ix) = DBLE(F2neu)




!-------THE FIXED |\vec{k}| RATIO WRT F2n --------------
!	   RAT_k(ix) = F2piK(ix) / F2n(ix)
!       kM(ix) = (mN/2.D0) * x * ( (2.D0-x)/(1.D0-x) ) !|k| for kT=0
	   kM(ix) = mN * x	!APPROX.
!---------------------------------------------------------
!
!  I don't think we need this calculations for our generator
!  perhaps removing them could faster the program
!
!____   WE WANT TO PLOT f_{pi-N}(x)  ______________________
	   ikTx = 0
	   kTminx = 0.D0
	   kTmaxx = 10.D0
	   kTintx = (kTmaxx - kTminx)/1000.D0

	   fpi_intx = 0.D0
	   fpiD_intx = 0.D0
	   frho_intx = 0.D0
	   frhoD_intx = 0.D0


	   DO kTx = kTminx, kTmaxx, kTintx
	      IF (ikTx.EQ.0) THEN
		 fpi_intx = fpi_intx + 2.D0*kTx*fypiN(x,kTx,L,typ,0)
		 fpiD_intx = fpiD_intx + 2.D0*kTx*fypiD(x,kTx,Ld,typ,0)
		 frho_intx = frho_intx + 2.D0*kTx*f_rhoN(x,kTx,L,typ,0)
		 frhoD_intx = frhoD_intx + 2.D0*kTx*f_RhoDel(x,kTx,Ld,typ,0)
	      ELSE IF (ikTx/2*2.NE.ikTx) THEN
		 fpi_intx = fpi_intx + 8.D0*kTx*fypiN(x,kTx,L,typ,0)
		 fpiD_intx = fpiD_intx + 8.D0*kTx*fypiD(x,kTx,Ld,typ,0)
		 frho_intx = frho_intx + 8.D0*kTx*f_rhoN(x,kTx,L,typ,0)
		 frhoD_intx = frhoD_intx + 8.D0*kTx*f_RhoDel(x,kTx,Ld,typ,0)
	      ELSE IF (ikTx/2*2.EQ.ikTx) THEN
		 fpi_intx = fpi_intx + 4.D0*kTx*fypiN(x,kTx,L,typ,0)
		 fpiD_intx = fpiD_intx + 4.D0*kTx*fypiD(x,kTx,Ld,typ,0)
		 frho_intx = frho_intx + 4.D0*kTx*f_rhoN(x,kTx,L,typ,0)
		 frhoD_intx = frhoD_intx + 4.D0*kTx*f_RhoDel(x,kTx,Ld,typ,0)
	      ENDIF

	      ikTx = ikTx + 1
	   ENDDO !end DO kTx

	   fpi_x(ix) = (kTintx/3.D0)*fpi_intx
	   fpiD_x(ix) = (kTintx/3.D0)*fpiD_intx
	   frho_x(ix) = (kTintx/3.D0)*frho_intx
	   frhoD_x(ix) = (kTintx/3.D0)*frhoD_intx

!__________________________________________________________

	   IF (ix/2*2.NE.ix) THEN
	      F2pi0 = F2pi0 + 4.D0*F2piK(ix)
	      HSF0 = HSF0 + 4.D0 *fpiD_x(ix)
	   ELSE IF (ix/2*2.EQ.ix) THEN
	      F2pi0 = F2pi0 + 2.D0*F2piK(ix)
	      HSF0 = HSF0 + 2.D0 *fpiD_x(ix)
	   ENDIF
!	   ix = ix + 1 ! no x loop (CA)
!	ENDDO !end DO x
	F2pi0 = (xint/3.D0) * F2pi0
	HSF0 = (xint/3.D0) * HSF0

!	print*, 'First moment of F2piK(x)=',F2pi0 ! NO NEED TO PRINT (CA)
!	print*, 'First moment of f(y)=',HSF0 ! NO NEED TO PRINT (CA)
!________________________________________________________________________________________
!________________________________________________________________________________________
!...    WRITE DATA TO FILE

!_____  (WANT: km1 = 0.060, 0.080, 0.100, 0.130, 0.160)_________________________
!       DASSI
!	IF (FLAG.EQ.0) THEN
!	   OPEN (14,FILE='BONUS-AP14/F2-piN_kL=CK2_Ang-cuts.dat', & !F2piK(ix)
!	STATUS='UNKNOWN', FORM='FORMATTED')
!	ELSE IF (FLAG.EQ.1) THEN
!	   OPEN (14,FILE='BONUS-AP14/F2-rhoN_kL=CK1.dat', &
!	STATUS='UNKNOWN', FORM='FORMATTED')
!	ENDIF
!	DO ix=2,nx ! no x loop (CA)
! 	   WRITE (14,*) xpt(ix), F2piK(ix)
!	ENDDO ! no x loop (CA)
!	CLOSE (14)

!________________________________________________________________________________________
!	IF (FLAG.EQ.0) THEN
!	   OPEN (15,FILE='BONUS-AP14/RF2-piN_kL=CK1_yint.dat', & !F2piK(ix) / F2n(ix)
!	STATUS='UNKNOWN', FORM='FORMATTED')
!	ELSE IF (FLAG.EQ.1) THEN
!	   OPEN (15,FILE='BONUS-AP14/RF2-rhoN_kL=CK1.dat', &
!	STATUS='UNKNOWN', FORM='FORMATTED')
!	ENDIF
!	DO ix=2,nx ! no x loop (CA)
!	   WRITE (15,*) xpt(ix), RAT_k(ix)
!	ENDDO ! no x loop (CA)
!	CLOSE (15)

!	PRINT*, 'THE NEW DATA HAVE BEEN WRITTEN!'

	END
!       ***********************************************************************
!       ***********************************************************************
	FUNCTION fypiN (y,kT,L,typ,dis)
!       
!       FUNCTION GIVING f(y) FOR N-PION VERTEX, WHERE
!       y IS THE IMF MOMENTUM OF THE INTERMEDIATE MESON.
!       
!       WRITTEN: W. Melnitchouk (1999)
!       MODIFIED: T. HOBBS (2013)
!       ***********************************************************************
	IMPLICIT NONE
	INTEGER typ,dis
	REAL*8  ss,kT,kT2,SpiN
	REAL*8  y,fypiN,t
	REAL*8  pi,mN,mpi,mP,g_piNN,gg,FF,L,sM

	pi = 4*DATAN(1.D0)
	mN  = 0.93891897D0	!masses in GeV!!
	IF (dis.EQ.0) THEN
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!       !**** THE DISSOCIATION P --> Dbar0 + LAMBDA_c^+ 
	   mpi = 0.13957018D0	!THE MASS OF THE pi^-
	   mP = mN		!THE MASS OF THE PROTON
!       !*** WE USE THE COUPLINGS INFERRED FROM HAIDENBAUER ET AL. ***
	   g_piNN = DSQRT (14.40D0 * 4*pi) ! g_{pi NN}
	   gg = 2.D0 * g_piNN**2 / (16.D0 * pi**2) !2 ISOSPIN FACTOR
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	ELSE IF (dis.EQ.1) THEN
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!       !**** THE DISSOCIATION P --> Dbar0 + SIGMA_c^+ 
	   mpi = 1.865D0	!THE MASS OF THE Dbar0
	   mP = 2.4529D0	!THE MASS OF THE CHARMED SIGMA
!       !*** WE USE THE COUPLINGS INFERRED FROM HAIDENBAUER ET AL. ***
	   g_piNN = DSQRT (0.576D0 * 4*pi) ! as N-D-Sigma_c
	   gg = g_piNN**2 / (16.D0 * pi**2) !1 ISOSPIN FACTOR
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	ENDIF
	IF (y.LE.0.D0 .OR. y.GE.1.D0) THEN
	   fypiN = 0.D0
	   RETURN
	ENDIF

	kT2 = kT**2
	SpiN = (kT2 + mpi**2)/y + (kT2 + mP**2)/(1.D0-y)

	IF (typ.EQ.0) THEN
	   FF = ((L**2 + mN**2) / (L**2 + SpiN)) ! monopole
	ELSE IF (typ.EQ.1) THEN
	   FF = ((L**2 + mN**2) / (L**2 + SpiN))**2 ! dipole
	ELSE IF (typ.EQ.2) THEN
	   FF = DEXP( (mN**2 - SpiN)/(L**2) ) ! expon
	ELSE IF (typ.EQ.3) THEN
	   t = (- kT2 - mN**2*y**2) /(1.D0-y)
	   FF = ((L**2 - mpi**2) / (L**2 - t))**2 ! cov dip
	ELSE IF (typ.EQ.4) THEN	! DIPOLE -- s-channel Lambda exchange
	   sM = (kT2 + (1.D0+y)*mpi**2.D0)/y + (kT2 &
	+ y*mP**2.D0)/(1.D0-y) + mN**2.D0
	   FF = (L**4 + mP**4)/(L**4 + sM**2) 
	ENDIF

        ss = ( kT2 + (mP - (1.D0-y) * mN)**2.D0) / (1.D0-y) &
             / ( (1.D0-y)*(SpiN - mN**2.D0) )**2.D0 * FF**2.D0

!       ss = ( kT2 + (mP - (1.D0-y) * mN)**2.D0) / (1.D0-y)  !NO FORM FACTOR!
!       &       / ( (1.D0-y)*(SpiN - mN**2.D0) )**2.D0

        fypiN =  gg * (1.D0-y) / y * ss
	RETURN   
	END
!       ***************************************************************************
        FUNCTION f_rhoN (y,kT,L,typ,dis)
!       Function giving numerical value of f(y) for N-rho-N
!       OUTPUT IS THE NEUTRON ---> rho^- + PROTON SPLITTING FUNCTION
!       By: T. Hobbs on NOV 12, 2013
!       taken from notes "spin-1, m_B /= m_N," May 17, 2012
!       ***************************************************************************
        IMPLICIT NONE
        INTEGER typ,dis
        REAL*8  kT,kT2,SRoN,P_k,pl_k,P_p,t,sv,st,si,ss
        REAL*8  y,f_rhoN,sM
        REAL*8  pi,mN,mrho,mP,g_RoNN,f_RoNN,gg,fff,fg,FF,L
	
	pi = 4*DATAN(1.D0)
	mN  = 0.93891897D0	!masses in GeV!!
	IF (dis.EQ.0) THEN
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!****   THE DISSOCIATION N --> rho^- + PROTON
	   mrho = 0.7754D0	!THE MASS OF THE rho^-
	   mP = mN		!THE MASS OF THE PROTON
!***    WE USE THE COUPLINGS FROM SU(2) SYMMETRY FOR rho-N-N***
	   g_RoNN = DSQRT (2.D0 * 0.55D0 * 4.D0 * pi) !rho-N-N (Hohler and Pieteranin)
	   f_RoNN = 6.1D0 * g_RoNN
	   gg = g_RoNN**2 / (16.D0 * pi**2)
	   fff = f_RoNN**2 / (16.D0 * pi**2)
	   fg = f_RoNN*g_RoNN / (16.D0 * pi**2) !FACTOR OF 2: ISOSPIN
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	ELSEIF (dis.EQ.1) THEN
!****   THE DISSOCIATION P --> Dbar*0 + LAMBDA_c^+ 
	   mrho = 0.7754D0	!THE MASS OF THE rho^-
	   mP = mN		!THE MASS OF THE PROTON
!***    WE USE THE COUPLINGS FROM SU(2) SYMMETRY FOR rho-N-N***
	   g_RoNN = DSQRT (0.55D0 * 4.D0 * pi) !rho-N-N (Hohler and Pieteranin)
	   f_RoNN = 6.1D0 * g_RoNN
	   gg = g_RoNN**2 / (16.D0 * pi**2)
	   fff = f_RoNN**2 / (16.D0 * pi**2)
	   fg = f_RoNN*g_RoNN / (16.D0 * pi**2)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	ENDIF

        IF (y.LE.0.D0 .OR. y.GE.1.D0) THEN
	   f_rhoN = 0.D0
	   RETURN
        ENDIF

	kT2 = kT**2
	SRoN = (kT2 + mrho**2)/y + (kT2 + mP**2)/(1.D0-y)
	
	IF (typ.EQ.0) THEN
	   FF = ((L**2 + mN**2) / (L**2 + sRoN)) ! MONOPOLE
	ELSE IF (typ.EQ.1) THEN
	   FF = ((L**2 + mN**2) / (L**2 + sRoN))**2 ! DIPOLE
	ELSE IF (typ.EQ.2) THEN
	   FF = DEXP( (mN**2 - SRoN)/(L**2) ) ! EXPONENTIAL
	ELSE IF (typ.EQ.3) THEN
	   t = (- kT2 - mN**2*y**2) /(1.D0-y)
	   FF = ((L**2 - mrho**2) / (L**2 - t))**2 ! COV. DIPOLE
	ELSE IF (typ.EQ.4) THEN	! DIPOLE -- s-CHANNEL LAMBDA EXCHANGE
	   sM = (kT2 + (1.D0+y)*mrho**2.D0)/y + (kT2 &
	+ y*mP**2.D0)/(1.D0-y) + mN**2.D0
	   FF = (L**4 + mP**4)/(L**4 + sM**2) 
	ENDIF

!C...    TOPT WITH P[alpha] - p[alpha] DERIVATIVE COUPLING
	P_k = (mrho**2 + y**2*mN**2 + kT2)/2.D0/y
	P_p = (mP**2 + (1.D0-y)**2*mN**2 + kT2)/2.D0/(1.D0-y)
	pl_k = (mP**2+kT2)*y/2.D0/(1.D0-y) &
      + (mrho**2+kT2)*(1.D0-y)/2.D0/y + kT2  

	sv = -6.D0*mN*mP + 4.D0*P_k*pl_k/mrho**2 + 2.D0*P_p

	st = -(P_p)**2 + P_p*(mP+mN)**2 - mP*mN*(mP**2+mN**2+mP*mN) &
      + 1.D0/(2.D0*mrho**2) * ( (P_p - mP*mN)*(P_k-pl_k)**2 &
      - 2.D0*(P_k-pl_k)*(mP**2*P_k - mN**2*pl_k) &
      + 2.D0*P_k*pl_k*(2.D0*P_p-mP**2-mN**2) ) 

	si = -4.D0*(mP+mN)*(mP*mN - P_p) &
      -2.D0*(mP*P_k**2 - (mP+mN)*P_k*pl_k + mN*pl_k**2)/mrho**2
	
	IF (kt.gt.0.0 .and. (sv.lt.0.0 .or. st.lt.0.0)) THEN
	   PRINT *,'CS1 -- ##### kT,y,sv,st =',kt,y,sv,st
	   STOP 
	ENDIF
        
	ss = (gg*sv + fff*st/mN**2 + fg*si/mN) & !THE FULL
      / ( y*(SRoN - mN**2) )**2 * FF**2 !EXPRESSION (UN-INT.)**
!__________________________________________________________________________
!************TO STUDY SEPARATE TERMS' CONTRIBUTIONS!!!! ******************
!       ss = gg*sv                                         
!       &       / ( y*(SRoN - mN**2) )**2 * FF**2  * (2.D0*kT)  !VECTOR
!       ss = fg*si/mN                                      
!       &       / ( y*(SRoN - mN**2) )**2 * FF**2  * (2.D0*kT)  !INTERFERENCE
!       ss = fff*st/mN**2                                  
!       &       / ( y*(SRoN - mN**2) )**2 * FF**2  * (2.D0*kT)  !TENSOR
!__________________________________________________________________________
        f_rhoN = y / (1.D0-y) * ss
        RETURN
        END
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!       ***************************************************************************
!       ***************************************************************************
        FUNCTION f_RhoDel (y,kT,L,typ,dis)
!       Function giving numerical value of f(y) for N-Rho-Delta
!       OUTPUT IS THE NEUTRON ---> RHO + DELTA SPLITTING FUNCTION
!       By: T. Hobbs on JAN 7, 2014
!       ***************************************************************************
        IMPLICIT NONE
        INTEGER typ,dis
        REAL*8  kT,kT2,SRoD,P_k,pl_k,P_p,t
        REAL*8  y,f_RhoDel,sr,ss
        REAL*8  pi,mN,mD,mRo,g_NDS,gg,FF,L,sM
	
	pi = 4*DATAN(1.D0)
	mN  = 0.93891897D0	!masses in GeV!!
	IF (dis.EQ.0) THEN
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!       !**** THE DISSOCIATION N --> (RHO^0 + DELTA^0) + (RHO^- + DELTA^+)
	   mRo = 0.7754D0	!THE MASS OF THE RHO MESON
	   mD = 1.232D0		!Mass OF THE DELTA ISOBAR
!       !*** WE USE THE COUPLINGS INFERRED FROM HOLZENKAMP ET AL. ***
	   g_NDS = DSQRT (20.448D0 * 4.D0*pi) ! as N-D*-Sigma*_c
	   gg = (2.D0/3.D0) * g_NDS**2 / (16.D0 * pi**2 * mRo**2)
!       WE INCLUDE AN OVERALL FACTOR OF 2/3 TO ACCOUNT FOR ISOSPIN!!!!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	ELSE IF (dis.EQ.1) THEN
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!       !**** A BLANK SPACE FOR OTHER SEPARATE MODES
	   mRo = 0.7754D0	!THE MASS OF THE RHO MESON
	   mD = 1.232D0		!Mass OF THE DELTA BARYON
!       !*** TERMINATE THE OUTPUT ***
	   g_NDS = 0.D0		! ZERO OUTPUT
	   gg = g_NDS**2 / (16.D0 * pi**2 * mRO**2)
!       WE INCLUDE AN OVERALL FACTOR OF 1 TO ACCOUNT FOR ISOSPIN!!!!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	ENDIF
        IF (y.LE.0.D0 .OR. y.GE.0.999D0) THEN
	   f_RhoDel = 0.D0
	   RETURN
        ENDIF
	
	kT2 = kT**2
	SRoD = (kT2 + mRo**2)/y + (kT2 + mD**2)/(1.D0-y)
        
	IF (typ.EQ.0) THEN
	   FF = ((L**2 + mN**2) / (L**2 + sRoD)) ! monopole
	ELSE IF (typ.EQ.1) THEN
	   FF = ((L**2 + mN**2) / (L**2 + sRoD))**2 ! dipole
	ELSE IF (typ.EQ.2) THEN
	   FF = DEXP( (mN**2 - SRoD)/(L**2) ) ! expon
	ELSE IF (typ.EQ.3) THEN
	   t = ( -kT2 - (1.D0-y)*(mRo**2 - y*mN**2) ) / y
	   FF = ((L**2 - mD**2) / (L**2 - t))**2 ! t-channel dip
	ELSE IF (typ.EQ.4) THEN
	   sM = (kT2 + (2.D0-y)*mRo**2.D0)/(1.D0-y) + (kT2 &
         + (1.D0-y)*mD**2.D0)/y + mN**2.D0
	   FF = (L**4 + mD**4)/(L**4 + sM**2) !DIPOLE -- s-channel Lambda exchange
	ENDIF

!...    TOPT with P[alpha] - p[alpha] derivative coupling
	P_k = (mRo**2 + y**2*mN**2 + kT2)/2.D0/y
	P_p = (mD**2 + (1.D0-y)**2*mN**2 + kT2)/2.D0/(1.D0-y)
	pl_k = (mD**2+kT2)*y/2.D0/(1.D0-y) &
      + (mRo**2+kT2)*(1.D0-y)/2.D0/y + kT2  

	sr = -4.D0*mN*mD/3.D0*(2.D0*mD**2.D0+mN*mD+2.D0*mN**2.D0) &
      -4.D0*mN*mD/(3.D0*mRo**2.D0)*(P_k-pl_k)**2.D0 &
      -4.D0/(3.D0*mRo**2.D0)*(mD**2.D0*P_k**2.D0+mN**2.D0*pl_k**2.D0) &
      +4.D0*P_p/3.D0*(2.D0*mD**2.D0+4.D0*mN*mD+mN**2.D0) &
      +4.D0*P_p/(3.D0*mRo**2.D0)*pl_k**2.D0*(1.D0-mN**2.D0/mD**2.D0) &
      -4.D0*P_p**2.D0*(1.D0-2.D0*P_k*pl_k/(3.D0*mRo**2.D0*mD**2.D0) &
      -P_p/(3.D0*mD**2.D0))
	
	if (kt.gt.0.0 .and. (sr.lt.0.0)) then
	   print *,'CS2 -- ##### kT,y,sr =',kt,y,sr
	   stop
	endif
        
	ss = sr &
      /((1.D0-y)*(SRoD - mN**2))**2 * FF**2

        f_RhoDel = gg * (1.D0-y) / y * ss

        RETURN
        END
!       ***************************************************************************
!       ***********************************************************************
	FUNCTION fypiD (y,kT,L,typ,dis)
!       
!       Function giving numerical value of f(y) for the SU(4) analogue of the pi-Delta
!       interaction, as usual y is IMF momentum fraction of the charmed BARYON.
!       ***********************************************************************
	IMPLICIT NONE
	INTEGER typ,dis
	REAL*8  ss,kT,kT2,SpiD
	REAL*8  y,yR,fypiD,t,sM
	REAL*8  pi,mN,mD,mpi,g_DLcN,gg,FF,L

	pi = 4*DATAN(1.D0)
	mN  = 0.93891897D0	!masses in GeV!!
	IF (dis.EQ.0) THEN
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!       !**** THE DISSOCIATION N ---> (pi^0 + DELTA^0) + (pi^- + DELTA^-)
	   mD = 1.232D0		!THE MASS OF THE DELTA BARYON
	   mpi = 0.13957018D0	!Mass OF THE PION 
!       !*** WE USE THE COUPLINGS INFERRED FROM HAIDENBAUER ET AL. ***
	   g_DLcN = DSQRT (0.2237D0 * 4*pi) ! as N-pi-Del from Holzenkamp et al.
	   gg =  (2.D0/3.D0) * g_DLcN**2 / (16.D0 * pi**2) !2/3 ISOSPIN FACTOR
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	ELSE IF (dis.EQ.1) THEN
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!       !**** THE NULL DISSOCIATION
	   mD = 1.232D0		!THE MASS OF THE DELTA
	   mpi = 0.13957018D0	!Mass OF THE PION
!       !*** WE USE ZERO COUPLINGS TO RETURN A NULL OUTPUT
	   g_DLcN = 0.D0    
	   gg = g_DLcN**2 / (16.D0 * pi**2)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	ENDIF

	IF (y.LE.0.D0 .OR. y.GE.0.999D0) THEN
	   fypiD = 0.D0
	   RETURN
	ENDIF

	kT2 = kT**2
	SpiD = (kT2 + mpi**2)/y + (kT2 + mD**2)/(1.D0-y)
!********************************************************************************
	IF (typ.EQ.0) THEN
	   FF = ((L**2 + mN**2) / (L**2 + SpiD)) ! monopole
	ELSE IF (typ.EQ.1) THEN
	   FF = ((L**2 + mN**2) / (L**2 + SpiD))**2 ! dipole
	ELSE IF (typ.EQ.2) THEN
	   FF = DEXP( (mN**2 - SpiD)/(L**2) ) ! expon
	ELSE IF (typ.EQ.3) THEN
	   t = (- kT2 - mN**2*y**2) /(1.D0-y)
	   FF = ((L**2 - mpi**2) / (L**2 - t))**2 ! cov dip
	ELSE IF (typ.EQ.4) THEN	! DIPOLE -- s-channel Lambda exchange
	   sM = (kT2 + (1.D0+y)*mpi**2.D0)/y + (kT2 &
	+ y*mD**2.D0)/(1.D0-y) + mN**2.D0
	   FF = (L**4 + mD**4)/(L**4 + sM**2) 
	ENDIF

	yR = 1.D0 - y
	ss = ( kT2 + (mD-yR*mN)**2) * ( kT2 + (mD+yR*mN)**2 )**2 &
      / ( 6.D0*mD**2*yR**3) / ( (1.D0-yR)*(SpiD - mN**2) )**2 &
     * FF**2

	fypiD =  gg/(mpi**2) * (1.D0-yR) / yR * ss

	RETURN   
	END
!       ***************************************************************************
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!       THE PION STRUCTURE FUNCTION IS TAKEN FROM THIS SUBROUTINE
!       ***************************************************************************
        SUBROUTINE GRV (x,Q2,xVpi,xSpi)
!       
!       Subroutine giving x-dependence of parametrization of LO pion valence
!       and sea distribution functions xVpi, xSpi, valid for 0.25 < Q2 < 10^8 GeV^2,
!       and 10^-5 < x < 1.
!       
!       Gluck, Reya, Vogt: Z.Phys. C53 (1992) 651 (appendix 1).
!       ******************************************************************************
!       IMPLICIT UNDEFINED (A-Z)
!       IMPLICIT NONE (A-Z)
        REAL*8  x,xVpi,a,AA,D,Nv, &
             xSpi,alpha,as,AAs,Bs,Ds,E,Epr,beta
        REAL*8  Q2,Q02,L,s

        xVpi = 0.D0
        xSpi = 0.D0

	IF (x.LT.1.D-5) x=1.01D-5
	IF (Q2.LT.0.25D0) Q2=0.2501D0

        Q02 = 0.25D0
        L = 0.232D0
        IF (Q2.LE.Q02) RETURN
        s = DLOG( DLOG(Q2/L**2) / DLOG(Q02/L**2) )

!...    Valence distribution
        Nv = 0.519D0 + 0.180D0*s - 0.011D0*s**2
        a = 0.499D0 - 0.027D0*s
        AA = 0.381D0 - 0.419D0*s
        D = 0.367D0 + 0.563D0*s
        xVpi = 0.D0
        IF (x.LT.1.D0) &
             xVpi = Nv * x**a * (1.D0+AA*DSQRT(x)) * (1.D0-x)**D

!...    Sea distribution (SU(3) symmetric)
        alpha = 0.55D0
        as = 2.538D0 - 0.763D0*s
        AAs = -0.748D0
        Bs = 0.313D0 + 0.935D0*s
        Ds = 3.359D0
        E = 4.433D0 + 1.301D0*s
        Epr = 9.30D0 - 0.887D0*s
        beta = 0.56D0
        xSpi = 0.D0
        IF (x.LT.1.D0) &
             xSpi = s**alpha / (DLOG(1.D0/x))**as &
             * (1.D0 + AAs*DSQRT(x) + Bs*x) * (1.D0 - x)**Ds &
             * DEXP(-E + DSQRT(Epr*s**beta*DLOG(1.D0/x)))

  !      IF (xVpi.EQ.0) THEN
  !         print*, 'GRV (sub):', x, Q2, xVpi, xSpi
  !      ENDIF


        RETURN
        END
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

