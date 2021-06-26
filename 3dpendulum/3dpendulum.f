C**   The name of this program is 3dpendulum.f
C**   Created by AUTOLEV 3.2 on Fri Jun 11 17:49:27 2021

      IMPLICIT         DOUBLE PRECISION (A - Z)
      INTEGER          ILOOP, IPRINT, PRINTINT
      CHARACTER        MESSAGE(99)
      EXTERNAL         EQNS1
      DIMENSION        VAR(9)
      COMMON/CONSTNTS/ G,IX,IY,IZ,LA,LAO,M
      COMMON/VARIBLES/ Q1,Q2,Q3,Q4,Q5,Q6,U4,U5,U6
      COMMON/ALGBRAIC/ RX,RY,RZ,U1,U2,U3,Q1p,Q2p,Q3p,Q4p,Q5p,Q6p,U4p,U5p
     &,U6p,KE,P1X,P1Y,P1Z,P2X,P2Y,P2Z,PE,TE
      COMMON/MISCLLNS/ PI,DEGtoRAD,RADtoDEG,COEF(3,3),RHS(3)

C**   Open input and output files
      OPEN(UNIT=20, FILE='3dpendulum.in', STATUS='OLD')
      OPEN(UNIT=21, FILE='3dpendulum.1',  STATUS='UNKNOWN')
      OPEN(UNIT=22, FILE='3dpendulum.2',  STATUS='UNKNOWN')

C**   Read message from input file
      READ(20,7000,END=7100,ERR=7101) (MESSAGE(ILOOP),ILOOP = 1,99)

C**   Read values of constants from input file
      READ(20,7010,END=7100,ERR=7101) G,IX,IY,IZ,LA,LAO,M

C**   Read the initial value of each variable from input file
      READ(20,7010,END=7100,ERR=7101) Q1,Q2,Q3,Q4,Q5,Q6,U4,U5,U6

C**   Read integration parameters from input file
      READ(20,7011,END=7100,ERR=7101) TINITIAL,TFINAL,INTEGSTP,PRINTINT,
     &ABSERR,RELERR

C**   Write heading(s) to output file(s)
      WRITE(*, 6021) (MESSAGE(ILOOP), ILOOP = 1,99) 
      WRITE(21,6021) (MESSAGE(ILOOP), ILOOP = 1,99) 
      WRITE(22,6022) (MESSAGE(ILOOP), ILOOP = 1,99) 

C**   Degree to radian conversion
      PI       = 4*ATAN(1.0D0)
      DEGtoRAD = PI/180.0D0
      RADtoDEG = 180.0D0/PI

C**   Evaluate constants
      U1 = 0
      Q1p = U1
      U2 = 0
      Q2p = U2
      U3 = 0
      Q3p = U3

C**   Initialize time, print counter, variables array for integrator
      T      = TINITIAL
      IPRINT = 0
      VAR(1) = Q1
      VAR(2) = Q2
      VAR(3) = Q3
      VAR(4) = Q4
      VAR(5) = Q5
      VAR(6) = Q6
      VAR(7) = U4
      VAR(8) = U5
      VAR(9) = U6

C**   Initalize numerical integrator with call to EQNS1 at T=TINITIAL
      CALL KUTTA(EQNS1, 9, VAR, T, INTEGSTP, ABSERR, RELERR, 0, *5920)

C**   Numerically integrate; print results
5900  IF( TFINAL.GE.TINITIAL .AND. T+.01D0*INTEGSTP.GE.TFINAL) IPRINT=-7
      IF( TFINAL.LE.TINITIAL .AND. T+.01D0*INTEGSTP.LE.TFINAL) IPRINT=-7
      IF( IPRINT .LE. 0 ) THEN
        CALL IO(T)
        IF( IPRINT .EQ. -7 ) GOTO 5930
        IPRINT = PRINTINT
      ENDIF
      CALL KUTTA(EQNS1, 9, VAR, T, INTEGSTP, ABSERR, RELERR, 1, *5920)
      IPRINT = IPRINT - 1
      GOTO 5900

C**   Print message if numerical integration fails to converge
5920  CALL IO(T)
      WRITE(*, 6997)
      WRITE(21,6997)
      WRITE(22,6997)

C**   Inform user of input and output filename(s)
5930  WRITE(*,6999)

6021  FORMAT(1X,'FILE: 3dpendulum.1 ',//1X,'*** ',99A1,///,8X,'T',13X,'P
     &1X',12X,'P1Y',12X,'P1Z',12X,'P2X',12X,'P2Y',12X,'P2Z',/,5X,'(UNITS
     &)',8X,'(UNITS)',8X,'(UNITS)',8X,'(UNITS)',8X,'(UNITS)',8X,'(UNITS)
     &',8X,'(UNITS)',/)
6022  FORMAT(1X,'FILE: 3dpendulum.2 ',//1X,'*** ',99A1,///,8X,'T',13X,'R
     &X',13X,'RY',13X,'RZ',13X,'KE',13X,'PE',13X,'TE',/,5X,'(UNITS)',8X,
     &'(UNITS)',8X,'(UNITS)',8X,'(UNITS)',8X,'(UNITS)',8X,'(UNITS)',8X,'
     &(UNITS)',/)
6997  FORMAT(/7X,'Error: Numerical integration failed to converge',/)
6999  FORMAT(//1X,'Input is in the file 3dpendulum.in',//1X,'Output is i
     &n the file(s) 3dpendulum.i  (i=1,2)',//1X,'The output quantities a
     &nd associated files are listed in file 3dpendulum.dir',/)
7000  FORMAT(//,99A1,///)
7010  FORMAT( 1000(59X,E30.0,/) )
7011  FORMAT( 3(59X,E30.0,/), 1(59X,I30,/), 2(59X,E30.0,/) )
      STOP
7100  WRITE(*,*) 'Premature end of file while reading 3dpendulum.in '
7101  WRITE(*,*) 'Error while reading file 3dpendulum.in'
      STOP
      END


C**********************************************************************
      SUBROUTINE       EQNS1(T, VAR, VARp, BOUNDARY)
      IMPLICIT         DOUBLE PRECISION (A - Z)
      INTEGER          BOUNDARY
      DIMENSION        VAR(*), VARp(*)
      COMMON/CONSTNTS/ G,IX,IY,IZ,LA,LAO,M
      COMMON/VARIBLES/ Q1,Q2,Q3,Q4,Q5,Q6,U4,U5,U6
      COMMON/ALGBRAIC/ RX,RY,RZ,U1,U2,U3,Q1p,Q2p,Q3p,Q4p,Q5p,Q6p,U4p,U5p
     &,U6p,KE,P1X,P1Y,P1Z,P2X,P2Y,P2Z,PE,TE
      COMMON/MISCLLNS/ PI,DEGtoRAD,RADtoDEG,COEF(3,3),RHS(3)

C**   Update variables after integration step
      Q1 = VAR(1)
      Q2 = VAR(2)
      Q3 = VAR(3)
      Q4 = VAR(4)
      Q5 = VAR(5)
      Q6 = VAR(6)
      U4 = VAR(7)
      U5 = VAR(8)
      U6 = VAR(9)

      Q4p = U4 + TAN(Q5)*(SIN(Q4)*U5+COS(Q4)*U6)
      Q5p = COS(Q4)*U5 - SIN(Q4)*U6
      Q6p = (SIN(Q4)*U5+COS(Q4)*U6)/COS(Q5)

      COEF(1,1) = -IX
      COEF(1,2) = 0
      COEF(1,3) = 0
      COEF(2,1) = 0
      COEF(2,2) = -IY - M*LAO**2
      COEF(2,3) = 0
      COEF(3,1) = 0
      COEF(3,2) = 0
      COEF(3,3) = -IZ - M*LAO**2
      RHS(1) = -(IY-IZ)*U5*U6
      RHS(2) = G*LAO*M*COS(Q4)*COS(Q5) + (IX-IZ-M*LAO**2)*U4*U6
      RHS(3) = -G*LAO*M*SIN(Q4)*COS(Q5) - (IX-IY-M*LAO**2)*U4*U5
      CALL SOLVE(3,COEF,RHS,VARp)

C**   Update variables after uncoupling equations
      U4p = VARp(1)
      U5p = VARp(2)
      U6p = VARp(3)

C**   Update derivative array prior to integration step
      VARp(1) = Q1p
      VARp(2) = Q2p
      VARp(3) = Q3p
      VARp(4) = Q4p
      VARp(5) = Q5p
      VARp(6) = Q6p
      VARp(7) = U4p
      VARp(8) = U5p
      VARp(9) = U6p

      RETURN
      END


C**********************************************************************
      SUBROUTINE       IO(T)
      IMPLICIT         DOUBLE PRECISION (A - Z)
      INTEGER          ILOOP
      COMMON/CONSTNTS/ G,IX,IY,IZ,LA,LAO,M
      COMMON/VARIBLES/ Q1,Q2,Q3,Q4,Q5,Q6,U4,U5,U6
      COMMON/ALGBRAIC/ RX,RY,RZ,U1,U2,U3,Q1p,Q2p,Q3p,Q4p,Q5p,Q6p,U4p,U5p
     &,U6p,KE,P1X,P1Y,P1Z,P2X,P2Y,P2Z,PE,TE
      COMMON/MISCLLNS/ PI,DEGtoRAD,RADtoDEG,COEF(3,3),RHS(3)

C**   Evaluate output quantities
      RX = LAO*M*((SIN(Q4)*SIN(Q6)+SIN(Q5)*COS(Q4)*COS(Q6))*U4*U6-COS(Q5
     &)*COS(Q6)*(U5**2+U6**2)-(SIN(Q6)*COS(Q4)-SIN(Q4)*SIN(Q5)*COS(Q6))*
     &U4*U5-(SIN(Q4)*SIN(Q6)+SIN(Q5)*COS(Q4)*COS(Q6))*U5p-(SIN(Q6)*COS(Q
     &4)-SIN(Q4)*SIN(Q5)*COS(Q6))*U6p)
      RY = -LAO*M*(SIN(Q6)*COS(Q5)*(U5**2+U6**2)+(SIN(Q4)*COS(Q6)-SIN(Q5
     &)*SIN(Q6)*COS(Q4))*U4*U6-(COS(Q4)*COS(Q6)+SIN(Q4)*SIN(Q5)*SIN(Q6))
     &*U4*U5-(COS(Q4)*COS(Q6)+SIN(Q4)*SIN(Q5)*SIN(Q6))*U6p-(SIN(Q4)*COS(
     &Q6)-SIN(Q5)*SIN(Q6)*COS(Q4))*U5p)
      RZ = -M*(G-LAO*(SIN(Q4)*COS(Q5)*U4*U5+COS(Q4)*COS(Q5)*U4*U6+SIN(Q5
     &)*(U5**2+U6**2)+SIN(Q4)*COS(Q5)*U6p-COS(Q4)*COS(Q5)*U5p))
      KE = 0.5D0*IX*U4**2 + 0.5D0*IY*U5**2 + 0.5D0*IZ*U6**2 + 0.5D0*M*LA
     &O**2*(U5**2+U6**2)
      PE = -G*M*(Q3-LAO*SIN(Q5))
      TE = PE + KE
      P1X = Q1
      P1Y = Q2
      P1Z = Q3
      P2X = Q1 + LA*COS(Q5)*COS(Q6)
      P2Y = Q2 + LA*SIN(Q6)*COS(Q5)
      P2Z = Q3 - LA*SIN(Q5)

C**   Write output to screen and to output file(s)
      WRITE(*, 6020) T,P1X,P1Y,P1Z,P2X,P2Y,P2Z
      WRITE(21,6020) T,P1X,P1Y,P1Z,P2X,P2Y,P2Z
      WRITE(22,6020) T,RX,RY,RZ,KE,PE,TE

6020  FORMAT( 99(1X, 1PE14.6E3) )

      RETURN
      END


