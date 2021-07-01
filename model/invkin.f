C**   The name of this program is invkin.f
C**   Created by AUTOLEV 3.2 on Thu Jul 01 14:02:41 2021

      IMPLICIT         DOUBLE PRECISION (A - Z)
      INTEGER          ILOOP, IPRINT, PRINTINT
      CHARACTER        MESSAGE(99)
      EXTERNAL         EQNS1
      DIMENSION        VAR(12)
      COMMON/CONSTNTS/ T,G,IXA,IXB,IYA,IYB,IZA,IZB,LA,LAO,LB,LBO,MA,MB,P
     &X,PY,PZ,PXp,PYp,PZp
      COMMON/VARIBLES/ Q1,Q2,Q3,Q4,Q5,Q6,U1,U2,U3,U4,U5,U6
      COMMON/ALGBRAIC/ AMOMX,AMOMY,AMOMZ,KE,P1X,P1Y,P1Z,P2X,P2Y,P2Z,P3X,
     &P3Y,P3Z,PE,TE,VP2X,VP2Y,VP2Z,VP3X,VP3Y,VP3Z
      COMMON/MISCLLNS/ PI,DEGtoRAD,RADtoDEG,COEF(12,12),RHS(12),RHSGUESS
     &(12),GUESS(12)

C**   Open input and output files
      OPEN(UNIT=20, FILE='invkin.in', STATUS='OLD')
      OPEN(UNIT=21, FILE='invkin.1',  STATUS='UNKNOWN')

C**   Read message from input file
      READ(20,7000,END=7100,ERR=7101) (MESSAGE(ILOOP),ILOOP = 1,99)

C**   Read values of constants from input file
      READ(20,7010,END=7100,ERR=7101) T,G,IXA,IXB,IYA,IYB,IZA,IZB,LA,LAO
     &,LB,LBO,MA,MB,PX,PY,PZ,PXp,PYp,PZp

C**   Read the initial value of each variable from input file
      READ(20,7010,END=7100,ERR=7101) Q1,Q2,Q3,Q4,Q5,Q6,U1,U2,U3,U4,U5,U
     &6

C**   Read integration parameters from input file
      READ(20,7010,END=7100,ERR=7101) INTEGSTP,ABSERR,RELERR

C**   Write heading(s) to output file(s)
      WRITE(*, 6021) (MESSAGE(ILOOP), ILOOP = 1,99) 
      WRITE(21,6021) (MESSAGE(ILOOP), ILOOP = 1,99) 

C**   Degree to radian conversion
      PI       = 4*ATAN(1.0D0)
      DEGtoRAD = PI/180.0D0
      RADtoDEG = 180.0D0/PI

C**   Evaluate constants
      P1X = PX
      P1Y = PY
      P1Z = PZ


C**   Initialize independent variable and integrator array
      INTEGTAU = 0.0D0
      TAUFINAL = 1.0D0
      VAR(1) = Q1
      VAR(2) = Q2
      VAR(3) = Q3
      VAR(4) = Q4
      VAR(5) = Q5
      VAR(6) = Q6
      VAR(7) = U1
      VAR(8) = U2
      VAR(9) = U3
      VAR(10) = U4
      VAR(11) = U5
      VAR(12) = U6

C**   Store guess(es) and value(s) of function(s) at guess(es)
      CALL KUTTA(EQNS1,12,VAR,INTEGTAU,INTEGSTP,ABSERR,RELERR,0,*5920)
      DO 5800 ILOOP = 1, 12
         GUESS(ILOOP) = VAR(ILOOP)
      RHSGUESS(ILOOP) = RHS(ILOOP)
5800  CONTINUE

C**   Numerically integrate to obtain solution
      CALL KUTTA(EQNS1,12,VAR,INTEGTAU,INTEGSTP,ABSERR,RELERR,0,*5920)
5900  IF( INTEGTAU + .01D0*INTEGSTP .GE. TAUFINAL ) THEN
        CALL IO()
        GOTO 5930
      ENDIF
      CALL KUTTA(EQNS1,12,VAR,INTEGTAU,INTEGSTP,ABSERR,RELERR,1,*5920)
      GOTO 5900

C**   Print message if numerical integration fails to converge
5920  CALL IO()
      WRITE(*, 6997)
      WRITE(21,6997)

C**   Inform user of input and output filename(s)
5930  WRITE(*,6999)

6021  FORMAT(1X,'FILE: invkin.1',//1X,'*** ',99A1)
6997  FORMAT(/7X,'Error: Numerical integration failed to converge',/)
6999  FORMAT(//1X,'Input is in the file invkin.in',//1X,'Output is in th
     &e file(s) invkin.1',/)
7000  FORMAT(//,99A1,///)
7010  FORMAT( 1000(59X,E30.0,/) )
      STOP
7100  WRITE(*,*) 'Premature end of file while reading invkin.in '
7101  WRITE(*,*) 'Error while reading file invkin.in'
      STOP
      END


C**********************************************************************
      SUBROUTINE       EQNS1(INTEGTAU, VAR, VARp, BOUNDARY)
      IMPLICIT         DOUBLE PRECISION (A - Z)
      INTEGER          BOUNDARY
      DIMENSION        VAR(*), VARp(*)
      COMMON/CONSTNTS/ T,G,IXA,IXB,IYA,IYB,IZA,IZB,LA,LAO,LB,LBO,MA,MB,P
     &X,PY,PZ,PXp,PYp,PZp
      COMMON/VARIBLES/ Q1,Q2,Q3,Q4,Q5,Q6,U1,U2,U3,U4,U5,U6
      COMMON/ALGBRAIC/ AMOMX,AMOMY,AMOMZ,KE,P1X,P1Y,P1Z,P2X,P2Y,P2Z,P3X,
     &P3Y,P3Z,PE,TE,VP2X,VP2Y,VP2Z,VP3X,VP3Y,VP3Z
      COMMON/MISCLLNS/ PI,DEGtoRAD,RADtoDEG,COEF(12,12),RHS(12),RHSGUESS
     &(12),GUESS(12)

C**   Update variables after integration step
      Q1 = VAR(1)
      Q2 = VAR(2)
      Q3 = VAR(3)
      Q4 = VAR(4)
      Q5 = VAR(5)
      Q6 = VAR(6)
      U1 = VAR(7)
      U2 = VAR(8)
      U3 = VAR(9)
      U4 = VAR(10)
      U5 = VAR(11)
      U6 = VAR(12)

      P2X = PX + LA*COS(Q2)*COS(Q3)
      P2Y = PY + LA*(SIN(Q3)*COS(Q1)+SIN(Q1)*SIN(Q2)*COS(Q3))
      P2Z = PZ + LA*(SIN(Q1)*SIN(Q3)-SIN(Q2)*COS(Q1)*COS(Q3))
      P3X = PX + LA*COS(Q2)*COS(Q3) + LB*(COS(Q2)*COS(Q3)*COS(Q5)*COS(Q6
     &)+SIN(Q2)*(SIN(Q4)*SIN(Q6)-SIN(Q5)*COS(Q4)*COS(Q6))-SIN(Q3)*COS(Q2
     &)*(SIN(Q6)*COS(Q4)+SIN(Q4)*SIN(Q5)*COS(Q6)))
      P3Y = PY + LA*(SIN(Q3)*COS(Q1)+SIN(Q1)*SIN(Q2)*COS(Q3)) - LB*(SIN(
     &Q1)*COS(Q2)*(SIN(Q4)*SIN(Q6)-SIN(Q5)*COS(Q4)*COS(Q6))-COS(Q5)*COS(
     &Q6)*(SIN(Q3)*COS(Q1)+SIN(Q1)*SIN(Q2)*COS(Q3))-(SIN(Q6)*COS(Q4)+SIN
     &(Q4)*SIN(Q5)*COS(Q6))*(COS(Q1)*COS(Q3)-SIN(Q1)*SIN(Q2)*SIN(Q3)))
      P3Z = PZ + LA*(SIN(Q1)*SIN(Q3)-SIN(Q2)*COS(Q1)*COS(Q3)) + LB*(COS(
     &Q1)*COS(Q2)*(SIN(Q4)*SIN(Q6)-SIN(Q5)*COS(Q4)*COS(Q6))+COS(Q5)*COS(
     &Q6)*(SIN(Q1)*SIN(Q3)-SIN(Q2)*COS(Q1)*COS(Q3))+(SIN(Q1)*COS(Q3)+SIN
     &(Q2)*SIN(Q3)*COS(Q1))*(SIN(Q6)*COS(Q4)+SIN(Q4)*SIN(Q5)*COS(Q6)))
      VP2X = PXp - LA*SIN(Q2)*U2 - LA*SIN(Q3)*COS(Q2)*U3
      VP2Y = PYp + LA*SIN(Q1)*COS(Q2)*U2 + LA*(COS(Q1)*COS(Q3)-SIN(Q1)*S
     &IN(Q2)*SIN(Q3))*U3
      VP2Z = PZp + LA*(SIN(Q1)*COS(Q3)+SIN(Q2)*SIN(Q3)*COS(Q1))*U3 - LA*
     &COS(Q1)*COS(Q2)*U2
      VP3X = PXp + LB*(SIN(Q2)*COS(Q4)*COS(Q5)+SIN(Q5)*COS(Q2)*COS(Q3)+S
     &IN(Q3)*SIN(Q4)*COS(Q2)*COS(Q5))*(SIN(Q6)*COS(Q5)*U1-U5-(SIN(Q4)*CO
     &S(Q6)+SIN(Q5)*SIN(Q6)*COS(Q4))*U3-(COS(Q4)*COS(Q6)-SIN(Q4)*SIN(Q5)
     &*SIN(Q6))*U2) - LA*SIN(Q2)*U2 - LA*SIN(Q3)*COS(Q2)*U3 - LB*(SIN(Q2
     &)*(SIN(Q4)*COS(Q6)+SIN(Q5)*SIN(Q6)*COS(Q4))-SIN(Q6)*COS(Q2)*COS(Q3
     &)*COS(Q5)-SIN(Q3)*COS(Q2)*(COS(Q4)*COS(Q6)-SIN(Q4)*SIN(Q5)*SIN(Q6)
     &))*(SIN(Q4)*COS(Q5)*U2-U6-SIN(Q5)*U1-COS(Q4)*COS(Q5)*U3)
      VP3Y = PYp + LA*SIN(Q1)*COS(Q2)*U2 + LA*(COS(Q1)*COS(Q3)-SIN(Q1)*S
     &IN(Q2)*SIN(Q3))*U3 + LB*(SIN(Q1)*COS(Q2)*(SIN(Q4)*COS(Q6)+SIN(Q5)*
     &SIN(Q6)*COS(Q4))+SIN(Q6)*COS(Q5)*(SIN(Q3)*COS(Q1)+SIN(Q1)*SIN(Q2)*
     &COS(Q3))-(COS(Q1)*COS(Q3)-SIN(Q1)*SIN(Q2)*SIN(Q3))*(COS(Q4)*COS(Q6
     &)-SIN(Q4)*SIN(Q5)*SIN(Q6)))*(SIN(Q4)*COS(Q5)*U2-U6-SIN(Q5)*U1-COS(
     &Q4)*COS(Q5)*U3) + LB*(SIN(Q5)*(SIN(Q3)*COS(Q1)+SIN(Q1)*SIN(Q2)*COS
     &(Q3))-SIN(Q1)*COS(Q2)*COS(Q4)*COS(Q5)-SIN(Q4)*COS(Q5)*(COS(Q1)*COS
     &(Q3)-SIN(Q1)*SIN(Q2)*SIN(Q3)))*(SIN(Q6)*COS(Q5)*U1-U5-(SIN(Q4)*COS
     &(Q6)+SIN(Q5)*SIN(Q6)*COS(Q4))*U3-(COS(Q4)*COS(Q6)-SIN(Q4)*SIN(Q5)*
     &SIN(Q6))*U2)
      VP3Z = PZp + LA*(SIN(Q1)*COS(Q3)+SIN(Q2)*SIN(Q3)*COS(Q1))*U3 + LB*
     &(COS(Q1)*COS(Q2)*COS(Q4)*COS(Q5)+SIN(Q5)*(SIN(Q1)*SIN(Q3)-SIN(Q2)*
     &COS(Q1)*COS(Q3))-SIN(Q4)*COS(Q5)*(SIN(Q1)*COS(Q3)+SIN(Q2)*SIN(Q3)*
     &COS(Q1)))*(SIN(Q6)*COS(Q5)*U1-U5-(SIN(Q4)*COS(Q6)+SIN(Q5)*SIN(Q6)*
     &COS(Q4))*U3-(COS(Q4)*COS(Q6)-SIN(Q4)*SIN(Q5)*SIN(Q6))*U2) + LB*(SI
     &N(Q6)*COS(Q5)*(SIN(Q1)*SIN(Q3)-SIN(Q2)*COS(Q1)*COS(Q3))-COS(Q1)*CO
     &S(Q2)*(SIN(Q4)*COS(Q6)+SIN(Q5)*SIN(Q6)*COS(Q4))-(SIN(Q1)*COS(Q3)+S
     &IN(Q2)*SIN(Q3)*COS(Q1))*(COS(Q4)*COS(Q6)-SIN(Q4)*SIN(Q5)*SIN(Q6)))
     &*(SIN(Q4)*COS(Q5)*U2-U6-SIN(Q5)*U1-COS(Q4)*COS(Q5)*U3) - LA*COS(Q1
     &)*COS(Q2)*U2

      COEF(1,1) = 0
      COEF(1,2) = 0
      COEF(1,3) = 0
      COEF(1,4) = 0
      COEF(1,5) = 0
      COEF(1,6) = 0
      COEF(1,7) = 0
      COEF(1,8) = 0
      COEF(1,9) = 0
      COEF(1,10) = 0
      COEF(1,11) = 0
      COEF(1,12) = 0
      COEF(2,1) = 0
      COEF(2,2) = 0
      COEF(2,3) = 0
      COEF(2,4) = 0
      COEF(2,5) = 0
      COEF(2,6) = 0
      COEF(2,7) = 0
      COEF(2,8) = 0
      COEF(2,9) = 0
      COEF(2,10) = 0
      COEF(2,11) = 0
      COEF(2,12) = 0
      COEF(3,1) = 0
      COEF(3,2) = 0
      COEF(3,3) = 0
      COEF(3,4) = 0
      COEF(3,5) = 0
      COEF(3,6) = 0
      COEF(3,7) = 0
      COEF(3,8) = 0
      COEF(3,9) = 0
      COEF(3,10) = 0
      COEF(3,11) = 0
      COEF(3,12) = 0
      COEF(4,1) = 0
      COEF(4,2) = 0
      COEF(4,3) = 0
      COEF(4,4) = 0
      COEF(4,5) = 0
      COEF(4,6) = 0
      COEF(4,7) = 0
      COEF(4,8) = 0
      COEF(4,9) = 0
      COEF(4,10) = 0
      COEF(4,11) = 0
      COEF(4,12) = 0
      COEF(5,1) = 0
      COEF(5,2) = 0
      COEF(5,3) = 0
      COEF(5,4) = 0
      COEF(5,5) = 0
      COEF(5,6) = 0
      COEF(5,7) = 0
      COEF(5,8) = 0
      COEF(5,9) = 0
      COEF(5,10) = 0
      COEF(5,11) = 0
      COEF(5,12) = 0
      COEF(6,1) = 0
      COEF(6,2) = 0
      COEF(6,3) = 0
      COEF(6,4) = 0
      COEF(6,5) = 0
      COEF(6,6) = 0
      COEF(6,7) = 0
      COEF(6,8) = 0
      COEF(6,9) = 0
      COEF(6,10) = 0
      COEF(6,11) = 0
      COEF(6,12) = 0
      COEF(7,1) = 0
      COEF(7,2) = 0
      COEF(7,3) = 0
      COEF(7,4) = 0
      COEF(7,5) = 0
      COEF(7,6) = 0
      COEF(7,7) = 0
      COEF(7,8) = 0
      COEF(7,9) = 0
      COEF(7,10) = 0
      COEF(7,11) = 0
      COEF(7,12) = 0
      COEF(8,1) = 0
      COEF(8,2) = 0
      COEF(8,3) = 0
      COEF(8,4) = 0
      COEF(8,5) = 0
      COEF(8,6) = 0
      COEF(8,7) = 0
      COEF(8,8) = 0
      COEF(8,9) = 0
      COEF(8,10) = 0
      COEF(8,11) = 0
      COEF(8,12) = 0
      COEF(9,1) = 0
      COEF(9,2) = 0
      COEF(9,3) = 0
      COEF(9,4) = 0
      COEF(9,5) = 0
      COEF(9,6) = 0
      COEF(9,7) = 0
      COEF(9,8) = 0
      COEF(9,9) = 0
      COEF(9,10) = 0
      COEF(9,11) = 0
      COEF(9,12) = 0
      COEF(10,1) = 0
      COEF(10,2) = 0
      COEF(10,3) = 0
      COEF(10,4) = 0
      COEF(10,5) = 0
      COEF(10,6) = 0
      COEF(10,7) = 0
      COEF(10,8) = 0
      COEF(10,9) = 0
      COEF(10,10) = 0
      COEF(10,11) = 0
      COEF(10,12) = 0
      COEF(11,1) = 0
      COEF(11,2) = 0
      COEF(11,3) = 0
      COEF(11,4) = 0
      COEF(11,5) = 0
      COEF(11,6) = 0
      COEF(11,7) = 0
      COEF(11,8) = 0
      COEF(11,9) = 0
      COEF(11,10) = 0
      COEF(11,11) = 0
      COEF(11,12) = 0
      COEF(12,1) = 0
      COEF(12,2) = 0
      COEF(12,3) = 0
      COEF(12,4) = 0
      COEF(12,5) = 0
      COEF(12,6) = 0
      COEF(12,7) = 0
      COEF(12,8) = 0
      COEF(12,9) = 0
      COEF(12,10) = 0
      COEF(12,11) = 0
      COEF(12,12) = 0
      RHS(1) = PX + LA*COS(Q2)*COS(Q3) - P2X
      RHS(2) = PY + LA*(SIN(Q3)*COS(Q1)+SIN(Q1)*SIN(Q2)*COS(Q3)) - P2Y
      RHS(3) = PZ + LA*(SIN(Q1)*SIN(Q3)-SIN(Q2)*COS(Q1)*COS(Q3)) - P2Z
      RHS(4) = PX + LA*COS(Q2)*COS(Q3) + LB*(COS(Q2)*COS(Q3)*COS(Q5)*COS
     &(Q6)+SIN(Q2)*(SIN(Q4)*SIN(Q6)-SIN(Q5)*COS(Q4)*COS(Q6))-SIN(Q3)*COS
     &(Q2)*(SIN(Q6)*COS(Q4)+SIN(Q4)*SIN(Q5)*COS(Q6))) - P3X
      RHS(5) = PY + LA*(SIN(Q3)*COS(Q1)+SIN(Q1)*SIN(Q2)*COS(Q3)) - LB*(S
     &IN(Q1)*COS(Q2)*(SIN(Q4)*SIN(Q6)-SIN(Q5)*COS(Q4)*COS(Q6))-COS(Q5)*C
     &OS(Q6)*(SIN(Q3)*COS(Q1)+SIN(Q1)*SIN(Q2)*COS(Q3))-(SIN(Q6)*COS(Q4)+
     &SIN(Q4)*SIN(Q5)*COS(Q6))*(COS(Q1)*COS(Q3)-SIN(Q1)*SIN(Q2)*SIN(Q3))
     &) - P3Y
      RHS(6) = PZ + LA*(SIN(Q1)*SIN(Q3)-SIN(Q2)*COS(Q1)*COS(Q3)) + LB*(C
     &OS(Q1)*COS(Q2)*(SIN(Q4)*SIN(Q6)-SIN(Q5)*COS(Q4)*COS(Q6))+COS(Q5)*C
     &OS(Q6)*(SIN(Q1)*SIN(Q3)-SIN(Q2)*COS(Q1)*COS(Q3))+(SIN(Q1)*COS(Q3)+
     &SIN(Q2)*SIN(Q3)*COS(Q1))*(SIN(Q6)*COS(Q4)+SIN(Q4)*SIN(Q5)*COS(Q6))
     &) - P3Z
      RHS(7) = PXp - LA*SIN(Q2)*U2 - LA*SIN(Q3)*COS(Q2)*U3 - VP2X
      RHS(8) = PYp + LA*SIN(Q1)*COS(Q2)*U2 + LA*(COS(Q1)*COS(Q3)-SIN(Q1)
     &*SIN(Q2)*SIN(Q3))*U3 - VP2Y
      RHS(9) = PZp + LA*(SIN(Q1)*COS(Q3)+SIN(Q2)*SIN(Q3)*COS(Q1))*U3 - L
     &A*COS(Q1)*COS(Q2)*U2 - VP2Z
      RHS(10) = PXp + LB*(SIN(Q2)*COS(Q4)*COS(Q5)+SIN(Q5)*COS(Q2)*COS(Q3
     &)+SIN(Q3)*SIN(Q4)*COS(Q2)*COS(Q5))*(SIN(Q6)*COS(Q5)*U1-U5-(SIN(Q4)
     &*COS(Q6)+SIN(Q5)*SIN(Q6)*COS(Q4))*U3-(COS(Q4)*COS(Q6)-SIN(Q4)*SIN(
     &Q5)*SIN(Q6))*U2) - LA*SIN(Q2)*U2 - LA*SIN(Q3)*COS(Q2)*U3 - LB*(SIN
     &(Q2)*(SIN(Q4)*COS(Q6)+SIN(Q5)*SIN(Q6)*COS(Q4))-SIN(Q6)*COS(Q2)*COS
     &(Q3)*COS(Q5)-SIN(Q3)*COS(Q2)*(COS(Q4)*COS(Q6)-SIN(Q4)*SIN(Q5)*SIN(
     &Q6)))*(SIN(Q4)*COS(Q5)*U2-U6-SIN(Q5)*U1-COS(Q4)*COS(Q5)*U3) - VP3X
      RHS(11) = PYp + LA*SIN(Q1)*COS(Q2)*U2 + LA*(COS(Q1)*COS(Q3)-SIN(Q1
     &)*SIN(Q2)*SIN(Q3))*U3 + LB*(SIN(Q1)*COS(Q2)*(SIN(Q4)*COS(Q6)+SIN(Q
     &5)*SIN(Q6)*COS(Q4))+SIN(Q6)*COS(Q5)*(SIN(Q3)*COS(Q1)+SIN(Q1)*SIN(Q
     &2)*COS(Q3))-(COS(Q1)*COS(Q3)-SIN(Q1)*SIN(Q2)*SIN(Q3))*(COS(Q4)*COS
     &(Q6)-SIN(Q4)*SIN(Q5)*SIN(Q6)))*(SIN(Q4)*COS(Q5)*U2-U6-SIN(Q5)*U1-C
     &OS(Q4)*COS(Q5)*U3) + LB*(SIN(Q5)*(SIN(Q3)*COS(Q1)+SIN(Q1)*SIN(Q2)*
     &COS(Q3))-SIN(Q1)*COS(Q2)*COS(Q4)*COS(Q5)-SIN(Q4)*COS(Q5)*(COS(Q1)*
     &COS(Q3)-SIN(Q1)*SIN(Q2)*SIN(Q3)))*(SIN(Q6)*COS(Q5)*U1-U5-(SIN(Q4)*
     &COS(Q6)+SIN(Q5)*SIN(Q6)*COS(Q4))*U3-(COS(Q4)*COS(Q6)-SIN(Q4)*SIN(Q
     &5)*SIN(Q6))*U2) - VP3Y
      RHS(12) = PZp + LA*(SIN(Q1)*COS(Q3)+SIN(Q2)*SIN(Q3)*COS(Q1))*U3 + 
     &LB*(COS(Q1)*COS(Q2)*COS(Q4)*COS(Q5)+SIN(Q5)*(SIN(Q1)*SIN(Q3)-SIN(Q
     &2)*COS(Q1)*COS(Q3))-SIN(Q4)*COS(Q5)*(SIN(Q1)*COS(Q3)+SIN(Q2)*SIN(Q
     &3)*COS(Q1)))*(SIN(Q6)*COS(Q5)*U1-U5-(SIN(Q4)*COS(Q6)+SIN(Q5)*SIN(Q
     &6)*COS(Q4))*U3-(COS(Q4)*COS(Q6)-SIN(Q4)*SIN(Q5)*SIN(Q6))*U2) + LB*
     &(SIN(Q6)*COS(Q5)*(SIN(Q1)*SIN(Q3)-SIN(Q2)*COS(Q1)*COS(Q3))-COS(Q1)
     &*COS(Q2)*(SIN(Q4)*COS(Q6)+SIN(Q5)*SIN(Q6)*COS(Q4))-(SIN(Q1)*COS(Q3
     &)+SIN(Q2)*SIN(Q3)*COS(Q1))*(COS(Q4)*COS(Q6)-SIN(Q4)*SIN(Q5)*SIN(Q6
     &)))*(SIN(Q4)*COS(Q5)*U2-U6-SIN(Q5)*U1-COS(Q4)*COS(Q5)*U3) - LA*COS
     &(Q1)*COS(Q2)*U2 - VP3Z
      CALL SOLVE(12,COEF,RHSGUESS,VARp)

      RETURN
      END


C**********************************************************************
      SUBROUTINE       IO()
      IMPLICIT         DOUBLE PRECISION (A - Z)
      INTEGER          ILOOP
      COMMON/CONSTNTS/ T,G,IXA,IXB,IYA,IYB,IZA,IZB,LA,LAO,LB,LBO,MA,MB,P
     &X,PY,PZ,PXp,PYp,PZp
      COMMON/VARIBLES/ Q1,Q2,Q3,Q4,Q5,Q6,U1,U2,U3,U4,U5,U6
      COMMON/ALGBRAIC/ AMOMX,AMOMY,AMOMZ,KE,P1X,P1Y,P1Z,P2X,P2Y,P2Z,P3X,
     &P3Y,P3Z,PE,TE,VP2X,VP2Y,VP2Z,VP3X,VP3Y,VP3Z
      COMMON/MISCLLNS/ PI,DEGtoRAD,RADtoDEG,COEF(12,12),RHS(12),RHSGUESS
     &(12),GUESS(12)

C**   Evaluate output quantities
      KE = 0.5D0*IXA*U1**2 + 0.5D0*IYA*U2**2 + 0.5D0*IZA*U3**2 + 0.5D0*I
     &ZB*SIN(Q4)*COS(Q5)*U2*(SIN(Q4)*COS(Q5)*U2-U6-SIN(Q5)*U1-COS(Q4)*CO
     &S(Q5)*U3) + 0.5D0*IXB*U4*(U4+COS(Q5)*COS(Q6)*U1+(SIN(Q6)*COS(Q4)+S
     &IN(Q4)*SIN(Q5)*COS(Q6))*U2+(SIN(Q4)*SIN(Q6)-SIN(Q5)*COS(Q4)*COS(Q6
     &))*U3) + 0.5D0*IXB*COS(Q5)*COS(Q6)*U1*(U4+COS(Q5)*COS(Q6)*U1+(SIN(
     &Q6)*COS(Q4)+SIN(Q4)*SIN(Q5)*COS(Q6))*U2+(SIN(Q4)*SIN(Q6)-SIN(Q5)*C
     &OS(Q4)*COS(Q6))*U3) + 0.5D0*IYB*SIN(Q6)*COS(Q5)*U1*(SIN(Q6)*COS(Q5
     &)*U1-U5-(SIN(Q4)*COS(Q6)+SIN(Q5)*SIN(Q6)*COS(Q4))*U3-(COS(Q4)*COS(
     &Q6)-SIN(Q4)*SIN(Q5)*SIN(Q6))*U2) + 0.5D0*IXB*(SIN(Q6)*COS(Q4)+SIN(
     &Q4)*SIN(Q5)*COS(Q6))*U2*(U4+COS(Q5)*COS(Q6)*U1+(SIN(Q6)*COS(Q4)+SI
     &N(Q4)*SIN(Q5)*COS(Q6))*U2+(SIN(Q4)*SIN(Q6)-SIN(Q5)*COS(Q4)*COS(Q6)
     &)*U3) + 0.5D0*IXB*(SIN(Q4)*SIN(Q6)-SIN(Q5)*COS(Q4)*COS(Q6))*U3*(U4
     &+COS(Q5)*COS(Q6)*U1+(SIN(Q6)*COS(Q4)+SIN(Q4)*SIN(Q5)*COS(Q6))*U2+(
     &SIN(Q4)*SIN(Q6)-SIN(Q5)*COS(Q4)*COS(Q6))*U3) - 0.5D0*IZB*U6*(SIN(Q
     &4)*COS(Q5)*U2-U6-SIN(Q5)*U1-COS(Q4)*COS(Q5)*U3) - 0.5D0*IZB*SIN(Q5
     &)*U1*(SIN(Q4)*COS(Q5)*U2-U6-SIN(Q5)*U1-COS(Q4)*COS(Q5)*U3) - 0.5D0
     &*IZB*COS(Q4)*COS(Q5)*U3*(SIN(Q4)*COS(Q5)*U2-U6-SIN(Q5)*U1-COS(Q4)*
     &COS(Q5)*U3) - 0.5D0*IYB*U5*(SIN(Q6)*COS(Q5)*U1-U5-(SIN(Q4)*COS(Q6)
     &+SIN(Q5)*SIN(Q6)*COS(Q4))*U3-(COS(Q4)*COS(Q6)-SIN(Q4)*SIN(Q5)*SIN(
     &Q6))*U2) - 0.5D0*IYB*(SIN(Q4)*COS(Q6)+SIN(Q5)*SIN(Q6)*COS(Q4))*U3*
     &(SIN(Q6)*COS(Q5)*U1-U5-(SIN(Q4)*COS(Q6)+SIN(Q5)*SIN(Q6)*COS(Q4))*U
     &3-(COS(Q4)*COS(Q6)-SIN(Q4)*SIN(Q5)*SIN(Q6))*U2) - 0.5D0*IYB*(COS(Q
     &4)*COS(Q6)-SIN(Q4)*SIN(Q5)*SIN(Q6))*U2*(SIN(Q6)*COS(Q5)*U1-U5-(SIN
     &(Q4)*COS(Q6)+SIN(Q5)*SIN(Q6)*COS(Q4))*U3-(COS(Q4)*COS(Q6)-SIN(Q4)*
     &SIN(Q5)*SIN(Q6))*U2) - 0.5D0*MA*(2*LAO*PXp*SIN(Q2)*U2+2*LAO*PXp*SI
     &N(Q3)*COS(Q2)*U3+2*LAO*PZp*COS(Q1)*COS(Q2)*U2-PXp**2-PYp**2-PZp**2
     &-LAO**2*U2**2-LAO**2*U3**2-2*LAO*PYp*SIN(Q1)*COS(Q2)*U2-2*LAO*PZp*
     &(SIN(Q1)*COS(Q3)+SIN(Q2)*SIN(Q3)*COS(Q1))*U3-2*LAO*PYp*(COS(Q1)*CO
     &S(Q3)-SIN(Q1)*SIN(Q2)*SIN(Q3))*U3) - 0.5D0*MB*(2*LA*PXp*SIN(Q2)*U2
     &+2*LA*PXp*SIN(Q3)*COS(Q2)*U3+2*LA*PZp*COS(Q1)*COS(Q2)*U2+2*LA*LBO*
     &(COS(Q4)*COS(Q6)-SIN(Q4)*SIN(Q5)*SIN(Q6))*U3*(SIN(Q4)*COS(Q5)*U2-U
     &6-SIN(Q5)*U1-COS(Q4)*COS(Q5)*U3)+2*LA*LBO*SIN(Q4)*COS(Q5)*U3*(SIN(
     &Q6)*COS(Q5)*U1-U5-(SIN(Q4)*COS(Q6)+SIN(Q5)*SIN(Q6)*COS(Q4))*U3-(CO
     &S(Q4)*COS(Q6)-SIN(Q4)*SIN(Q5)*SIN(Q6))*U2)+2*LA*LBO*COS(Q4)*COS(Q5
     &)*U2*(SIN(Q6)*COS(Q5)*U1-U5-(SIN(Q4)*COS(Q6)+SIN(Q5)*SIN(Q6)*COS(Q
     &4))*U3-(COS(Q4)*COS(Q6)-SIN(Q4)*SIN(Q5)*SIN(Q6))*U2)+2*LBO*PXp*(SI
     &N(Q2)*(SIN(Q4)*COS(Q6)+SIN(Q5)*SIN(Q6)*COS(Q4))-SIN(Q6)*COS(Q2)*CO
     &S(Q3)*COS(Q5)-SIN(Q3)*COS(Q2)*(COS(Q4)*COS(Q6)-SIN(Q4)*SIN(Q5)*SIN
     &(Q6)))*(SIN(Q4)*COS(Q5)*U2-U6-SIN(Q5)*U1-COS(Q4)*COS(Q5)*U3)-PXp**
     &2-PYp**2-PZp**2-LA**2*U2**2-LA**2*U3**2-2*LA*PYp*SIN(Q1)*COS(Q2)*U
     &2-2*LA*PZp*(SIN(Q1)*COS(Q3)+SIN(Q2)*SIN(Q3)*COS(Q1))*U3-2*LA*PYp*(
     &COS(Q1)*COS(Q3)-SIN(Q1)*SIN(Q2)*SIN(Q3))*U3-LBO**2*(SIN(Q4)*COS(Q5
     &)*U2-U6-SIN(Q5)*U1-COS(Q4)*COS(Q5)*U3)**2-2*LA*LBO*(SIN(Q4)*COS(Q6
     &)+SIN(Q5)*SIN(Q6)*COS(Q4))*U2*(SIN(Q4)*COS(Q5)*U2-U6-SIN(Q5)*U1-CO
     &S(Q4)*COS(Q5)*U3)-LBO**2*(SIN(Q6)*COS(Q5)*U1-U5-(SIN(Q4)*COS(Q6)+S
     &IN(Q5)*SIN(Q6)*COS(Q4))*U3-(COS(Q4)*COS(Q6)-SIN(Q4)*SIN(Q5)*SIN(Q6
     &))*U2)**2-2*LBO*PXp*(SIN(Q2)*COS(Q4)*COS(Q5)+SIN(Q5)*COS(Q2)*COS(Q
     &3)+SIN(Q3)*SIN(Q4)*COS(Q2)*COS(Q5))*(SIN(Q6)*COS(Q5)*U1-U5-(SIN(Q4
     &)*COS(Q6)+SIN(Q5)*SIN(Q6)*COS(Q4))*U3-(COS(Q4)*COS(Q6)-SIN(Q4)*SIN
     &(Q5)*SIN(Q6))*U2)-2*LBO*PYp*(SIN(Q1)*COS(Q2)*(SIN(Q4)*COS(Q6)+SIN(
     &Q5)*SIN(Q6)*COS(Q4))+SIN(Q6)*COS(Q5)*(SIN(Q3)*COS(Q1)+SIN(Q1)*SIN(
     &Q2)*COS(Q3))-(COS(Q1)*COS(Q3)-SIN(Q1)*SIN(Q2)*SIN(Q3))*(COS(Q4)*CO
     &S(Q6)-SIN(Q4)*SIN(Q5)*SIN(Q6)))*(SIN(Q4)*COS(Q5)*U2-U6-SIN(Q5)*U1-
     &COS(Q4)*COS(Q5)*U3)-2*LBO*PZp*(COS(Q1)*COS(Q2)*COS(Q4)*COS(Q5)+SIN
     &(Q5)*(SIN(Q1)*SIN(Q3)-SIN(Q2)*COS(Q1)*COS(Q3))-SIN(Q4)*COS(Q5)*(SI
     &N(Q1)*COS(Q3)+SIN(Q2)*SIN(Q3)*COS(Q1)))*(SIN(Q6)*COS(Q5)*U1-U5-(SI
     &N(Q4)*COS(Q6)+SIN(Q5)*SIN(Q6)*COS(Q4))*U3-(COS(Q4)*COS(Q6)-SIN(Q4)
     &*SIN(Q5)*SIN(Q6))*U2)-2*LBO*PYp*(SIN(Q5)*(SIN(Q3)*COS(Q1)+SIN(Q1)*
     &SIN(Q2)*COS(Q3))-SIN(Q1)*COS(Q2)*COS(Q4)*COS(Q5)-SIN(Q4)*COS(Q5)*(
     &COS(Q1)*COS(Q3)-SIN(Q1)*SIN(Q2)*SIN(Q3)))*(SIN(Q6)*COS(Q5)*U1-U5-(
     &SIN(Q4)*COS(Q6)+SIN(Q5)*SIN(Q6)*COS(Q4))*U3-(COS(Q4)*COS(Q6)-SIN(Q
     &4)*SIN(Q5)*SIN(Q6))*U2)-2*LBO*PZp*(SIN(Q6)*COS(Q5)*(SIN(Q1)*SIN(Q3
     &)-SIN(Q2)*COS(Q1)*COS(Q3))-COS(Q1)*COS(Q2)*(SIN(Q4)*COS(Q6)+SIN(Q5
     &)*SIN(Q6)*COS(Q4))-(SIN(Q1)*COS(Q3)+SIN(Q2)*SIN(Q3)*COS(Q1))*(COS(
     &Q4)*COS(Q6)-SIN(Q4)*SIN(Q5)*SIN(Q6)))*(SIN(Q4)*COS(Q5)*U2-U6-SIN(Q
     &5)*U1-COS(Q4)*COS(Q5)*U3))
      PE = -G*((MA+MB)*(LAO+LBO+PZ)+(LA*MB+LAO*MA)*(SIN(Q1)*SIN(Q3)-SIN(
     &Q2)*COS(Q1)*COS(Q3))+LBO*MB*(COS(Q1)*COS(Q2)*(SIN(Q4)*SIN(Q6)-SIN(
     &Q5)*COS(Q4)*COS(Q6))+COS(Q5)*COS(Q6)*(SIN(Q1)*SIN(Q3)-SIN(Q2)*COS(
     &Q1)*COS(Q3))+(SIN(Q1)*COS(Q3)+SIN(Q2)*SIN(Q3)*COS(Q1))*(SIN(Q6)*CO
     &S(Q4)+SIN(Q4)*SIN(Q5)*COS(Q6))))
      TE = PE + KE
      AMOMX = LAO*MA*(PZp*(SIN(Q3)*COS(Q1)+SIN(Q1)*SIN(Q2)*COS(Q3))-PYp*
     &(SIN(Q1)*SIN(Q3)-SIN(Q2)*COS(Q1)*COS(Q3))) + MB*(LA*PZp*(SIN(Q3)*C
     &OS(Q1)+SIN(Q1)*SIN(Q2)*COS(Q3))-LA*PYp*(SIN(Q1)*SIN(Q3)-SIN(Q2)*CO
     &S(Q1)*COS(Q3))-LBO*PYp*(COS(Q1)*COS(Q2)*(SIN(Q4)*SIN(Q6)-SIN(Q5)*C
     &OS(Q4)*COS(Q6))+COS(Q5)*COS(Q6)*(SIN(Q1)*SIN(Q3)-SIN(Q2)*COS(Q1)*C
     &OS(Q3))+(SIN(Q1)*COS(Q3)+SIN(Q2)*SIN(Q3)*COS(Q1))*(SIN(Q6)*COS(Q4)
     &+SIN(Q4)*SIN(Q5)*COS(Q6)))-LBO*PZp*(SIN(Q1)*COS(Q2)*(SIN(Q4)*SIN(Q
     &6)-SIN(Q5)*COS(Q4)*COS(Q6))-COS(Q5)*COS(Q6)*(SIN(Q3)*COS(Q1)+SIN(Q
     &1)*SIN(Q2)*COS(Q3))-(SIN(Q6)*COS(Q4)+SIN(Q4)*SIN(Q5)*COS(Q6))*(COS
     &(Q1)*COS(Q3)-SIN(Q1)*SIN(Q2)*SIN(Q3)))) + SIN(Q2)*(IZA+MA*LAO**2+L
     &A*MB*(LA+LBO*COS(Q5)*COS(Q6)))*U3 + COS(Q2)*COS(Q3)*(IXA*U1-LA*LBO
     &*MB*((SIN(Q6)*COS(Q4)+SIN(Q4)*SIN(Q5)*COS(Q6))*U2+(SIN(Q4)*SIN(Q6)
     &-SIN(Q5)*COS(Q4)*COS(Q6))*U3)) + (SIN(Q2)*COS(Q4)*COS(Q5)+SIN(Q5)*
     &COS(Q2)*COS(Q3)+SIN(Q3)*SIN(Q4)*COS(Q2)*COS(Q5))*(IZB*U6+IZB*SIN(Q
     &5)*U1+IZB*COS(Q4)*COS(Q5)*U3-IZB*SIN(Q4)*COS(Q5)*U2-LBO*MB*(LBO+LA
     &*COS(Q5)*COS(Q6))*(SIN(Q4)*COS(Q5)*U2-U6-SIN(Q5)*U1-COS(Q4)*COS(Q5
     &)*U3)) + (SIN(Q2)*(SIN(Q4)*COS(Q6)+SIN(Q5)*SIN(Q6)*COS(Q4))-SIN(Q6
     &)*COS(Q2)*COS(Q3)*COS(Q5)-SIN(Q3)*COS(Q2)*(COS(Q4)*COS(Q6)-SIN(Q4)
     &*SIN(Q5)*SIN(Q6)))*(IYB*U5+IYB*(SIN(Q4)*COS(Q6)+SIN(Q5)*SIN(Q6)*CO
     &S(Q4))*U3+IYB*(COS(Q4)*COS(Q6)-SIN(Q4)*SIN(Q5)*SIN(Q6))*U2-IYB*SIN
     &(Q6)*COS(Q5)*U1-LBO*MB*(LBO+LA*COS(Q5)*COS(Q6))*(SIN(Q6)*COS(Q5)*U
     &1-U5-(SIN(Q4)*COS(Q6)+SIN(Q5)*SIN(Q6)*COS(Q4))*U3-(COS(Q4)*COS(Q6)
     &-SIN(Q4)*SIN(Q5)*SIN(Q6))*U2)) + (COS(Q2)*COS(Q3)*COS(Q5)*COS(Q6)+
     &SIN(Q2)*(SIN(Q4)*SIN(Q6)-SIN(Q5)*COS(Q4)*COS(Q6))-SIN(Q3)*COS(Q2)*
     &(SIN(Q6)*COS(Q4)+SIN(Q4)*SIN(Q5)*COS(Q6)))*(IXB*U4+IXB*COS(Q5)*COS
     &(Q6)*U1+IXB*(SIN(Q6)*COS(Q4)+SIN(Q4)*SIN(Q5)*COS(Q6))*U2+IXB*(SIN(
     &Q4)*SIN(Q6)-SIN(Q5)*COS(Q4)*COS(Q6))*U3+LA*LBO*MB*(SIN(Q5)*(SIN(Q4
     &)*COS(Q5)*U2-U6-SIN(Q5)*U1-COS(Q4)*COS(Q5)*U3)-SIN(Q6)*COS(Q5)*(SI
     &N(Q6)*COS(Q5)*U1-U5-(SIN(Q4)*COS(Q6)+SIN(Q5)*SIN(Q6)*COS(Q4))*U3-(
     &COS(Q4)*COS(Q6)-SIN(Q4)*SIN(Q5)*SIN(Q6))*U2))) - SIN(Q3)*COS(Q2)*(
     &IYA+MA*LAO**2+LA*MB*(LA+LBO*COS(Q5)*COS(Q6)))*U2
      AMOMY = (COS(Q1)*COS(Q3)-SIN(Q1)*SIN(Q2)*SIN(Q3))*(IYA+MA*LAO**2+L
     &A*MB*(LA+LBO*COS(Q5)*COS(Q6)))*U2 + (SIN(Q3)*COS(Q1)+SIN(Q1)*SIN(Q
     &2)*COS(Q3))*(IXA*U1-LA*LBO*MB*((SIN(Q6)*COS(Q4)+SIN(Q4)*SIN(Q5)*CO
     &S(Q6))*U2+(SIN(Q4)*SIN(Q6)-SIN(Q5)*COS(Q4)*COS(Q6))*U3)) + (SIN(Q5
     &)*(SIN(Q3)*COS(Q1)+SIN(Q1)*SIN(Q2)*COS(Q3))-SIN(Q1)*COS(Q2)*COS(Q4
     &)*COS(Q5)-SIN(Q4)*COS(Q5)*(COS(Q1)*COS(Q3)-SIN(Q1)*SIN(Q2)*SIN(Q3)
     &))*(IZB*U6+IZB*SIN(Q5)*U1+IZB*COS(Q4)*COS(Q5)*U3-IZB*SIN(Q4)*COS(Q
     &5)*U2-LBO*MB*(LBO+LA*COS(Q5)*COS(Q6))*(SIN(Q4)*COS(Q5)*U2-U6-SIN(Q
     &5)*U1-COS(Q4)*COS(Q5)*U3)) - LAO*MA*(PZp*COS(Q2)*COS(Q3)-PXp*(SIN(
     &Q1)*SIN(Q3)-SIN(Q2)*COS(Q1)*COS(Q3))) - MB*(LA*PZp*COS(Q2)*COS(Q3)
     &+LBO*PZp*(COS(Q2)*COS(Q3)*COS(Q5)*COS(Q6)+SIN(Q2)*(SIN(Q4)*SIN(Q6)
     &-SIN(Q5)*COS(Q4)*COS(Q6))-SIN(Q3)*COS(Q2)*(SIN(Q6)*COS(Q4)+SIN(Q4)
     &*SIN(Q5)*COS(Q6)))-LA*PXp*(SIN(Q1)*SIN(Q3)-SIN(Q2)*COS(Q1)*COS(Q3)
     &)-LBO*PXp*(COS(Q1)*COS(Q2)*(SIN(Q4)*SIN(Q6)-SIN(Q5)*COS(Q4)*COS(Q6
     &))+COS(Q5)*COS(Q6)*(SIN(Q1)*SIN(Q3)-SIN(Q2)*COS(Q1)*COS(Q3))+(SIN(
     &Q1)*COS(Q3)+SIN(Q2)*SIN(Q3)*COS(Q1))*(SIN(Q6)*COS(Q4)+SIN(Q4)*SIN(
     &Q5)*COS(Q6)))) - SIN(Q1)*COS(Q2)*(IZA+MA*LAO**2+LA*MB*(LA+LBO*COS(
     &Q5)*COS(Q6)))*U3 - (SIN(Q1)*COS(Q2)*(SIN(Q4)*COS(Q6)+SIN(Q5)*SIN(Q
     &6)*COS(Q4))+SIN(Q6)*COS(Q5)*(SIN(Q3)*COS(Q1)+SIN(Q1)*SIN(Q2)*COS(Q
     &3))-(COS(Q1)*COS(Q3)-SIN(Q1)*SIN(Q2)*SIN(Q3))*(COS(Q4)*COS(Q6)-SIN
     &(Q4)*SIN(Q5)*SIN(Q6)))*(IYB*U5+IYB*(SIN(Q4)*COS(Q6)+SIN(Q5)*SIN(Q6
     &)*COS(Q4))*U3+IYB*(COS(Q4)*COS(Q6)-SIN(Q4)*SIN(Q5)*SIN(Q6))*U2-IYB
     &*SIN(Q6)*COS(Q5)*U1-LBO*MB*(LBO+LA*COS(Q5)*COS(Q6))*(SIN(Q6)*COS(Q
     &5)*U1-U5-(SIN(Q4)*COS(Q6)+SIN(Q5)*SIN(Q6)*COS(Q4))*U3-(COS(Q4)*COS
     &(Q6)-SIN(Q4)*SIN(Q5)*SIN(Q6))*U2)) - (SIN(Q1)*COS(Q2)*(SIN(Q4)*SIN
     &(Q6)-SIN(Q5)*COS(Q4)*COS(Q6))-COS(Q5)*COS(Q6)*(SIN(Q3)*COS(Q1)+SIN
     &(Q1)*SIN(Q2)*COS(Q3))-(SIN(Q6)*COS(Q4)+SIN(Q4)*SIN(Q5)*COS(Q6))*(C
     &OS(Q1)*COS(Q3)-SIN(Q1)*SIN(Q2)*SIN(Q3)))*(IXB*U4+IXB*COS(Q5)*COS(Q
     &6)*U1+IXB*(SIN(Q6)*COS(Q4)+SIN(Q4)*SIN(Q5)*COS(Q6))*U2+IXB*(SIN(Q4
     &)*SIN(Q6)-SIN(Q5)*COS(Q4)*COS(Q6))*U3+LA*LBO*MB*(SIN(Q5)*(SIN(Q4)*
     &COS(Q5)*U2-U6-SIN(Q5)*U1-COS(Q4)*COS(Q5)*U3)-SIN(Q6)*COS(Q5)*(SIN(
     &Q6)*COS(Q5)*U1-U5-(SIN(Q4)*COS(Q6)+SIN(Q5)*SIN(Q6)*COS(Q4))*U3-(CO
     &S(Q4)*COS(Q6)-SIN(Q4)*SIN(Q5)*SIN(Q6))*U2)))
      AMOMZ = LAO*MA*(PYp*COS(Q2)*COS(Q3)-PXp*(SIN(Q3)*COS(Q1)+SIN(Q1)*S
     &IN(Q2)*COS(Q3))) + COS(Q1)*COS(Q2)*(IZA+MA*LAO**2+LA*MB*(LA+LBO*CO
     &S(Q5)*COS(Q6)))*U3 + (SIN(Q1)*COS(Q3)+SIN(Q2)*SIN(Q3)*COS(Q1))*(IY
     &A+MA*LAO**2+LA*MB*(LA+LBO*COS(Q5)*COS(Q6)))*U2 + (SIN(Q1)*SIN(Q3)-
     &SIN(Q2)*COS(Q1)*COS(Q3))*(IXA*U1-LA*LBO*MB*((SIN(Q6)*COS(Q4)+SIN(Q
     &4)*SIN(Q5)*COS(Q6))*U2+(SIN(Q4)*SIN(Q6)-SIN(Q5)*COS(Q4)*COS(Q6))*U
     &3)) + (COS(Q1)*COS(Q2)*COS(Q4)*COS(Q5)+SIN(Q5)*(SIN(Q1)*SIN(Q3)-SI
     &N(Q2)*COS(Q1)*COS(Q3))-SIN(Q4)*COS(Q5)*(SIN(Q1)*COS(Q3)+SIN(Q2)*SI
     &N(Q3)*COS(Q1)))*(IZB*U6+IZB*SIN(Q5)*U1+IZB*COS(Q4)*COS(Q5)*U3-IZB*
     &SIN(Q4)*COS(Q5)*U2-LBO*MB*(LBO+LA*COS(Q5)*COS(Q6))*(SIN(Q4)*COS(Q5
     &)*U2-U6-SIN(Q5)*U1-COS(Q4)*COS(Q5)*U3)) + (COS(Q1)*COS(Q2)*(SIN(Q4
     &)*SIN(Q6)-SIN(Q5)*COS(Q4)*COS(Q6))+COS(Q5)*COS(Q6)*(SIN(Q1)*SIN(Q3
     &)-SIN(Q2)*COS(Q1)*COS(Q3))+(SIN(Q1)*COS(Q3)+SIN(Q2)*SIN(Q3)*COS(Q1
     &))*(SIN(Q6)*COS(Q4)+SIN(Q4)*SIN(Q5)*COS(Q6)))*(IXB*U4+IXB*COS(Q5)*
     &COS(Q6)*U1+IXB*(SIN(Q6)*COS(Q4)+SIN(Q4)*SIN(Q5)*COS(Q6))*U2+IXB*(S
     &IN(Q4)*SIN(Q6)-SIN(Q5)*COS(Q4)*COS(Q6))*U3+LA*LBO*MB*(SIN(Q5)*(SIN
     &(Q4)*COS(Q5)*U2-U6-SIN(Q5)*U1-COS(Q4)*COS(Q5)*U3)-SIN(Q6)*COS(Q5)*
     &(SIN(Q6)*COS(Q5)*U1-U5-(SIN(Q4)*COS(Q6)+SIN(Q5)*SIN(Q6)*COS(Q4))*U
     &3-(COS(Q4)*COS(Q6)-SIN(Q4)*SIN(Q5)*SIN(Q6))*U2))) - MB*(LA*PXp*(SI
     &N(Q3)*COS(Q1)+SIN(Q1)*SIN(Q2)*COS(Q3))-LA*PYp*COS(Q2)*COS(Q3)-LBO*
     &PYp*(COS(Q2)*COS(Q3)*COS(Q5)*COS(Q6)+SIN(Q2)*(SIN(Q4)*SIN(Q6)-SIN(
     &Q5)*COS(Q4)*COS(Q6))-SIN(Q3)*COS(Q2)*(SIN(Q6)*COS(Q4)+SIN(Q4)*SIN(
     &Q5)*COS(Q6)))-LBO*PXp*(SIN(Q1)*COS(Q2)*(SIN(Q4)*SIN(Q6)-SIN(Q5)*CO
     &S(Q4)*COS(Q6))-COS(Q5)*COS(Q6)*(SIN(Q3)*COS(Q1)+SIN(Q1)*SIN(Q2)*CO
     &S(Q3))-(SIN(Q6)*COS(Q4)+SIN(Q4)*SIN(Q5)*COS(Q6))*(COS(Q1)*COS(Q3)-
     &SIN(Q1)*SIN(Q2)*SIN(Q3)))) - (SIN(Q6)*COS(Q5)*(SIN(Q1)*SIN(Q3)-SIN
     &(Q2)*COS(Q1)*COS(Q3))-COS(Q1)*COS(Q2)*(SIN(Q4)*COS(Q6)+SIN(Q5)*SIN
     &(Q6)*COS(Q4))-(SIN(Q1)*COS(Q3)+SIN(Q2)*SIN(Q3)*COS(Q1))*(COS(Q4)*C
     &OS(Q6)-SIN(Q4)*SIN(Q5)*SIN(Q6)))*(IYB*U5+IYB*(SIN(Q4)*COS(Q6)+SIN(
     &Q5)*SIN(Q6)*COS(Q4))*U3+IYB*(COS(Q4)*COS(Q6)-SIN(Q4)*SIN(Q5)*SIN(Q
     &6))*U2-IYB*SIN(Q6)*COS(Q5)*U1-LBO*MB*(LBO+LA*COS(Q5)*COS(Q6))*(SIN
     &(Q6)*COS(Q5)*U1-U5-(SIN(Q4)*COS(Q6)+SIN(Q5)*SIN(Q6)*COS(Q4))*U3-(C
     &OS(Q4)*COS(Q6)-SIN(Q4)*SIN(Q5)*SIN(Q6))*U2))

C**   Write output to screen and to output file(s)
      WRITE(*, 6020) GUESS(1),Q1,GUESS(2),Q2,GUESS(3),Q3,GUESS(4),Q4,GUE
     &SS(5),Q5,GUESS(6),Q6,GUESS(7),U1,GUESS(8),U2,GUESS(9),U3,GUESS(10)
     &,U4,GUESS(11),U5,GUESS(12),U6,T,P1X,P1Y,P1Z,P2X,P2Y,P2Z,P3X,P3Y,P3
     &Z,T,KE,PE,TE,AMOMX,AMOMY,AMOMZ
      WRITE(21,6020) GUESS(1),Q1,GUESS(2),Q2,GUESS(3),Q3,GUESS(4),Q4,GUE
     &SS(5),Q5,GUESS(6),Q6,GUESS(7),U1,GUESS(8),U2,GUESS(9),U3,GUESS(10)
     &,U4,GUESS(11),U5,GUESS(12),U6,T,P1X,P1Y,P1Z,P2X,P2Y,P2Z,P3X,P3Y,P3
     &Z,T,KE,PE,TE,AMOMX,AMOMY,AMOMZ
      WRITE( *,6022) (ILOOP,RHSGUESS(ILOOP),RHS(ILOOP), ILOOP = 1,12)
      WRITE(21,6022) (ILOOP,RHSGUESS(ILOOP),RHS(ILOOP), ILOOP = 1,12)

6020  FORMAT(//1X,'QUANTITY',23X,'GUESS',16X,'OUTPUT RESULT',/4X,'Q1',21
     &X,1PE14.6E3,11X,1PE14.6E3,/4X,'Q2',21X,1PE14.6E3,11X,1PE14.6E3,/4X
     &,'Q3',21X,1PE14.6E3,11X,1PE14.6E3,/4X,'Q4',21X,1PE14.6E3,11X,1PE14
     &.6E3,/4X,'Q5',21X,1PE14.6E3,11X,1PE14.6E3,/4X,'Q6',21X,1PE14.6E3,1
     &1X,1PE14.6E3,/4X,'U1',21X,1PE14.6E3,11X,1PE14.6E3,/4X,'U2',21X,1PE
     &14.6E3,11X,1PE14.6E3,/4X,'U3',21X,1PE14.6E3,11X,1PE14.6E3,/4X,'U4'
     &,21X,1PE14.6E3,11X,1PE14.6E3,/4X,'U5',21X,1PE14.6E3,11X,1PE14.6E3,
     &/4X,'U6',21X,1PE14.6E3,11X,1PE14.6E3,/4X,'T',47X,1PE14.6E3,/4X,'P1
     &X',45X,1PE14.6E3,/4X,'P1Y',45X,1PE14.6E3,/4X,'P1Z',45X,1PE14.6E3,/
     &4X,'P2X',45X,1PE14.6E3,/4X,'P2Y',45X,1PE14.6E3,/4X,'P2Z',45X,1PE14
     &.6E3,/4X,'P3X',45X,1PE14.6E3,/4X,'P3Y',45X,1PE14.6E3,/4X,'P3Z',45X
     &,1PE14.6E3,/4X,'T',47X,1PE14.6E3,/4X,'KE',46X,1PE14.6E3,/4X,'PE',4
     &6X,1PE14.6E3,/4X,'TE',46X,1PE14.6E3,/4X,'AMOMX',43X,1PE14.6E3,/4X,
     &'AMOMY',43X,1PE14.6E3,/4X,'AMOMZ',43X,1PE14.6E3)
6022  FORMAT(//1X,'FUNCTION',17X,'EVALUATED AT GUESS',6X,'EVALUATED AT R
     &ESULT',99(/1X,I4,22X,1PE14.6E3,11X,1PE14.6E3))

      RETURN
      END


C*****************************************************************************
C**                                                                         **
C** PURPOSE  Solves a set of first order ordinary differential equations    **
C**          of the form dy(i)/dt = F(t,y(1), ..., y(numeqns) (i = 1,       **
C**          ..., numeqns)                                                  **
C**                                                                         **
C** INPUT                                                                   **
C**    eqns: Subroutine that evaluates dy(i)/dt (i = 1, ..., numeqns), the  **
C**          first derivatives of y(1), ..., y(numeqns) with respect to t   **
C**                                                                         **
C** numeqns: The number of differential equations to be solved              **
C**                                                                         **
C**       y: One-dimensional array whose elements are y(1), ..., y(numeqns) **
C**                                                                         **
C**       t: Independent variable                                           **
C**                                                                         **
C** integstp: Maximum integration stepsize                                  **
C**                                                                         **
C**  abserr: Allowable absolute error in y(i)  (i=1, ..., numeqns)          **
C**                                                                         **
C**  relerr: Allowable relative error in y(i)  (i=1, ..., numeqns)          **
C**                                                                         **
C**     com: When com = 2, the Kutta-Merson algorithm (L. Fox, Numerical    **
C**          Solutions of Ordinary and Partial Differential Equations,      **
C**          Palo Alto: Addison-Wesley, 1962, pp. 24-25) is employed to     **
C**          perform the numerical solution of the differential equations.  **
C**          Accordingly, dy(i)/dt (i = 1, ..., numeqns) are evaluated at   **
C**          every integration boundary, including those at Tinitial,       **
C**          Tfinal, and ones created when integstp is halved to satisfy    **
C**          the requirements imposed by abserr and relerr.  Integration    **
C**          is self-starting at each boundary, and the occurrence, at      **
C**          boundaries, of discontinuities in derivatives does not lead    **
C**          to failure of the integration process.                         **
C**                                                                         **
C**          When com = 1, a modified form of the Kutta-Merson algorithm    **
C**          is employed.  It is nearly 20% faster than the one used when   **
C**          com = 2 because no recalculation of derivatives at inte-       **
C**          gration boundaries between Tinitial and Tfinal takes place.    **
C**          Integration is self-starting at Tinitial and Tfinal only.      **
C**          Integration may fail if any of dy(i)/dt (i = 1, ..., numeqns)  **
C**          is discontinuous between Tinitial and Tfinal.                  **
C**                                                                         **
C**          When com = 0, the function eqns is called and dy(i)/dt         **
C**          (i = 1, ..., numeqns) are evaluated, but no integration        **
C**          is performed.                                                  **
C**                                                                         **
C** OUTPUT                                                                  **
C**          The value of t+integstp is returned in t, and the values of    **
C**          y(i) at t+integstp are returned in y.                          **
C**                                                                         **
C** SOURCE                                                                  **
C**          Copyright 1995 by Paul C. Mitiguy, Thomas R. Kane, David A.    **
C**          Levinson, and David B. Schaechter.  Permission is granted      **
C**          to copy, modify, and distribute this subroutine, provided      **
C**          that this copyright notice appear.                             **
C**                                                                         **
C*****************************************************************************
      SUBROUTINE KUTTA (EQNS,NUMY,Y,T,INTEGSTP,ABSERR,RELERR,COM,*)
      EXTERNAL         EQNS
      INTEGER          NUMY, COM, NUMCUTS, I
      LOGICAL          STEPDBL, ENTRY
      DOUBLE PRECISION Y(NUMY), F0, F1, F2, Y1, Y2
      DOUBLE PRECISION T, INTEGSTP, ABSERR, RELERR, ERROR, TEST
      DOUBLE PRECISION TFINAL, TT, HC, H, H2, H3, H6, H8
      COMMON/CKUTTA/   F0(100),F1(100),F2(100),Y1(100),Y2(100)
      DATA             HC, NUMCUTS / 0.0D0, 20 /

C**   If COM=0, call EQNS subroutine and return.
      IF( COM .EQ. 0) THEN
        CALL EQNS(T, Y, F0, 1)
        RETURN
      ENDIF

C**   Check for initial entry and adjust current value of stepsize.
      IF(NUMY .EQ. 0) THEN
        HC = INTEGSTP
        RETURN
      ENDIF
      IF(INTEGSTP .EQ. 0) RETURN 1
      IF(HC*INTEGSTP .LT. 0) HC = -HC
      IF(HC .EQ. 0)          HC = INTEGSTP

C**   Set local variables
      H = HC
      TT = T + H
      TFINAL = T + INTEGSTP
      T  = TFINAL
      ENTRY = .TRUE.

C**   Check round-off problems.
100   IF( TT+H .EQ. TT ) THEN
        T = TT
        WRITE(*,2010) H, T
        CALL EQNS(T, Y, F0, 0)
        RETURN 1
      ENDIF
C**   Main Kutta-Merson step
      H2 = H * 0.5D0
      H3 = H / 3.0D0
      H6 = H / 6.0D0
      H8 = H * 0.125D0
      IF( COM .EQ. 2 .OR. ENTRY )  CALL EQNS(TT-H, Y, F0, 1)
      ENTRY = .FALSE.
      DO 110  I=1,NUMY
110     Y1(I) = Y(I) + H3*F0(I)
      CALL EQNS(TT-2.0*H3, Y1, F1, 0)
      DO 120  I=1,NUMY
120     Y1(I) = Y(I) + H6*(F0(I) + F1(I))
      CALL EQNS(TT-2.0*H3, Y1, F1, 0)
      DO 130  I=1,NUMY
130     Y1(I) = Y(I) + H8*(F0(I) + 3.0D0*F1(I) )
      CALL EQNS(TT-H2,     Y1, F2, 0)
      DO 140  I=1,NUMY
140     Y1(I) = Y(I) + H2*(F0(I) - 3.0D0*F1(I)+ 4.0D0*F2(I) )
      CALL EQNS(TT,        Y1, F1, 0)
      DO 150  I=1,NUMY
150     Y2(I) = Y(I) + H6*(F0(I) +  4.0D0*F2(I) + F1(I) )
C**   Assume that step needs to be doubled.  Check error criterion
      STEPDBL = .TRUE.
      DO 160 I=1,NUMY
        ERROR = DABS(Y1(I) - Y2(I)) * 0.2D0
        TEST  = DABS(Y1(I)) * RELERR
        IF(ERROR .GE. TEST .AND. ERROR .GE. ABSERR) THEN
          HC = H2
          H  = HC
          TT = TT - H2
          NUMCUTS = NUMCUTS - 1
          IF(NUMCUTS .GE. 0) GO TO 100
          T = TT - H
          WRITE(*,2000) T
          CALL EQNS(T, Y, F0, 0)
          RETURN 1
        ENDIF
      IF(STEPDBL .AND. 64.0D0*ERROR .GT. TEST
     &           .AND. 64.0D0*ERROR .GT. ABSERR) STEPDBL=.FALSE.
160   CONTINUE
      DO 170  I = 1,NUMY
170     Y(I) = Y2(I)
C**   Double the STEPSIZE, maybe.
      IF( STEPDBL .AND. DABS(H+H) .LE. DABS(INTEGSTP) .AND.
     &     DABS(TT+H+H) .LE. DABS(TFINAL) )  THEN
        HC = H + H
        H  = HC
        NUMCUTS = NUMCUTS + 1
      ENDIF
      IF( TT .EQ. TFINAL ) THEN
        CALL EQNS(TFINAL, Y, F0, 2)
        RETURN
      ENDIF
      TT = TT + H
      IF( (H .GT. 0 .AND. TT .GT. TFINAL-0.1D0*H) .OR.
     &    (H .LT. 0 .AND. TT .LT. TFINAL-0.1D0*H)  )  THEN
        H  = TFINAL - (TT-H)
        TT = TFINAL
      ENDIF
      IF( COM .EQ. 1 ) THEN
        DO 180  I = 1,NUMY
180       F0(I) = F1(I)
      ENDIF
      GOTO 100 

2000  FORMAT(/1X,'THE STEPSIZE HAS BEEN HALVED TOO MANY TIMES; T = ',
     &1PD12.4,/1X,'ERROR: NUMERICAL INTEGRATION FAILED TO CONVERGE.',//)
2010  FORMAT(/1X,'THE STEPSIZE OF ',1PD22.14,' IS TOO SMALL RELATIVE ', 
     &'TO THE TERMINAL TIME OF',/1PD22.14,'.  INTEGRATION HALTED BECA',
     &'USE OF NUMERICAL ROUND-OFF.',/,'THE STEPSIZE MAY HAVE BEEN CUT ',
     &'TOO MANY TIMES.'//)
      END



C**************************************************************************** 
C**                                                                        ** 
C** PURPOSE  The matrix equation a x = b is solved for x, where a is an    ** 
C**          n by n matrix, and x and b are n by 1 matrices.               ** 
C**                                                                        ** 
C** INPUT                                                                  **
C**       N: n                                                             ** 
C**                                                                        ** 
C**       A: an N by N double precision array whose elements are those     **      
C**          of the matrix a                                               ** 
C**                                                                        ** 
C**       B: an N by 1 double precision array whose elements are those     **      
C**          of the matrix b                                               ** 
C**                                                                        ** 
C** OUTPUT                                                                 ** 
C**       X: an N by 1 double precision array whose elements are those     **
C**          of the matrix x                                               ** 
C**                                                                        ** 
C**************************************************************************** 
        SUBROUTINE SOLVE(N, A, B, X)
        IMPLICIT DOUBLE PRECISION (A - Z)
        INTEGER N,IPS(100),I,J,K,IP,KP,KP1,NM1,IDXPIV,IP1,IM1,NP1,IBACK
        DIMENSION A(N,N),SCALES(100),B(N),X(N)

C*************** Beginning of LU decomposition of A ********************
        ZERO = 0.0D0
        DO 5 I=1,N
        IPS(I) = I
        ROWNRM = 0.0D0
        DO 20 J=1,N
        ROWNRM = DMAX1(ROWNRM,DABS(A(I,J)))
   20   CONTINUE
        IF(ROWNRM.EQ.ZERO) GOTO 500
        SCALES(I) = 1.0D0 / ROWNRM
    5   CONTINUE
        NM1 = N-1
        DO 17 K=1,NM1
        BIG = 0.0D0
        DO 11 I=K,N
        IP = IPS(I)
        SIZE = DABS(A(IP,K))*SCALES(IP)
        IF(SIZE .LE. BIG) GO TO 11
        BIG = SIZE
        IDXPIV = I
   11   CONTINUE
        IF(BIG .EQ. ZERO) GOTO 520
        IF(IDXPIV .EQ. K) GO TO 15
        J = IPS(K)
        IPS(K) = IPS(IDXPIV)
        IPS(IDXPIV) = J
   15   KP = IPS(K)
        PIVOT = A(KP,K)
        KP1 = K+1
        DO 16 I=KP1,N
        IP = IPS(I)
        EM = A(IP,K)/PIVOT
        A(IP,K) = EM
        DO 16 J = KP1,N
        A(IP,J) = A(IP,J) - EM*A(KP,J)
   16   CONTINUE
   17   CONTINUE
        IF(A(IPS(N),N) .EQ. ZERO) GOTO 520

C**     Note: The LU decomposition of A is returned in A
C***************** Beginning of back substitution **********************
        NP1 = N+1
        X(1) = B(IPS(1))
        DO 2 I=2,N
        IP = IPS(I)
        IM1 = I-1
        SUM = 0.0D0
        DO 1 J=1,IM1
        SUM = SUM + A(IP,J)*X(J)
    1   CONTINUE
        X(I) = B(IP) - SUM
    2   CONTINUE
        X(N) = X(N)/A(IPS(N),N)
        DO 4 IBACK=2,N
        I = NP1-IBACK
        IP = IPS(I)
        IP1 = I+1
        SUM = 0.0D0
        DO 3 J=IP1,N
        SUM = SUM + A(IP,J)*X(J)
    3   CONTINUE
    4   X(I) = (X(I)-SUM)/A(IP,I)
        RETURN

  500  WRITE(*,600) I
       STOP
  520  WRITE(*,620)
       STOP
  600  FORMAT(/1X,'ALL ELEMENTS IN ROW ',I3,'   OF COEF ARE ZEROS'/)
  620  FORMAT(/1X,'A PIVOT ELEMENT ENCOUNTERED IN THE DECOMPOSITION',
     & ' OF COEF IS ZERO',/15X,'COEFFICIENT MATRIX IS SINGULAR')
        END



