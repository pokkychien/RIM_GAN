C                  
      SUBROUTINE CONSX(N,M,K,ITMAX,ALPHA,BETA,GAMMA,DELTA,X,R,F
     1,IT,IEV2,NO,G,H,XC,IPRINT)
      implicit real*8 (a-h,o-z)
C
C      COORDINATES SPECIAL PURPOSE SUBROUTINES
C
C      ARGUMENT LIST
C
C      IT    = ITERATION INDEX
C      IEV1  = INDEX OF POENT WITH MINIMUM FUNCTION VALUE.
C      IEV2  = INDEX OF POINT WITH MAXIMUM FUNCTION VALUE
C      I     = POINT INDEX.
C      KODE  = CONTROL KEY USED TO DETERMINE IF IMPLICIT CONSTRAINTS
C              ARE PROVIDED
C
C      ALL OTHERS PREVIOUSLY DEFINED IN MAIN LINE.
C
      DIMENSION X(K,M),R(K,N),F(K),G(M),H(M),XC(N)
      INTEGER GAMMA
C
C----  IQ =0 PRINT OUT RESULTS FOR EVERY 50 ITERATIONS
	IQ=0
      IT=1
      KODE=0
      IF(M-N)20,20,10
10    KODE=1
20     CONTINUE
      DO 40 II=2,K
      DO 30 J=1,N
30    X(II,J)=0.0
40    CONTINUE
C
C      CALCULATE COMPLEX POINTS AND CHECK AGAINST CONSTRAINTS
C
      DO 65 II=2,K
      DO 50 J=1,N
      I=II
      CALL CONST(N,M,K,X,G,H,I)
      X(II,J)=G(J)+R(II,J)*(H(J)-G(J))
50    CONTINUE
      K1=II
      CALL  CHECK(N,M,K,X,G,H,I,KODE,XC,DELTA,K1)
       IF(II-2)51,51,55
51    IF(IPRINT)52,65,52
52     WRITE(NO,18)
18    FORMAT(//,2X,30HCOORDINATES OF INITIAL COMPLEX)
      IO=1
       WRITE(NO,19)(IO,J,X(IO,J),J=1,N)
19    FORMAT(/,3(1X,2HX(,I2,1H,,I2,3H) =,F9.4))
55    IF(IPRINT)56,65,56
56     WRITE(NO,19)(II,J,X(II,J),J=1,N)
65    CONTINUE
      K1=K
      DO 70 I=1,K
	WRITE(6,*)'I=',I
      CALL FUNC(N,M,K,X,F,I)
70    CONTINUE
        ITR=MOD(IT,50)
        IF(IQ.EQ.1) ITR=0
      KOUNT=1
      IA=0
C
C      FIND POINT WITH LOWEST FUNCTION VALUE
C
      IF(IPRINT)72,80,72
72    IF(ITR.EQ.0) WRITE(NO,21)
21    FORMAT(/,2X,22HVALUES OF THE FUNCTION )
      IF(ITR.EQ.0) WRITE(NO,22)(J,F(J),J=1,K)
22    FORMAT(/,3(2X,2HF(,I2,4H) = ,E10.3))
80    IEV1=1
        ITR=MOD(IT,50)
        IF(IQ.EQ.1) ITR=0
      DO 100  ICM=2,K
      IF(F(IEV1)-F(ICM))100,100,90
90    IEV1=ICM
100   CONTINUE
C
C      FIND POINT WITH HIGHEST FUNCTION VALUE
C
      IEV2=1
      DO 120 ICM=2,K
      IF(F(IEV2)-F(ICM))110,110,120
110   IEV2=ICM
120   CONTINUE
C
C      CHECK CONVERGENCE CRITERIA
C
      IF (F(IEV2)-(F(IEV1)+BETA))140,130,130
130   KOUNT=1
      GO TO 150
140   KOUNT=KOUNT+1
      IF(KOUNT-GAMMA)150,240,240
C
C      REPLACE POINT WITH LOWEST FUNCTION VALUE
C
150   CALL  CENTR(N,M,K,IEV1,I,XC,X,K1)
      DO 160 JJ=1,N
160   X(IEV1,JJ)=(1.0+ALPHA)*(XC(JJ))-ALPHA*(X(IEV1,JJ))
      I=IEV1
      CALL CHECK(N,M,K,X,G,H,I,KODE,XC,DELTA,K1)
      CALL FUNC(N,M,K,X,F,I)
C
C      REPLACE NEW POINT IF IT REPEATS AS LOWEST FUNCTION VALUE
C
      ICNT=0
170   ICNT=ICNT+1
      IF (ICNT.LT.10) GO TO 190
      ICNT=0
      IF (IPRINT.NE.1) GO TO 183
      IF(ITR.EQ.0) WRITE(NO,180)
180   FORMAT(/,2X,'CENTROID HAS MINIMUM VALUE')
      IF(ITR.EQ.0)WRITE(NO,1001)IEV1,IEV2
 1001	FORMAT(' IEV1,IEV2=',2I10)
      IF(ITR.EQ.0)WRITE(NO,1002)IEV1,IEV2,F(IEV1),F(IEV2)
 1002	FORMAT('F(',I10,'),F(',I10,')=',2F12.4)
      IF(ITR.EQ.0)WRITE(NO,1003)(X(IEV1,JJ),JJ=1,N)
 1003	FORMAT('X(IEV1)',5F10.4)
      IF(ITR.EQ.0)WRITE(NO,1004)(X(IEV2,JJ),JJ=1,N)
 1004	FORMAT('X(IEV2)=',5F10.4)
      IF(ITR.EQ.0)WRITE(NO,1005)(XC(JJ),JJ=1,N)
 1005	FORMAT('  XC=',5F10.4)
183   DO 185 JJ=1,N
185      XC(JJ)=X(IEV2,JJ)
190   IMIN=1
      DO 195 II=2,K
195      IF (F(II).LT.F(IMIN)) IMIN=II
      IF (IMIN.NE.IEV1) GO TO 220
200   DO 210 JJ=1,N
      X(IEV1,JJ)=(X(IEV1,JJ)+XC(JJ))/2.0
210   CONTINUE
      I=IEV1
      CALL CHECK(N,M,K,X,G,H,I,KODE,XC,DELTA,K1)
      CALL FUNC(N,M,K,X,F,I)
      GO TO 170
220   CONTINUE
      IF(F(IEV2).LT.F(IEV1))IEV2=IEV1
      IF(IPRINT)230,228,230
230   IF(ITR.EQ.0) WRITE(NO,23)IT
23    FORMAT(//,2X,17HITERATION NUMBER  ,I5)
      IF(ITR.EQ.0) WRITE(NO,24)
24    FORMAT(/,2X,30HCOORDINATES OF CORRECTED POINT)
      IF(ITR.EQ.0) WRITE(NO,19)(IEV1,JC,X(IEV1,JC),JC=1,N)
      IF(ITR.EQ.0) WRITE(NO,21)
      IF(ITR.EQ.0) WRITE(NO,22)(I,F(I),I=1,K)
      IF(ITR.EQ.0) WRITE(NO,25)
25    FORMAT(/,2X,27HCOORDINATES OF THE CENTROID)
      IF(ITR.EQ.0) WRITE(NO,26)(JC,XC(JC),JC=1,N)
26    FORMAT(/,3(1X,2HX(,I2,6H,C) = ,F9.4,1X))
228   IT=IT+1
      WRITE(40,*) 'IT=',IT
      IF(IT-ITMAX)80,80,240
240   RETURN
      END
C
C************************************************************
C
      SUBROUTINE CHECK(N,M,K,X,G,H,I,KODE,XC,DELTA,K1)
C
C      ARGUMENT LIST
C
C      ALL ARGUMENTS  DEFINED IN MAIN LINE AND CONSX
C
      implicit real*8 (a-h,o-z)
      DIMENSION X(K,M),G(M),H(M),XC(N)
C
10    KT=0
      CALL CONST(N,M,K,X,G,H,I)
C
C      CHECK AGAINST EXPLICIT CONSTRAINTS
C
      DO 50 J=1,N
      IF(X(I,J)-G(J))20,20,30
20    X(I,J)=G(J)+DELTA
      GO  TO 50
30    IF(H(J)-X(I,J))40,40,50
40    X(I,J)=H(J)-DELTA
50    CONTINUE
      DO 55 J=1,N
55       IF (G(J).EQ.H(J)) X(I,J)=G(J)
C
      IF(KODE)110,110,60
C
C      CHECK AGAINST THE  IMPLICIT  CONSTRAINTS
C
60    NN=N+1
      DO 100 J=NN,M
      CALL CONST(N,M,K,X,G,H,I)
      IF(X(I,J)-G(J))80,70,70
70    IF(H(J)-X(I,J))80,100,100
80    IEV1=I
      KT=1
      CALL  CENTR(N,M,K,IEV1,I,XC,X,K1)
      DO 90  JJ=1,N
      X(I,JJ)=(X(I,JJ)+XC(JJ))/2.0
90    CONTINUE
100   CONTINUE
      IF(KT)110,110,10
110   RETURN
      END
C
C***********************************************************
C
      SUBROUTINE CENTR(N,M,K,IEV1,I,XC,X,K1)
C
      implicit real*8 (a-h,o-z)
      DIMENSION X(K,M),XC(N)
C
      DO 20 J=1,N
      XC(J)=0.0
      DO 10 IL=1,K1
10    XC(J)=XC(J)+X(IL,J)
      RK=K1
20    XC(J)=(XC(J)-X(IEV1,J))/(RK-1.0)
      RETURN
      END
