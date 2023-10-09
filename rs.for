      SUBROUTINE RS(NM,N,A,W,MATZ,Z,FV1,FV2,IERR)
C
      INTEGER N,NM,IERR,MATZ
      REAL A(NM,N),W(N),Z(NM,N),FV1(N),FV2(N)
      IF (N .LE. NM) GO TO 10
      IERR = 10 * N
      GO TO 50
C
   10 IF (MATZ .NE. 0) GO TO 20
C     ********** FIND EIGENVALUES ONLY **********
      CALL  TRED1(NM,N,A,W,FV1,FV2)
      CALL  TQLRAT(N,W,FV2,IERR)
      GO TO 50
C     ********** FIND BOTH EIGENVALUES AND EIGENVECTORS **********
   20 CALL  TRED2(NM,N,A,W,FV1,Z)
      CALL  TQL2(NM,N,W,FV1,Z,IERR)
   50 RETURN
C     ********** LAST CARD OF RS **********
      END
C
      SUBROUTINE TQL1(N,D,E,IERR)
C
      INTEGER I,J,L,M,N,II,L1,MML,IERR
      REAL D(N),E(N)
      REAL B,C,F,G,H,P,R,S,MACHEP
C     REAL SQRT,ABS,SIGN
C
      MACHEP = 2.**(-23)
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
C
      F = 0.0
      B = 0.0
      E(N) = 0.0
C
      DO 290 L = 1, N
         J = 0
         H = MACHEP * (ABS(D(L)) + ABS(E(L)))
         IF (B .LT. H) B = H
C     ********** LOOK FOR SMALL SUB-DIAGONAL ELEMENT **********
         DO 110 M = L, N
            IF (ABS(E(M)) .LE. B) GO TO 120
C     ********** E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP **********
  110    CONTINUE
C
  120    IF (M .EQ. L) GO TO 210
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     ********** FORM SHIFT **********
         L1 = L + 1
         G = D(L)
         P = (D(L1) - G) / (2.0 * E(L))
         R = SQRT(P*P+1.0)
         D(L) = E(L) / (P + SIGN(R,P))
         H = G - D(L)
C
         DO 140 I = L1, N
  140    D(I) = D(I) - H
C
         F = F + H
C     ********** QL TRANSFORMATION **********
         P = D(M)
         C = 1.0
         S = 0.0
         MML = M - L
C     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********
         DO 200 II = 1, MML
            I = M - II
            G = C * E(I)
            H = C * P
            IF (ABS(P) .LT. ABS(E(I))) GO TO 150
            C = E(I) / P
            R = SQRT(C*C+1.0)
            E(I+1) = S * P * R
            S = C / R
            C = 1.0 / R
            GO TO 160
  150       C = P / E(I)
            R = SQRT(C*C+1.0)
            E(I+1) = S * E(I) * R
            S = 1.0 / R
            C = C * S
  160       P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
  200    CONTINUE
C
         E(L) = S * P
         D(L) = C * P
         IF (ABS(E(L)) .GT. B) GO TO 130
  210    P = D(L) + F
C     ********** ORDER EIGENVALUES **********
         IF (L .EQ. 1) GO TO 250
C     ********** FOR I=L STEP -1 UNTIL 2 DO -- **********
         DO 230 II = 2, L
            I = L + 2 - II
            IF (P .GE. D(I-1)) GO TO 270
            D(I) = D(I-1)
  230    CONTINUE
C
  250    I = 1
  270    D(I) = P
  290 CONTINUE
C
      GO TO 1001
C     ********** SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS **********
 1000 IERR = L
 1001 RETURN
C     ********** LAST CARD OF TQL1 **********
      END
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)
C
      INTEGER I,J,K,L,M,N,II,L1,NM,MML,IERR
      REAL D(N),E(N),Z(NM,N)
      REAL B,C,F,G,H,P,R,S,MACHEP
C     REAL SQRT,ABS,SIGN
      MACHEP = 2.**(-23)
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
C
      F = 0.0
      B = 0.0
      E(N) = 0.0
C
      DO 240 L = 1, N
         J = 0
         H = MACHEP * (ABS(D(L)) + ABS(E(L)))
         IF (B .LT. H) B = H
C     ********** LOOK FOR SMALL SUB-DIAGONAL ELEMENT **********
         DO 110 M = L, N
            IF (ABS(E(M)) .LE. B) GO TO 120
C     ********** E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP **********
  110    CONTINUE
C
  120    IF (M .EQ. L) GO TO 220
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     ********** FORM SHIFT **********
         L1 = L + 1
         G = D(L)
         P = (D(L1) - G) / (2.0 * E(L))
         R = SQRT(P*P+1.0)
         D(L) = E(L) / (P + SIGN(R,P))
         H = G - D(L)
C
         DO 140 I = L1, N
  140    D(I) = D(I) - H
C
         F = F + H
C     ********** QL TRANSFORMATION **********
         P = D(M)
         C = 1.0
         S = 0.0
         MML = M - L
C     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********
         DO 200 II = 1, MML
            I = M - II
            G = C * E(I)
            H = C * P
            IF (ABS(P) .LT. ABS(E(I))) GO TO 150
            C = E(I) / P
            R = SQRT(C*C+1.0)
            E(I+1) = S * P * R
            S = C / R
            C = 1.0 / R
            GO TO 160
  150       C = P / E(I)
            R = SQRT(C*C+1.0)
            E(I+1) = S * E(I) * R
            S = 1.0 / R
            C = C * S
  160       P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
C     ********** FORM VECTOR **********
            DO 180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
  180       CONTINUE
C
  200    CONTINUE
C
         E(L) = S * P
         D(L) = C * P
         IF (ABS(E(L)) .GT. B) GO TO 130
  220    D(L) = D(L) + F
  240 CONTINUE
C     ********** ORDER EIGENVALUES AND EIGENVECTORS **********
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
C
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
C
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
C
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
C
  300 CONTINUE
C
      GO TO 1001
C     ********** SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS **********
 1000 IERR = L
 1001 RETURN
C     ********** LAST CARD OF TQL2 **********
      END
C
      SUBROUTINE TQLRAT(N,D,E2,IERR)
C
      INTEGER I,J,L,M,N,II,L1,MML,IERR
      REAL D(N),E2(N)
      REAL B,C,F,G,H,P,R,S,MACHEP
C     REAL SQRT,ABS,SIGN
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQLRAT,
      MACHEP = 2.**(-23)
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E2(I-1) = E2(I)
C
      F = 0.0
      B = 0.0
      E2(N) = 0.0
C
      DO 290 L = 1, N
         J = 0
         H = MACHEP * (ABS(D(L)) + SQRT(E2(L)))
         IF (B .GT. H) GO TO 105
         B = H
         C = B * B
C     ********** LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT **********
  105    DO 110 M = L, N
            IF (E2(M) .LE. C) GO TO 120
C     ********** E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP **********
  110    CONTINUE
C
  120    IF (M .EQ. L) GO TO 210
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     ********** FORM SHIFT **********
         L1 = L + 1
         S = SQRT(E2(L))
         G = D(L)
         P = (D(L1) - G) / (2.0 * S)
         R = SQRT(P*P+1.0)
         D(L) = S / (P + SIGN(R,P))
         H = G - D(L)
C
         DO 140 I = L1, N
  140    D(I) = D(I) - H
C
         F = F + H
C     ********** RATIONAL QL TRANSFORMATION **********
         G = D(M)
         IF (G .EQ. 0.0) G = B
         H = G
         S = 0.0
         MML = M - L
C     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********
         DO 200 II = 1, MML
            I = M - II
            P = G * H
            R = P + E2(I)
            E2(I+1) = S * R
            S = E2(I) / R
            D(I+1) = H + S * (H + D(I))
            G = D(I) - E2(I) / G
            IF (G .EQ. 0.0) G = B
            H = G * P / R
  200    CONTINUE
C
         E2(L) = S * G
         D(L) = H
C     ********** GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST **********
         IF (H .EQ. 0.0) GO TO 210
         IF (ABS(E2(L)) .LE. ABS(C/H)) GO TO 210
         E2(L) = H * E2(L)
         IF (E2(L) .NE. 0.0) GO TO 130
  210    P = D(L) + F
C     ********** ORDER EIGENVALUES **********
         IF (L .EQ. 1) GO TO 250
C     ********** FOR I=L STEP -1 UNTIL 2 DO -- **********
         DO 230 II = 2, L
            I = L + 2 - II
            IF (P .GE. D(I-1)) GO TO 270
            D(I) = D(I-1)
  230    CONTINUE
C
  250    I = 1
  270    D(I) = P
  290 CONTINUE
C
      GO TO 1001
C     ********** SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS **********
 1000 IERR = L
 1001 RETURN
C     ********** LAST CARD OF TQLRAT **********
      END
C
      SUBROUTINE TRED1(NM,N,A,D,E,E2)
C
      INTEGER I,J,K,L,N,II,NM,JP1
      REAL A(NM,N),D(N),E(N),E2(N)
      REAL F,G,H,SCALE
C     REAL SQRT,ABS,SIGN
      DO 100 I = 1, N
  100 D(I) = A(I,I)
C     ********** FOR I=N STEP -1 UNTIL 1 DO -- **********
      DO 300 II = 1, N
         I = N + 1 - II
         L = I - 1
         H = 0.0
         SCALE = 0.0
         IF (L .LT. 1) GO TO 130
C     ********** SCALE ROW (ALGOL TOL THEN NOT NEEDED) **********
         DO 120 K = 1, L
  120    SCALE = SCALE + ABS(A(I,K))
C
         IF (SCALE .NE. 0.0) GO TO 140
  130    E(I) = 0.0
         E2(I) = 0.0
         GO TO 290
C
  140    DO 150 K = 1, L
            A(I,K) = A(I,K) / SCALE
            H = H + A(I,K) * A(I,K)
  150    CONTINUE
C
         E2(I) = SCALE * SCALE * H
         F = A(I,L)
         G = -SIGN(SQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         A(I,L) = F - G
         IF (L .EQ. 1) GO TO 270
         F = 0.0
C
         DO 240 J = 1, L
            G = 0.0
C     ********** FORM ELEMENT OF A*U **********
            DO 180 K = 1, J
  180       G = G + A(J,K) * A(I,K)
C
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
C
            DO 200 K = JP1, L
  200       G = G + A(K,J) * A(I,K)
C     ********** FORM ELEMENT OF P **********
  220       E(J) = G / H
            F = F + E(J) * A(I,J)
  240    CONTINUE
C
         H = F / (H + H)
C     ********** FORM REDUCED A **********
         DO 260 J = 1, L
            F = A(I,J)
            G = E(J) - H * F
            E(J) = G
C
            DO 260 K = 1, J
               A(J,K) = A(J,K) - F * E(K) - G * A(I,K)
  260    CONTINUE
C
  270    DO 280 K = 1, L
  280    A(I,K) = SCALE * A(I,K)
C
  290    H = D(I)
         D(I) = A(I,I)
         A(I,I) = H
  300 CONTINUE
C
      RETURN
C     ********** LAST CARD OF TRED1 **********
      END
C
      SUBROUTINE TRED2(NM,N,A,D,E,Z)
C
      INTEGER I,J,K,L,N,II,NM,JP1
      REAL A(NM,N),D(N),E(N),Z(NM,N)
      REAL F,G,H,HH,SCALE
C     REAL SQRT,ABS,SIGN
C
      DO 100 I = 1, N
C
         DO 100 J = 1, I
            Z(I,J) = A(I,J)
  100 CONTINUE
C
      IF (N .EQ. 1) GO TO 320
C     ********** FOR I=N STEP -1 UNTIL 2 DO -- **********
      DO 300 II = 2, N
         I = N + 2 - II
         L = I - 1
         H = 0.0
         SCALE = 0.0
         IF (L .LT. 2) GO TO 130
C     ********** SCALE ROW (ALGOL TOL THEN NOT NEEDED) **********
         DO 120 K = 1, L
  120    SCALE = SCALE + ABS(Z(I,K))
C
         IF (SCALE .NE. 0.0) GO TO 140
  130    E(I) = Z(I,L)
         GO TO 290
C
  140    DO 150 K = 1, L
            Z(I,K) = Z(I,K) / SCALE
            H = H + Z(I,K) * Z(I,K)
  150    CONTINUE
C
         F = Z(I,L)
         G = -SIGN(SQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         Z(I,L) = F - G
         F = 0.0
C
         DO 240 J = 1, L
            Z(J,I) = Z(I,J) / H
            G = 0.0
C     ********** FORM ELEMENT OF A*U **********
            DO 180 K = 1, J
  180       G = G + Z(J,K) * Z(I,K)
C
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
C
            DO 200 K = JP1, L
  200       G = G + Z(K,J) * Z(I,K)
C     ********** FORM ELEMENT OF P **********
  220       E(J) = G / H
            F = F + E(J) * Z(I,J)
  240    CONTINUE
C
         HH = F / (H + H)
C     ********** FORM REDUCED A **********
         DO 260 J = 1, L
            F = Z(I,J)
            G = E(J) - HH * F
            E(J) = G
C
            DO 260 K = 1, J
               Z(J,K) = Z(J,K) - F * E(K) - G * Z(I,K)
  260    CONTINUE
C
  290    D(I) = H
  300 CONTINUE
C
  320 D(1) = 0.0
      E(1) = 0.0
C     ********** ACCUMULATION OF TRANSFORMATION MATRICES **********
      DO 500 I = 1, N
         L = I - 1
         IF (D(I) .EQ. 0.0) GO TO 380
C
         DO 360 J = 1, L
            G = 0.0
C
            DO 340 K = 1, L
  340       G = G + Z(I,K) * Z(K,J)
C
            DO 360 K = 1, L
               Z(K,J) = Z(K,J) - G * Z(K,I)
  360    CONTINUE
C
  380    D(I) = Z(I,I)
         Z(I,I) = 1.0
         IF (L .LT. 1) GO TO 500
C
         DO 400 J = 1, L
            Z(I,J) = 0.0
            Z(J,I) = 0.0
  400    CONTINUE
C
  500 CONTINUE
C
      RETURN
C     ********** LAST CARD OF TRED2 **********
      END
C
      SUBROUTINE TRED3(N,NV,A,D,E,E2)
C
      INTEGER I,J,K,L,N,II,IZ,JK,NV
      REAL A(NV),D(N),E(N),E2(N)
      REAL F,G,H,HH,SCALE
C     REAL SQRT,ABS,SIGN
C     ********** FOR I=N STEP -1 UNTIL 1 DO -- **********
      DO  300 II = 1, N
         I = N + 1 - II
         L = I - 1
         IZ = (I * L) / 2
         H = 0.0
         SCALE = 0.0
         IF (L .LT. 1) GO TO 130
C     ********** SCALE ROW (ALGOL TOL THEN NOT NEEDED) **********
         DO 120 K = 1, L
            IZ = IZ + 1
            D(K) = A(IZ)
            SCALE = SCALE + ABS(D(K))
  120    CONTINUE
C
         IF (SCALE .NE. 0.0) GO TO 140
  130    E(I) = 0.0
         E2(I) = 0.0
         GO TO 290
C
  140    DO 150 K = 1, L
            D(K) = D(K) / SCALE
            H = H + D(K) * D(K)
  150    CONTINUE
C
         E2(I) = SCALE * SCALE * H
         F = D(L)
         G = -SIGN(SQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         D(L) = F - G
         A(IZ) = SCALE * D(L)
         IF (L .EQ. 1) GO TO 290
         F = 0.0
C
         DO 240 J = 1, L
            G = 0.0
            JK = (J * (J-1)) / 2
C     ********** FORM ELEMENT OF A*U **********
            DO 180 K = 1, L
               JK = JK + 1
               IF (K .GT. J) JK = JK + K - 2
               G = G + A(JK) * D(K)
  180       CONTINUE
C     ********** FORM ELEMENT OF P **********
            E(J) = G / H
            F = F + E(J) * D(J)
  240    CONTINUE
C
         HH = F / (H + H)
         JK = 0
C     ********** FORM REDUCED A **********
         DO 260 J = 1, L
            F = D(J)
            G = E(J) - HH * F
            E(J) = G
C
            DO 260 K = 1, J
               JK = JK + 1
               A(JK) = A(JK) - F * E(K) - G * D(K)
  260    CONTINUE
C
  290    D(I) = A(IZ+1)
         A(IZ+1) = SCALE * SQRT(H)
  300 CONTINUE
C
      RETURN
C     ********** LAST CARD OF TRED3 **********
      END

