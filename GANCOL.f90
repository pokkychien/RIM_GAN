SUBROUTINE GANCOL(JN,JH,QX,QY,QZ,aa_wz,cc_wz,RI,CC,CC1)
IMPLICIT NONE
INTEGER :: K1,K2,KK2,ii,jj,JN,JH,N1,N2,N3,M1,M2,M3
REAL*8 :: PI, a1(3),a2(3),a3(3),aa_wz,cc_wz,PL(3,3),b1(3),b2(3),b3(3),PRL(3,3),VA,ETA,Q(3),QQ,VQ(3),Q1,Q2,RR(3),RI(40,3)
REAL*8 :: ECUT,RCUT,G(3),GQ(3),G2,GQ2,SS,LAT(3),RL(3),RL2,X(3),X2,X1,H,H1,H2,VVQ,SS1
REAL*8 :: QX,QY,QZ
COMPLEX*16 :: CI,CC(4,4,3,3),CC1(4,3,3),DD,PHI1,PHI2,HAB
!-------------------------------------------------------------------
PI = dacos(-1.d0) !3.1415926
CI=(0.d0,1.d0)
!---- primitive lattice vectors
a1(1)=aa_wz
a1(2)=0.
a1(3)=0.

a2(1)=aa_wz*1./2.
a2(2)=aa_wz*sqrt(3.)/2.
a2(3)=0.

a3(1)=0.
a3(2)=0.
a3(3)=cc_wz

PL(1:3,1)=a1
PL(1:3,2)=a2
PL(1:3,3)=a3
!---- repiprical lattice vectors
b1(1)=2*pi/aa_wz
b1(2)=-2*pi/aa_wz*1./sqrt(3.)
b1(3)=0.

b2(1)=0.
b2(2)=2*pi/aa_wz*2./sqrt(3.)
b2(3)=0.

b3(1)=0.
b3(2)=0.
b3(3)=2*pi/cc_wz

PRL(1:3,1)=b1
PRL(1:3,2)=b2
PRL(1:3,3)=b3

VA=sqrt(3.)*aa_wz**2*cc_wz/2.! volume of unit cell
ETA=sqrt(2.*pi/aa_wz*2.*PI/cc_wz) ! ETA=ALPHA^2. ALPHA is the adjustable parameters in Ewald's method

Q(1)=2*PI/aa_wz*Qx
Q(2)=2*PI/aa_wz*Qy
Q(3)=2*PI/cc_wz*Qz

QQ=sqrt(Qx**2+Qy**2+Qz**2)
!-----------------dont understand------------------------------------
IF (QQ.NE.0) THEN
    VQ(1)=Qx/QQ
    VQ(2)=Qy/QQ
    VQ(3)=Qz/QQ
ELSE
    VQ(1)=0
    VQ(2)=0
    VQ(3)=1
ENDIF
!-----------------dont understand------------------------------------
!--------------calcualte coulumb matrix 
CC(:,:,:,:)=0.
DO K1=1,4
!DO K1=1,1
    IF ((K1.EQ.1).OR.( K1.EQ.3)) THEN
        Q1=1
    ELSE
        Q1=-1
    ENDIF
    DO K2=K1,4
	!DO K2=4,4
        IF ((K2.EQ.1).OR.( K2.EQ.3)) THEN
            Q2=1
        ELSE
            Q2=-1
        ENDIF
    RR=RI(K2,:)-RI(K1,:)
    DO ii=1,3
        DO jj=1,ii
              !--CALCULATE THE SUM IN RECIPROCAL SPACE
              ECUT=4.*(JN*2.*PI/aa_wz)**2/3.
              PHI1=0.
              DO N1= -JN,JN
                  DO N2= -JN,JN
                      DO N3= -2*JN,2*JN
						  G(1)=N1*PRL(1,1)+N2*PRL(1,2)+N3*PRL(1,3)
						  G(2)=N1*PRL(2,1)+N2*PRL(2,2)+N3*PRL(2,3)
						  G(3)=N1*PRL(3,1)+N2*PRL(3,2)+N3*PRL(3,3)
                          G2=G(1)**2+G(2)**2
                          IF (G2.LT.ECUT) THEN
                              GQ=G+Q
                              GQ2=GQ(1)**2+GQ(2)**2+GQ(3)**2
                              IF (GQ2.GT.1E-10) THEN
                                  SS=G(1)*RR(1)+G(2)*RR(2)+G(3)*RR(3)
                                  PHI1=PHI1+GQ(ii)*GQ(jj)*4*PI/GQ2*EXP(-GQ2/4./ETA)*EXP(-CI*SS);
                              ELSE
                                  VVQ = VQ(1)**2 +VQ(2)**2 + VQ(3)**2
                                  PHI1=PHI1+VQ(ii)*VQ(jj)*4*PI/VVQ
                              ENDIF
                          ENDIF ! ECUT
                      ENDDO ! N3
                  ENDDO ! N2
              ENDDO ! N1
              ! EVALUATE CORRECTION SUM IN REAL SPACE
              RCUT=(JH*aa_wz)**2
              PHI2=0
              DO M1= -JH,JH
                  DO M2= -JH,JH
                      DO M3= -2*JH,2*JH
						  ! LAT=M1*a1+M2*a2+M3*a3;
						  LAT(1)=M1*PL(1,1)+M2*PL(1,2)+M3*PL(1,3)
						  LAT(2)=M1*PL(2,1)+M2*PL(2,2)+M3*PL(2,3)
						  LAT(3)=M1*PL(3,1)+M2*PL(3,2)+M3*PL(3,3)
                          RL2=LAT(1)**2+LAT(2)**2
                          IF (RL2.LT.RCUT) THEN
                            RL=-LAT+RR
                            X=sqrt(ETA)*RL
                            X2=X(1)**2+X(2)**2+X(3)**2
                            X1=sqrt(X2)
                            IF (X1.GT.sqrt(ETA)*aa_wz/100.) THEN
                                H1=3/(X1**3)*(1-erf(X1))+2./sqrt(PI)*(3./X2+2)*EXP(-X2)
                                H1=X(ii)*X(jj)*H1/X2
                                IF (ii.EQ.jj) THEN
                                    H2=1./(X1**3)*(1-erf(X1))+2./sqrt(PI)/X2*EXP(-X2)
                                ELSE
                                    H2=0.
                                ENDIF
                                H=H1-H2
                            ELSEIF (ii.EQ.jj) THEN
                                H=4./3./sqrt(PI)
                            ELSE
                                H=0.
                            ENDIF
                            SS1=Q(1)*RL(1)+Q(2)*RL(2)+Q(3)*RL(3)
                            HAB=-ETA**(1.5)*H*CDEXP(CI*SS1)
                            PHI2=PHI2+HAB
                          ENDIF ! RCUT
                      ENDDO ! M3
                  ENDDO ! M2
              ENDDO ! M1 
              CC(K1,K2,ii,jj) =(PHI1+PHI2*VA)*Q1*Q2
        ENDDO ! jj
    ENDDO ! ii
    CC(K1,K2,1,2)=CC(K1,K2,2,1)
	CC(K1,K2,1,3)=CC(K1,K2,3,1)
	CC(K1,K2,2,3)=CC(K1,K2,3,2)
    ENDDO ! K2
ENDDO ! K1
DO K2=1,4
    KK2=K2+1
    DO K1=KK2,4
        DO ii=1,3
            DO jj=1,3
                CC(K1,K2,ii,jj) = CONJG(CC(K2,K1,jj,ii))
            ENDDO
        ENDDO
    ENDDO
ENDDO
!------ CONSTRCUT CC1 AT ZONE-CENTER -----
CC1(:,:,:)=0.
IF (QQ.LT.0.00001) THEN
    DO K1=1,4
        DO ii=1,3
            DO jj=1,3
                DD=0.
                DO K2=1,4
                    DD=DD-CC(K1,K2,ii,jj)
                ENDDO
                CC1(K1,ii,jj)=DD
            ENDDO
        ENDDO
    ENDDO
ENDIF


END SUBROUTINE