! Fit phonon dispersion at special points
PROGRAM RIMGAN_band
IMPLICIT NONE
INTEGER, PARAMETER :: M=37,N=37,K=M+1,Nq=51
INTEGER :: IGAMMA,IPRINT,ITMAX,ISEED,ii,jj,kk,IT,IEV2,IEV1,NO,FOP
REAL*8 :: ALPHA,BETA,DELTA
REAL*8 :: X(K,M),R(K,M),F(K),VARY(M)
REAL*8 :: G(M),H(M),XC(M)
REAL*8 :: XIN(37),P1,P2,P3
REAL*8 :: PB_RIM(12,500),PB_HS(12,Nq)
COMMON/PARK/ XIN,QLs(3,7)
!-------------------------------------------------------------------
!WRITE(6,*) 'FOP= (1: YES, 0: NO)'
!READ(5,*) FOP
!IF (FOP.EQ.1) THEN
    !OPEN(UNIT=11,file='input_GaN.txt',status='old')
!ELSEIF (FOP.EQ.0) THEN
   OPEN(UNIT=11,file='pshr.txt',status='old')
!ENDIF
do i=1,M
READ(11,*) XIN(i)
enddo
WRITE(6,*) 'XIN=' , XIN
    read(11,*)aL    ! lattice constant (length of a2 for triangular lattice) in A
	read(11,*)aij(:,1)   ! a(i,j)  j=1,2 for a1 & a2 primitive vectors (in units of a)
    read(11,*)aij(:,2)
    read(11,*)bij(1,:)   ! b(i,j)  j=1,2 for b1 & b2 primitive vectors  (in units of 2*pi/a)
    read(11,*)bij(2,:)
!--------------------------
XC=XIN
isp=0
if(isp==1)CALL CONSX(N,M,K,ITMAX,ALPHA,BETA,IGAMMA,DELTA,X,R,F,IT,IEV2,NO,G,H,XC,IPRINT)
CALL GANPB(Nq,XC,PB_RIM,PB_HS,0)

DO ii=1,Nq
WRITE (22,7001) float(Nq-1) (PB_HS(jj,ii),jj=1,12) !  f (in cm^-1)
ENDDO

CLOSE(UNIT=11)


7001 FORMAT(13F10.5)
7002 FORMAT(1F10.5)

END PROGRAM
! =================================================================================================
SUBROUTINE FUNC(N,M,K,X,F,I)
IMPLICIT NONE
INTEGER :: N,M,K,I,J,IB,IK,ISK
REAL*8 :: X(K,M),F(K),PV(N),XIN(37)
REAL*8 :: D1,D2,D3,D4,D5,D6,DD,PB_HS(12,7),PB_RIM(12,5000)
REAL*8 :: D7,D8,D9,D10,D11,D12,D13,D14,D15,D16,D17,D18,D19,D20
COMMON/PARK/ XIN,QLs(3,7)

DO J=1,N
    PV(J)=X(I,J)
END DO

CALL GANPB(PV,PB_RIM,PB_HS,1)
DD=0.
!DO IK=3,3
    !------ Gamma point + ---------
    D1=(PB_HS(4,4)-144.)**2
    D2=(PB_HS(6,4)-329.)**2
    D3=(PB_HS(7,4)/2.+PB_HS(8,4)/2. -558.8)**2
    D4=(PB_HS(9,4)/2.+PB_HS(10,4)/2.-567.6)**2
    D5=(PB_HS(11,4)-692.)**2
    D6=(PB_HS(12,4)-734.)**2
	!------ Gamma point - ---------
	D7=(PB_HS(7,3)-531.8)**2
	D8=(PB_HS(12,3)-741)**2
	IF (PB_HS(1,3).LE.0.) THEN
        D9=1000
    ELSE
        D9=0.
    ENDIF
    !------ K point---------
    D10=(PB_HS(1,1)-215.)**2
    D11=(PB_HS(7,1)/2.+PB_HS(8,1)/2.-614.)**2
	D20=(PB_HS(7,1)-PB_HS(8,1))**2
	D19=(PB_HS(9,1)-642.)**2
    !------ M point---------
    D12=(PB_HS(1,2)-137.)**2
    D13=(PB_HS(2,2)-184.)**2
    D14=(PB_HS(3,2)-193.)**2
    D15=(PB_HS(4,2)-238.)**2
    !D16=(PB_HS(7,2)/2.+PB_HS(8,2)/2.-576.)**2
    !------ A point---------
    D17=(PB_HS(5,5)/2.+PB_HS(6,5)/2.-231.)**2
    D18=(PB_HS(11,5)/2.+PB_HS(12,5)/2.-711.)**2
	!---------------
!ENDDO
!write(6,*) 'DD=',DD
!DD=(D1*2+D2+D3+D4+D5+D6)*2+D7+D8+D9+D10+D11+D12+D13+(D14+D15)
DD=(D1*5+D2+D3+D4+D5+D6+D7+D8+D9)*2+(D10+D11+D19*0.75+D20*2)+(D12*0.75+D13+D14+D15)+(D17+D18)
F(I)=-DD

END SUBROUTINE

! =================================================================================================
SUBROUTINE CONST(N,M,K,X,G,H,L)
INTEGER :: ii
REAL*8 :: X(K,M),G(M),H(M)
COMMON/PARK/ PAR(37),QLs(3,7)

DO ii=1,N
    G(ii)=PAR(ii)-1.
    H(ii)=PAR(ii)+1.
ENDDO
END
