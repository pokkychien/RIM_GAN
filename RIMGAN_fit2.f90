! Fit phonon dispersion at special points
PROGRAM RIMGAN_fit
IMPLICIT NONE
INTEGER, PARAMETER :: M=37,N=37,K=M+1,Nq=97
INTEGER :: No,ii,QI_NUM,JN,JH,ISP,icoul,ion_ii
INTEGER :: Nk,i,ik,IGAMMA,IPRINT,ITMAX,ISEED,jj,kk,IT,IEV2,IEV1,FOP,iqs,IQ_SYM(7)
REAL*8 :: ALPHA,BETA,DELTA,aa_wz,cc_wz,dnom,Z1
REAL*8 :: X(K,M),R(K,M),F(K),VARY(M)
REAL*8 :: G(M),H(M),XC(M)
REAL*8 :: XIN(37),P1,P2,P3,aij(2,2),bij(2,2),vk(3),QL(3,Nq),Qx,Qy,Qz,Q1
REAL*8 :: PB_RIM(12,500),PB_HS(12,7),QLs(3,7),RI(40,3)
COMPLEX*16 :: DR(12,12),DSR(4,4,3,3),DDR(4,3,3),DDD(4,3,3)
COMPLEX*16 :: CI,CCQ(4,4,3,3),CC1(4,3,3),CCD(4,3,3),DD,CCL(4,4,3,3,7)
COMMON/PARK/ XIN,RI,QLs,CCL,DDR
!-------------------------------------------------------------------
!WRITE(6,*) 'FOP= (1: YES, 0: NO)'
!READ(5,*) FOP
!IF (FOP.EQ.1) THEN
    !OPEN(UNIT=11,file='input_GaN.txt',status='old')
!ELSEIF (FOP.EQ.0) THEN
   OPEN(UNIT=11,file='pshr.txt',status='old')
   OPEN(UNIT=12,file='RIM_wurzite.txt',status='old')
!ENDIF
do i=1,37
READ(11,*) XIN(i)
!XIN(i)=1.2/2.0*XIN(i)   ! rescale from GaN to ZnO
XIN(i)=XIN(i)   ! rescale from GaN to ZnO
enddo
READ(11,*)Z1
WRITE(6,*) 'XIN=' , XIN
    read(12,*)aa_wz    ! lattice constant (length of a2 for triangular lattice) in A
	read(12,*)aij(:,1)   ! a(i,j)  j=1,2 for a1 & a2 primitive vectors (in units of a)
    read(12,*)aij(:,2)
    read(12,*)bij(1,:)   ! b(i,j)  j=1,2 for b1 & b2 primitive vectors  (in units of 2*pi/a)
    read(12,*)bij(2,:)
aa_wz=3.25  ! 3.189
cc_wz=5.207  ! 5.185
!
    QL=0
    read(12,*)isp,Nk,dnom
    do ik=1,Nk
    read(12,*)vk
    QL(1:2,ik)=(vk(1)*bij(:,1)+vk(2)*bij(:,2))/dnom  ! in Cart. coord. 2*pi*(1/a,1/a,1/c)
    QL(3,ik)=vk(3)/dnom 
    enddo
NO=Nq
if(isp==1)No=7
   call get_RI(aa_wz,cc_wz,RI)
   if(isp==0)call Set_Qpnt(ISP,No,aa_wz,cc_wz,RI,QL)
   if(isp==1)call Set_Qpnt(ISP,No,aa_wz,cc_wz,RI,QLs)
!----- install DDR, CC1 to maintain translational inveriant ---------
Qx=0.
Qy=0.
Qz=0.
Q1=sqrt(Qx**2+Qy**2+Qz**2)
CALL GANSRD(aa_wz,cc_wz,QX,QY,QZ,RI,Xin,DSR,DDR)
    JN=13
    JH=13
CALL GANCOL(JN,JH,QX,QY,QZ,aa_wz,cc_wz,RI,CCQ,CC1)
if(isp==1)then
do iqs=1,7
    Qx=QLs(1,iqs)
    Qy=QLs(2,iqs)
    Qz=QLs(3,iqs)
    CALL GANCOL(JN,JH,QX,QY,QZ,aa_wz,cc_wz,RI,CCQ,CCD)	
!------- include self-interaciton----------
	DO ion_ii=1,4
		DO ii=1,3
			DO jj=1,3
				CCQ(ion_ii,ion_ii,ii,jj)=CCQ(ion_ii,ion_ii,ii,jj)+CC1(ion_ii,ii,jj)
			ENDDO
		ENDDO
    ENDDO
    CCL(:,:,:,:,iqs)=CCQ
enddo  ! iqs
!    QLs=QL(:,1:No)
ALPHA=1.3
BETA=0.0001
IGAMMA=5
DELTA=0.00001
! ------IPRINT : 0 (NO PRINT), 1(PRINT EVERY 50 ITERATION), 2(PRINT EVERY ITER.)
IPRINT=2
ITMAX=1200

! GENERATE RANDOM NUMBER ARRAY
ISEED=16009
DO ii=1,N
    DO jj=2,K
        R(jj,ii)=RAN(ISEED)
    END DO
END DO
! INITIAL GUESS OF THE INPUT PARAMETERS (RIM PARAMETER)
DO II=1,M
    X(1,II)=XIN(II)
END DO
!----------
if(isp==0)CALL GANPB(No,XIN,PB_RIM,PB_HS,ISP,RI,QL,CCL,CC1,DDR)
if(isp==1)CALL GANPB(No,XIN,PB_RIM,PB_HS,ISP,RI,QLs,CCL,CC1,DDR)
!write(6,*) 'PB_HS', PB_HS(:,3)
DO ii=1,NO
WRITE (21,7001)float(ii-1), (PB_RIM(jj,ii),jj=1,12) !  f (in cm^-1)
ENDDO
! DO ii=1,4
! WRITE (22,7001) (PB_HS(jj,ii),jj=1,12) !  f (in cm^-1)
! ENDDO

!STOP
!--------------------------
if(isp==1)CALL CONSX(N,M,K,ITMAX,ALPHA,BETA,IGAMMA,DELTA,X,R,F,IT,IEV2,6,G,H,XC,IPRINT)
end if
if(isp==0)XC=XIN
  if(isp==0)CALL GANPB(No,XC,PB_RIM,PB_HS,ISP,RI,QL,CCL,CC1,DDR)
  if(isp==1)CALL GANPB(No,XC,PB_RIM,PB_HS,ISP,RI,QLs,CCL,CC1,DDR)

 DO ii=1,Nq
 WRITE (21,7001) float(ii-1),(PB_RIM(jj,ii),jj=1,12) !  f (in cm^-1)
 ENDDO
if(isp==1)then
WRITE(6,*) 'IT=',IT-1
WRITE(6,*) F
DO ii=1,7
WRITE (22,7001) (PB_HS(jj,ii),jj=1,12) !  f (in cm^-1)
ENDDO
WRITE(30,7002) XC
endif

CLOSE(UNIT=11)
CLOSE(UNIT=12)

7001 FORMAT(13F10.5)
7002 FORMAT(F12.6)
    END PROGRAM
!*************************************************************
    subroutine Set_Qpnt(iSP,No,aa_wz,cc_wz,RI,QLs)
IMPLICIT NONE
INTEGER, PARAMETER :: Nq=97   ! # q points for phonon dispersion
INTEGER :: ISP,No,ii,iq,iqs,NN,PATH,PATH_I,PATH_L1,PATH_L2,PATH_L3,PATH_L4,QI_LOC,QI,QI_NUM,JN,JH,icoul
INTEGER :: KK_INX,MK_INX,GK_INX,AK_INX,HK_INX,IQ_SYM(No)
real*8 aa_wz,cc_wz,RI(40,3),QLs(3,No),Qx,Qy,Qz,Q1,A_P(3),L_P(3),M_P(3),G_P(3),H_P(3),K_P(3)
REAL*8 :: QS_X,QS_Y,QS_Z,QE_X,QE_Y,QE_Z,Q_ST_X,Q_ST_Y,Q_ST_Z,QL_X(Nq),QL_Y(Nq),QL_Z(Nq),QL(3,97)
A_P(1)=0.
A_P(2)=0.
A_P(3)=1./2.

L_P(1)=1./2.
L_P(2)=1./(2.*sqrt(3.))
L_P(3)=1./2.

M_P(1)=1./2.
M_P(2)=1./(2.*sqrt(3.))
M_P(3)=0.

G_P(1)=0.
G_P(2)=0.
G_P(3)=0.

H_P(1)=2./3.
H_P(2)=0.
H_P(3)=1./2.

K_P(1)=2./3.
K_P(2)=0.
K_P(3)=0.
NN=20   ! 100

PATH_L1=NN
PATH_L2=NINT(NN/2.)
PATH_L3=NINT(NN*sqrt(3.)/2.)
PATH_L4=NINT(NN*3./4.*aa_wz/cc_wz)
QI_NUM=PATH_L1*2+PATH_L2*2+PATH_L3*2+PATH_L4-6
QI_LOC=1
DO PATH_I=1,7    !  k pnts for bands
	SELECT CASE (PATH_I)
	CASE (1)
		QS_X=G_P(1)
		QE_X=K_P(1)
		QS_Y=G_P(2)
		QE_Y=K_P(2)
		QS_Z=G_P(3)
		QE_Z=K_P(3)
		PATH=PATH_L1
	CASE (2)
		QS_X=K_P(1)
		QE_X=M_P(1)
		QS_Y=K_P(2)
		QE_Y=M_P(2)
		QS_Z=K_P(3)
		QE_Z=M_P(3)
		PATH=PATH_L2
	CASE (3)
		QS_X=M_P(1)
		QE_X=G_P(1)
		QS_Y=M_P(2)
		QE_Y=G_P(2)
		QS_Z=M_P(3)
		QE_Z=G_P(3)
		PATH=PATH_L3
	CASE (4)
		QS_X=G_P(1)
		QE_X=A_P(1)
		QS_Y=G_P(2)
		QE_Y=A_P(2)
		QS_Z=G_P(3)
		QE_Z=A_P(3)
		PATH=PATH_L4
	CASE (5)
		QS_X=A_P(1)
		QE_X=H_P(1)
		QS_Y=A_P(2)
		QE_Y=H_P(2)
		QS_Z=A_P(3)
		QE_Z=H_P(3)
		PATH=PATH_L1
	CASE (6)
		QS_X=H_P(1)
		QE_X=L_P(1)
		QS_Y=H_P(2)
		QE_Y=L_P(2)
		QS_Z=H_P(3)
		QE_Z=L_P(3)
		PATH=PATH_L2
	CASE (7)
		QS_X=L_P(1)
		QE_X=A_P(1)
		QS_Y=L_P(2)
		QE_Y=A_P(2)
		QS_Z=L_P(3)
		QE_Z=A_P(3)
		PATH=PATH_L3
	END SELECT
	Q_ST_X=(QE_X-QS_X)/(PATH-1)
	Q_ST_Y=(QE_Y-QS_Y)/(PATH-1)
	Q_ST_Z=(QE_Z-QS_Z)/(PATH-1)
	IF (PATH_I.NE.7) PATH=PATH-1
	DO QI=1,PATH
	    QL_X(QI_LOC)=QS_X+(QI-1)*Q_ST_X
		QL_Y(QI_LOC)=QS_Y+(QI-1)*Q_ST_Y
		QL_Z(QI_LOC)=QS_Z+(QI-1)*Q_ST_Z
		QI_LOC=QI_LOC+1
	ENDDO
ENDDO
!----- SPECIAL POINT K-INDEX---------
KK_INX=PATH_L1
MK_INX=PATH_L1+PATH_L2-1
GK_INX=PATH_L1+PATH_L2+PATH_L3-2
AK_INX=PATH_L1+PATH_L2+PATH_L3+PATH_L4-3
HK_INX=2*PATH_L1+PATH_L2+PATH_L3+PATH_L4-4

IQ_SYM(1)=KK_INX
IQ_SYM(2)=MK_INX
IQ_SYM(3)=GK_INX-1  ! for discontinuity along Gamma to M
IQ_SYM(4)=GK_INX
IQ_SYM(5)=AK_INX
IQ_SYM(6)=HK_INX
IQ_SYM(7)=KK_INX/2
!WRITE(6,*) QL_X(1:QI_NUM)

QL(1,:)=QL_X
QL(2,:)=QL_Y
QL(3,:)=QL_Z

!------------------------------------  
JN=13
JH=13
!----- instal DDR, CC1 to maintain translational inveriant ---------
!--------installation completed-------------------------------------
QI_NUM=No     ! for phonon bands
!DO iq=1,QI_NUM
!DO iq=235,235
DO iqs=1,QI_NUM
	IF (iSP.EQ.1) THEN
	iq=IQ_SYM(iqs)
	ELSE
	iq=iqs
	ENDIF
	IF (iSP.EQ.1) THEN
		IF (iqs.EQ.5) THEN
			QLs(1,iqs)=1./2.*1E-4
			QLs(2,iqs)=1./(2.*sqrt(3.))*1E-4
			QLs(3,iqs)=0.
		ELSE
			QLs(1,iqs)=QL(1,iq)
			QLs(2,iqs)=QL(2,iq)
			QLs(3,iqs)=QL(3,iq)
        ENDIF
            Qx=QLs(1,iqs)
			Qy=QLs(2,iqs)
			Qz=QLs(3,iqs)
	ELSE
     QLs(:,iq)=QL(:,iq)
    ENDIF
END DO
     END
!*************************************************************
subroutine get_RI(aa_wz,cc_wz,RI)
IMPLICIT NONE
INTEGER :: ii
real*8 u,aa_wz,cc_wz,PR(3),RI(40,3)
u=0.3768
DO ii=1,40
    SELECT CASE (ii)
	CASE (1)
	    PR(1)=0.
		PR(2)=0.
		PR(3)=0.
	CASE (2)
		PR(1)=0.
		PR(2)=0.
		PR(3)=cc_wz*u
	CASE (3)
		PR(1)=0.
		PR(2)=sqrt(3.)/3.*aa_wz
		PR(3)=cc_wz*1./2.
	CASE (4)
		PR(1)=0.
		PR(2)=sqrt(3.)/3.*aa_wz
		PR(3)=cc_wz*(1./2.+u)
	CASE (10)
		PR(1)=0.
		PR(2)=0.
		PR(3)=cc_wz*u-cc_wz
	CASE (11)
		PR(1)=0.
		PR(2)=sqrt(3.)/3.*aa_wz
		PR(3)=cc_wz*1./2.-cc_wz
	CASE (12)
		PR(1)=0.
		PR(2)=sqrt(3.)/3.*aa_wz
		PR(3)=cc_wz*(1./2.+u)-cc_wz
	CASE (13)
		PR(1)=aa_wz
		PR(2)=0.
		PR(3)=0.
	CASE (16)
		PR(1)=aa_wz
		PR(2)=0.
		PR(3)=cc_wz*u
	CASE (19)
		PR(1)=aa_wz
		PR(2)=sqrt(3.)/3.*aa_wz
		PR(3)=cc_wz*1./2.
	CASE (22)
		PR(1)=aa_wz
		PR(2)=sqrt(3.)/3*aa_wz
		PR(3)=cc_wz*(1./2.+u)
	CASE (37)
	    PR(1)=3/2*aa_wz
		PR(2)=-sqrt(3.)/6.*aa_wz
		PR(3)=-cc_wz*(1./2.-u)
	CASE (38)
		PR(1)=3/2*aa_wz
		PR(2)=-sqrt(3.)/6.*aa_wz
		PR(3)=1./2.*cc_wz
	CASE  DEFAULT
		PR(1:3)=0.
	END SELECT
	RI(ii,:)=PR    !  Lattice vectors for CC
    ENDDO
    end
! =================================================================================================
    SUBROUTINE FUNC(N,M,K,X,F,I) ! for ZnO in units of meV
IMPLICIT NONE
INTEGER :: Nq,N,M,K,I,J,IB,IK,ISK,ii,jj
REAL*8 :: X(K,M),F(K),PV(N),XIN(37)
REAL*8 :: D1,D2,D3,D4,D5,D6,DD,QLs(3,7),PB_HS(12,7),PB_RIM(12,500)  ! ,QLs(3,6)
REAL*8 :: D7,D8,D9,D10,D11,D12,D13,D14,D15,D16,D17,D18,D19,D20,D21,D22,D23,D24,D25,RI(40,3)
complex*16:: CCL(4,4,3,3,7),CC1(4,3,3),DDR(4,3,3)
COMMON/PARK/ XIN,RI,QLs,CCL,DDR

DO J=1,N
    PV(J)=X(I,J)
END DO
CC1=0
CALL GANPB(7,PV,PB_RIM,PB_HS,1,RI,QLs,CCL,CC1,DDR)
DO ii=1,7
WRITE (40,7001) (PB_HS(jj,ii),jj=1,12) !  f (in cm^-1)
ENDDO
7001 FORMAT(13F10.5)
DD=0.
!DO IK=3,3
    !------ Gamma point + ---------
    D1=(PB_HS(4,4)/2+PB_HS(5,4)/2-11.22)**2
    D2=(PB_HS(6,4)-32.21)**2
    D3=(PB_HS(7,4)/2.+PB_HS(8,4)/2. -50.71)**2
    D4=(PB_HS(9,4)/2.+PB_HS(10,4)/2.-54.38)**2
    D5=(PB_HS(11,4)-68.32)**2
    D6=(PB_HS(12,4)-68.91)**2
	!------ Gamma point - ---------
	D7=(PB_HS(7,3)-48.36)**2
	D8=(PB_HS(12,3)-70.08)**2
	IF (PB_HS(1,3).LE.0.) THEN
        D9=10000
    ELSE
        D9=0.
    ENDIF
	IF (PB_HS(1,7).LE.0.) THEN
        D21=10000
    ELSE
        D21=0.
    ENDIF
    !------ K point---------
    D10=(PB_HS(1,1)-18.06)**2
    D11=(PB_HS(7,1)/2.+PB_HS(8,1)/2.-59.01)**2
    !------ M point---------
    D12=(PB_HS(1,2)-10.45)**2
    D13=(PB_HS(2,2)-15.45)**2
    D14=(PB_HS(3,2)-16.33)**2
    D15=(PB_HS(4,2)-20.73)**2
    !D16=(PB_HS(7,2)/2.+PB_HS(8,2)/2.-576.)**2
    !------ A point---------
    D17=(PB_HS(5,5)/2.+PB_HS(6,5)/2.-23.29)**2
    D18=(PB_HS(11,5)/2.+PB_HS(12,5)/2.-68.79)**2
	!---------------
!ENDDO
!write(6,*) 'DD=',DD
DD=D1+D2+D3+D4+D5+D6+D7+D8+D10+D11+D12+D13+D14+D15+D17+D18+D9 
F(I)=-DD
if(I==1)write(6,*)PB_HS(1:3,4)
return
END SUBROUTINE
! =================================================================================================
SUBROUTINE CONST(N,M,K,X,G,H,L)
INTEGER :: ii
REAL*8 :: X(K,M),G(M),H(M)
COMMON/PARK/ PAR(37)

DO ii=1,N
    G(ii)=PAR(ii)*0.9-0.5
    H(ii)=PAR(ii)*1.1+0.5
ENDDO
    END