! phonon diseprsion in GaN in general direcitons
!PROGRAM MAIN
SUBROUTINE GANPB(No,PSHR,FBD,FBD_HS,ISP,RI,QLs,CCL,CC1,DDR)
! ISP=0  (phonon dispersion),  1 (box fitting)
IMPLICIT NONE
INTEGER, PARAMETER :: Nq=97   ! # q points for phonon dispersion
INTEGER :: No,II,NN,PATH,PATH_I,PATH_L1,PATH_L2,PATH_L3,PATH_L4,QI_LOC,QI,QI_NUM,JN,JH,ISP,icoul
INTEGER :: iq,ion_ii,jj,K1,K2,I1,I2,J1,J2,J,IER,IS,KK_INX,MK_INX,GK_INX,AK_INX,HK_INX,iqs,IQ_SYM(7)
REAL*8 :: meV,aa_wz,cc_wz,u,eps,seps,VA,Z1,ZA(3),PR(3),RI(40,3),A_P(3),L_P(3),M_P(3),G_P(3),H_P(3),K_P(3)
REAL*8 :: QS_X,QS_Y,QS_Z,QE_X,QE_Y,QE_Z,Q_ST_X,Q_ST_Y,Q_ST_Z,QL_X(Nq),QL_Y(Nq),QL_Z(Nq),QL(3,Nq),QLs(3,No)
REAL*8 :: FREQ(12),FBD(12,Nq),FBD_HS(12,7),M(2),Qx,Qy,Qz,Q1,THZ,INVCM,OMEGA(12),PI,PSHR(37)
REAL*8 :: HR(12,12),HI(12,12),ZR(12,12),ZI(12,12),FV1(12),FV2(12),FM1(2,12)
COMPLEX*16 :: DR(12,12),DSR(4,4,3,3),DDR(4,3,3),DDD(4,3,3),CCQ(4,4,3,3),CC1(4,3,3),CCD(4,3,3),CCL(4,4,3,3,No),CI
!-------------------------------------------------------------------
!OPEN(UNIT=11,FILE='pshr.txt',STATUS='OLD')
!READ(11,*) (PSHR(II), II=1,38)
!WRITE(6,*) (PSHR(II), II=1,38)
!-------------------------------------------------------------------
PI = dacos(-1.d0) !3.1415926
CI=(0.d0,1.d0)

!aa_wz=3.189
!cc_wz=5.185
!u=0.3768
aa_wz=3.25  ! 3.189
cc_wz=5.207  ! 5.185
u=0.3817

eps=5.5
seps=sqrt(eps)
VA=sqrt(3.)*aa_wz**2*cc_wz/2.


Z1=1.1550   ! GaN
!Z1=0.0    ! ZnO
!Z1=PSHR(38)

!ZA(1)=Z1/seps/2.
!ZA(2)=Z1/seps/2.
!ZA(3)=Z1/seps/2.

ZA(1)=Z1
ZA(2)=Z1
ZA(3)=Z1
icoul=isp-1

! M(1)=69.723 !Ga
! M(2)=14.0067!N
M(1)=65.409 !Zn
M(2)=15.994 !O

!M(1)=1 !Ga
!M(2)=1 !N
!----- install DDR, CC1 to maintain translational inveriant ---------
Qx=0.
Qy=0.
Qz=0.
Q1=sqrt(Qx**2+Qy**2+Qz**2)
CALL GANSRD(aa_wz,cc_wz,QX,QY,QZ,RI,PSHR,DSR,DDR)
    JN=13
    JH=13
    do iqs=1,No
            Qx=QLs(1,iqs)
			Qy=QLs(2,iqs)
			Qz=QLs(3,iqs)
    if(isp==1)CCQ=CCL(:,:,:,:,iqs)      
	CALL GANSRD(aa_wz,cc_wz,QX,QY,QZ,RI,PSHR,DSR,DDD)
    if(isp==0)CALL GANCOL(JN,JH,QX,QY,QZ,aa_wz,cc_wz,RI,CCQ,CCD)	
!------- include self-interaciton----------
	DO ion_ii=1,4
		DO ii=1,3
			DO jj=1,3
				DSR(ion_ii,ion_ii,ii,jj)=DSR(ion_ii,ion_ii,ii,jj)+DDR(ion_ii,ii,jj)
				if(isp==0)CCQ(ion_ii,ion_ii,ii,jj)=CCQ(ion_ii,ion_ii,ii,jj)+CC1(ion_ii,ii,jj)
			ENDDO
		ENDDO
	ENDDO
	DO K1=1,4
        DO K2=1,4
            IF ((K1.EQ.1).OR.(K1.EQ.3)) THEN
                I1=1
            ELSE
                I1=2
            ENDIF
            IF ((K2.EQ.1).OR.(K2.EQ.3)) THEN
                I2=1
            ELSE
                I2=2
            ENDIF
            DO ii=1,3
                DO jj=1,3
                    J1=ii+(K1-1)*3
                    J2=jj+(K2-1)*3
                    DR(J1,J2)=DSR(K1,K2,ii,jj)+ZA(ii)*ZA(jj)*CCQ(K1,K2,ii,jj)
					!DR(J1,J2)=ZA(ii)*ZA(jj)*CCQ(K1,K2,ii,jj)
					!DR(J1,J2)=DSR(K1,K2,ii,jj)
                    DR(J1,J2)=DR(J1,J2)/sqrt(M(I1)*M(I2))
                ENDDO
            ENDDO
        ENDDO
    ENDDO
	!------------------ CH------------------------------------------
	DO I1=1,12
		DO I2=1,12
			J1=I1
			J2=I2
			HR(I1,I2) =  DR(J1,J2)
			HI(I1,I2) = -CI*DR(J1,J2)
		ENDDO
	ENDDO
    IS=0
    if(isp==0)write(6,*)' iq=',iqs
	CALL CH(12,12,HR,HI,OMEGA,IS,ZR,ZI,FV1,FV2,FM1,IER)
	IF(IER.NE.0)WRITE(6,*)'  EISPACK ERROR, IER=',IER 
	THZ=sqrt(14.4*9./(9.35*VA))*100./(2.*PI)
    INVCM=THZ*33.3556
    meV=THz/0.2418d0
    DO J=1,12
        IF (OMEGA(J).LT.0) THEN
            FREQ(J)=-SQRT(ABS(OMEGA(J)))*meV ! THz  ! INVCM
        ELSE
            FREQ(J)=SQRT(ABS(OMEGA(J)))*meV ! THz  ! INVCM
        ENDIF
        FBD(J,iqs)=FREQ(J)  !  for phonon dispersion
    ENDDO
	IF (ISP.EQ.1) THEN
		FBD_HS(:,iqs)=FBD(:,iqs)
	ENDIF
ENDDO ! iqs

7001 FORMAT(12F10.5)
7002 FORMAT(3F10.5)
     return
!END PROGRAM
END SUBROUTINE
