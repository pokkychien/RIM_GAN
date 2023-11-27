SUBROUTINE GANSRD(AA_WZ,CC_WZ,Qx,Qy,Qz,RI,PSHR,DSR,DDR)
!   RI= k-vecs sampling 
IMPLICIT NONE
INTEGER :: NEAR_II,II,JJ,ION_II,ION_JJ,K1,K2
REAL*8 :: AA_WZ,CC_WZ,RI(40,3),RJ(3),RL(3),XI(3),XL(3),PI
REAL*8 :: A1,B1,A2,B2,C2,D2,A3,B3,C3,A4,B4,C4,D4,E4,A5,B5,C5,D5
REAL*8 :: A6,A7,B7,C7,D7,E7,F7,A8,B8,C8,D8,E8,F8,A9,B9,C9,D9,E9,F9
REAL*8 :: PSHR(37),Q2,Q1,Qx,Qy,Qz,DOT,DOT_RL
REAL*8 :: II_OP(3,3),C3_OP(3,3),SIGMA_V1(3,3),SIGMA_V2(3,3),SIGMA_V3(3,3),SYM_OP(3,3)
REAL*8 :: DD_TEMP(3,3),D_12(3,3),D_34(3,3),D_13(3,3),D_14(3,3),D_23b(3,3),D_23(3,3),D_14b(3,3),D_24(3,3),D_12b(3,3)
REAL*8 :: D_11p(3,3),D_22p(3,3),D_33p(3,3),D_44p(3,3),D_14_COP(3,3),D_23_COP(3,3)
COMPLEX*16 :: CPH,CPH_RL,Dsr_12(3,3),Dsr_34(3,3),Dsr_13(3,3),Dsr_13b(3,3),Dsr_23b(3,3),Dsr_23(3,3),Dsr_14b(3,3),Dsr_14(3,3)
COMPLEX*16 :: Dsr_24(3,3),Dsr_12b(3,3),Dsr_34b(3,3),Dsr_11p(3,3),Dsr_22p(3,3),Dsr_33p(3,3),Dsr_44p(3,3),Dsr_24b(3,3)
COMPLEX*16 :: Dsr_14_COP(3,3),Dsr_23_COP(3,3),DSR(4,4,3,3),DDR(4,3,3),DD
COMPLEX*16 :: CI
!-------------------------------------------------------------------
PI = dacos(-1.d0) !3.1415926
CI=(0.d0,1.d0)
!-------------------------------------------------------------------
A1=PSHR(1)
B1=PSHR(2)

A2=PSHR(3)
B2=PSHR(4)
C2=PSHR(5)
D2=PSHR(6)


A3=PSHR(7)
B3=PSHR(8)
C3=PSHR(9)

A4=PSHR(10)
B4=PSHR(11)
C4=PSHR(12)
D4=PSHR(13)
E4=PSHR(14)

A5=PSHR(15)
B5=PSHR(16)
C5=PSHR(17)
D5=PSHR(18)


A6=PSHR(19)

A7=PSHR(20)
B7=PSHR(21)
C7=PSHR(22)
D7=PSHR(23)
E7=PSHR(24)
F7=PSHR(25)

A8=PSHR(26)
B8=PSHR(27)
C8=PSHR(28)
D8=PSHR(29)
E8=PSHR(30)
F8=PSHR(31)

A9=PSHR(32)
B9=PSHR(33)
C9=PSHR(34)
D9=PSHR(35)
E9=PSHR(36)
F9=PSHR(37)
!-------------------------------------------------------------------
Q2=Qx**2+Qy**2+Qz**2
Q1=sqrt(Q2)
DSR(:,:,:,:)=0.
II_OP(:,:)=0.
C3_OP(:,:)=0.
sigma_v1(:,:)=0.
sigma_v2(:,:)=0.
sigma_v3(:,:)=0.

II_OP(1,1)= 1
II_OP(2,2)= 1
II_OP(3,3)= 1

C3_OP(1,1)= -1./2.
C3_OP(2,2)= -1./2.
C3_OP(3,3)= 1.
C3_OP(1,2)= -sqrt(3.)/2
C3_OP(2,1)=  sqrt(3.)/2
!-------C3
!      -1/2,   -sqrt(3)/2,     0
! sqrt(3)/2,         -1/2,     0
!         0,            0,     1

!-------SGV1
!sigma_v1=[-1,0,0;0,1,0;0,0,1]
sigma_v1(1,1)= -1
sigma_v1(2,2)=  1
sigma_v1(3,3)=  1
!-------SGV2
!sigma_v1=[1/2,sqrt(3)/2,0;sqrt(3)/2,-1/2,0;0,0,1]
sigma_v2(1,1)=1./2.
sigma_v2(1,2)=sqrt(3.)/2.
sigma_v2(2,1)=sqrt(3.)/2.
sigma_v2(2,2)=-1./2.
sigma_v2(3,3)= 1.
!-------SGV3
!sigma_v1=[1/2,-sqrt(3)/2,0;-sqrt(3)/2,-1/2,0;0,0,1]
sigma_v3(1,1)=1./2.
sigma_v3(1,2)=-sqrt(3.)/2.
sigma_v3(2,1)=-sqrt(3.)/2.
sigma_v3(2,2)=-1./2.
sigma_v3(3,3)= 1.
!===================================================================
Dsr_12(:,:)=0.
D_12(1,1)=A1
D_12(2,2)=A1
D_12(3,3)=B1
!D_12= A1,  0,  0,
!      0,  A1,  0,
!      0,   0,  B1,
RJ=RI(2,:)-RI(1,:)
DO near_ii=1,1
    SELECT CASE (near_ii)
        CASE (1) 
            SYM_OP=II_op
    END SELECT
    XI(1)=SYM_OP(1,1)*RJ(1)+SYM_OP(1,2)*RJ(2)+SYM_OP(1,3)*RJ(3)
	XI(2)=SYM_OP(2,1)*RJ(1)+SYM_OP(2,2)*RJ(2)+SYM_OP(2,3)*RJ(3)
	XI(3)=SYM_OP(3,1)*RJ(1)+SYM_OP(3,2)*RJ(2)+SYM_OP(3,3)*RJ(3)
    DOT=Qx*XI(1)/aa_wz+Qy*XI(2)/aa_wz+Qz*XI(3)/cc_wz
    CPH=EXP(CI*DOT*2*PI)
    
    DD_temp=MATMUL(MATMUL(SYM_OP,D_12),TRANSPOSE(SYM_OP))
    Dsr_12=Dsr_12+DD_temp*CPH
ENDDO
!-------------------------------------------------------------------
Dsr_34(:,:)=0.
D_34(1,1)=A1
D_34(2,2)=A1
D_34(3,3)=B1
!D_34= A1,  0,  0,
!      0,  A1,  0,
!      0,   0,  B1,
RJ=RI(4,:)-RI(3,:);
DO near_ii=1,1
    SELECT CASE (near_ii)
        CASE (1) 
            SYM_OP=II_op;
    END SELECT
    XI(1)=SYM_OP(1,1)*RJ(1)+SYM_OP(1,2)*RJ(2)+SYM_OP(1,3)*RJ(3)
	XI(2)=SYM_OP(2,1)*RJ(1)+SYM_OP(2,2)*RJ(2)+SYM_OP(2,3)*RJ(3)
	XI(3)=SYM_OP(3,1)*RJ(1)+SYM_OP(3,2)*RJ(2)+SYM_OP(3,3)*RJ(3)
	
    DOT=Qx*XI(1)/aa_wz+Qy*XI(2)/aa_wz+Qz*XI(3)/cc_wz
    CPH=EXP(CI*DOT*2*PI)
    
    DD_temp=MATMUL(MATMUL(SYM_OP,D_34),TRANSPOSE(SYM_OP))
    Dsr_34=Dsr_34+DD_temp*CPH
ENDDO
!===================================================================
Dsr_13(:,:)=0.
D_13(2,2)=A2
D_13(3,3)=B2
D_13(2,3)=C2
D_13(3,2)=D2
!D_13=[0,   0,  0;...
!      0,  A2, C2;...
!      0,  D2, B2];
RJ=RI(3,:)-RI(1,:)
DO near_ii=1,3
    SELECT CASE  (near_ii)
        CASE (1) 
            SYM_OP=II_op
        CASE (2) 
            SYM_OP=C3_op
        CASE (3) 
            SYM_OP=TRANSPOSE(C3_op)
    END SELECT
    XI(1)=SYM_OP(1,1)*RJ(1)+SYM_OP(1,2)*RJ(2)+SYM_OP(1,3)*RJ(3)
	XI(2)=SYM_OP(2,1)*RJ(1)+SYM_OP(2,2)*RJ(2)+SYM_OP(2,3)*RJ(3)
	XI(3)=SYM_OP(3,1)*RJ(1)+SYM_OP(3,2)*RJ(2)+SYM_OP(3,3)*RJ(3)
	
    DOT=Qx*XI(1)/aa_wz+Qy*XI(2)/aa_wz+Qz*XI(3)/cc_wz
    CPH=EXP(CI*DOT*2*PI)
    
	DD_temp=MATMUL(MATMUL(SYM_OP,D_13),TRANSPOSE(SYM_OP))
    Dsr_13=Dsr_13+DD_temp*CPH
ENDDO
!-------------------------------------------------------------------
Dsr_13b(:,:)=0.
D_13(2,2)=A2
D_13(3,3)=B2
D_13(2,3)=-D2
D_13(3,2)=-C2
!D_13=[0,   0,  0;...
!      0,  A2, C2;...
!      0,  D2, B2];
RL=RI(11,:)-RI(3,:)
RJ=RI(3,:)-RI(1,:)
DO near_ii=1,3
    SELECT CASE  (near_ii)
        CASE (1) 
            SYM_OP=II_op
        CASE (2) 
            SYM_OP=C3_op
        CASE (3) 
            SYM_OP=TRANSPOSE(C3_op)
    END SELECT
    XI(1)=SYM_OP(1,1)*RJ(1)+SYM_OP(1,2)*RJ(2)+SYM_OP(1,3)*RJ(3)
	XI(2)=SYM_OP(2,1)*RJ(1)+SYM_OP(2,2)*RJ(2)+SYM_OP(2,3)*RJ(3)
	XI(3)=SYM_OP(3,1)*RJ(1)+SYM_OP(3,2)*RJ(2)+SYM_OP(3,3)*RJ(3)
	
    DOT=Qx*XI(1)/aa_wz+Qy*XI(2)/aa_wz+Qz*XI(3)/cc_wz
	DOT_RL=Qx*RL(1)/aa_wz+Qy*RL(2)/aa_wz+Qz*RL(3)/cc_wz
    CPH=EXP(CI*DOT*2*PI)
	CPH_RL=EXP(CI*DOT_RL*2*PI)
    
	DD_temp=MATMUL(MATMUL(SYM_OP,D_13),TRANSPOSE(SYM_OP))
    Dsr_13b=Dsr_13b+DD_temp*CPH*CPH_RL
ENDDO
!===================================================================
Dsr_14(:,:)=0.
D_14(2,2)=A3
D_14(3,3)=B3
D_14(2,3)=C3
D_14(3,2)=C3
!D_14=[0,   0,  0;...
!      0,  A3,  C3;...
!      0,  C3,  B3];
RJ=RI(4,:)-RI(1,:)
DO near_ii=1,3
    SELECT CASE  (near_ii)
        CASE (1) 
            SYM_OP=II_op
        CASE (2) 
            SYM_OP=C3_op
        CASE (3) 
            SYM_OP=TRANSPOSE(C3_op)
    END SELECT
    XI(1)=SYM_OP(1,1)*RJ(1)+SYM_OP(1,2)*RJ(2)+SYM_OP(1,3)*RJ(3)
	XI(2)=SYM_OP(2,1)*RJ(1)+SYM_OP(2,2)*RJ(2)+SYM_OP(2,3)*RJ(3)
	XI(3)=SYM_OP(3,1)*RJ(1)+SYM_OP(3,2)*RJ(2)+SYM_OP(3,3)*RJ(3)
	
    DOT=Qx*XI(1)/aa_wz+Qy*XI(2)/aa_wz+Qz*XI(3)/cc_wz
    CPH=EXP(CI*DOT*2*PI)
    
	DD_temp=MATMUL(MATMUL(SYM_OP,D_14),TRANSPOSE(SYM_OP))
    Dsr_14=Dsr_14+DD_temp*CPH
ENDDO
!-------------------------------------------------------------------
Dsr_23b(:,:)=0.
D_23b(2,2)=A3
D_23b(3,3)=B3
D_23b(2,3)=-C3
D_23b(3,2)=-C3
!D_23b=[0,   0,  0;...
!      0,  A3,  -C3;...
!      0, -C3,  B3];
RL=RI(11,:)-RI(3,:)
RJ=RI(3,:)-RI(2,:)
DO near_ii=1,3
    SELECT CASE  (near_ii)
        CASE (1) 
            SYM_OP=II_op
        CASE (2) 
            SYM_OP=C3_op
        CASE (3) 
            SYM_OP=TRANSPOSE(C3_op)
    END SELECT
    XI(1)=SYM_OP(1,1)*RJ(1)+SYM_OP(1,2)*RJ(2)+SYM_OP(1,3)*RJ(3)
	XI(2)=SYM_OP(2,1)*RJ(1)+SYM_OP(2,2)*RJ(2)+SYM_OP(2,3)*RJ(3)
	XI(3)=SYM_OP(3,1)*RJ(1)+SYM_OP(3,2)*RJ(2)+SYM_OP(3,3)*RJ(3)
	
    DOT=Qx*XI(1)/aa_wz+Qy*XI(2)/aa_wz+Qz*XI(3)/cc_wz
	DOT_RL=Qx*RL(1)/aa_wz+Qy*RL(2)/aa_wz+Qz*RL(3)/cc_wz
    CPH=EXP(CI*DOT*2*PI)
	CPH_RL=EXP(CI*DOT_RL*2*PI)
    
	DD_temp=MATMUL(MATMUL(SYM_OP,D_23b),TRANSPOSE(SYM_OP))
    Dsr_23b=Dsr_23b+DD_temp*CPH*CPH_RL
ENDDO
!===================================================================
Dsr_23(:,:)=0.
D_23(1,1)=A4
D_23(2,2)=B4
D_23(3,3)=C4
D_23(2,3)=D4
D_23(3,2)=E4
!D_23=[A4, 0,    0;...
!      0,  B4,  D4;...
!      0,  E4,  C4];
RJ=RI(3,:)-RI(2,:)
DO near_ii=1,3
    SELECT CASE  (near_ii)
        CASE (1) 
            SYM_OP=II_op
        CASE (2) 
            SYM_OP=C3_op
        CASE (3) 
            SYM_OP=TRANSPOSE(C3_op)
    END SELECT
    XI(1)=SYM_OP(1,1)*RJ(1)+SYM_OP(1,2)*RJ(2)+SYM_OP(1,3)*RJ(3)
	XI(2)=SYM_OP(2,1)*RJ(1)+SYM_OP(2,2)*RJ(2)+SYM_OP(2,3)*RJ(3)
	XI(3)=SYM_OP(3,1)*RJ(1)+SYM_OP(3,2)*RJ(2)+SYM_OP(3,3)*RJ(3)
	
    DOT=Qx*XI(1)/aa_wz+Qy*XI(2)/aa_wz+Qz*XI(3)/cc_wz
    CPH=EXP(CI*DOT*2*PI)
    
	DD_temp=MATMUL(MATMUL(SYM_OP,D_23),TRANSPOSE(SYM_OP))
    Dsr_23=Dsr_23+DD_temp*CPH
ENDDO
!-------------------------------------------------------------------
Dsr_14b(:,:)=0.
D_14b(1,1)=A4
D_14b(2,2)=B4
D_14b(3,3)=C4
D_14b(2,3)=-E4
D_14b(3,2)=-D4
!D_14b=[A4,   0,   0;...
!        0,  B4, -E4;...
!        0, -D4,  C4];
RL=RI(12,:)-RI(4,:)
RJ=RI(4,:)-RI(1,:)
DO near_ii=1,3
    SELECT CASE  (near_ii)
        CASE (1) 
            SYM_OP=II_op
        CASE (2) 
            SYM_OP=C3_op
        CASE (3) 
            SYM_OP=TRANSPOSE(C3_op)
    END SELECT
    XI(1)=SYM_OP(1,1)*RJ(1)+SYM_OP(1,2)*RJ(2)+SYM_OP(1,3)*RJ(3)
	XI(2)=SYM_OP(2,1)*RJ(1)+SYM_OP(2,2)*RJ(2)+SYM_OP(2,3)*RJ(3)
	XI(3)=SYM_OP(3,1)*RJ(1)+SYM_OP(3,2)*RJ(2)+SYM_OP(3,3)*RJ(3)
	
    DOT=Qx*XI(1)/aa_wz+Qy*XI(2)/aa_wz+Qz*XI(3)/cc_wz
    DOT_RL=Qx*RL(1)/aa_wz+Qy*RL(2)/aa_wz+Qz*RL(3)/cc_wz
    CPH=EXP(CI*DOT*2*PI)
	CPH_RL=EXP(CI*DOT_RL*2*PI)
    
	DD_temp=MATMUL(MATMUL(SYM_OP,D_14b),TRANSPOSE(SYM_OP))
    Dsr_14b=Dsr_14b+DD_temp*CPH*CPH_RL
ENDDO
!===================================================================
Dsr_24(:,:)=0.
D_24(2,2)=A5
D_24(3,3)=B5
D_24(2,3)=C5
D_24(3,2)=D5
!D_24=[0,   0,  0;...
!      0,  A5, C5;...
!      0,  D5, B5];
RJ=RI(4,:)-RI(2,:)
DO near_ii=1,3
    SELECT CASE  (near_ii)
        CASE (1) 
            SYM_OP=II_op
        CASE (2) 
            SYM_OP=C3_op
        CASE (3) 
            SYM_OP=TRANSPOSE(C3_op)
    END SELECT
    XI(1)=SYM_OP(1,1)*RJ(1)+SYM_OP(1,2)*RJ(2)+SYM_OP(1,3)*RJ(3)
	XI(2)=SYM_OP(2,1)*RJ(1)+SYM_OP(2,2)*RJ(2)+SYM_OP(2,3)*RJ(3)
	XI(3)=SYM_OP(3,1)*RJ(1)+SYM_OP(3,2)*RJ(2)+SYM_OP(3,3)*RJ(3)
	
    DOT=Qx*XI(1)/aa_wz+Qy*XI(2)/aa_wz+Qz*XI(3)/cc_wz
    CPH=EXP(CI*DOT*2*PI)
    
	DD_temp=MATMUL(MATMUL(SYM_OP,D_24),TRANSPOSE(SYM_OP))
    Dsr_24=Dsr_24+DD_temp*CPH
	!WRITE(6,7001) REAL(TRANSPOSE(Dsr_24))
ENDDO
!-------------------------------------------------------------------
Dsr_24b(:,:)=0.
D_24(2,2)=A5
D_24(3,3)=B5
D_24(2,3)=-D5
D_24(3,2)=-C5
!D_24=[0,   0,  0;...
!      0,  A5, -D5;...
!      0, -C5,  B5];
RL=RI(12,:)-RI(4,:)
RJ=RI(4,:)-RI(2,:)
DO near_ii=1,3
    SELECT CASE  (near_ii)
        CASE (1) 
            SYM_OP=II_op
        CASE (2) 
            SYM_OP=C3_op
        CASE (3) 
            SYM_OP=TRANSPOSE(C3_op)
    END SELECT
    XI(1)=SYM_OP(1,1)*RJ(1)+SYM_OP(1,2)*RJ(2)+SYM_OP(1,3)*RJ(3)
	XI(2)=SYM_OP(2,1)*RJ(1)+SYM_OP(2,2)*RJ(2)+SYM_OP(2,3)*RJ(3)
	XI(3)=SYM_OP(3,1)*RJ(1)+SYM_OP(3,2)*RJ(2)+SYM_OP(3,3)*RJ(3)
	
    DOT=Qx*XI(1)/aa_wz+Qy*XI(2)/aa_wz+Qz*XI(3)/cc_wz
    DOT_RL=Qx*RL(1)/aa_wz+Qy*RL(2)/aa_wz+Qz*RL(3)/cc_wz
    CPH=EXP(CI*DOT*2*PI)
	CPH_RL=EXP(CI*DOT_RL*2*PI)
    
	DD_temp=MATMUL(MATMUL(SYM_OP,D_24),TRANSPOSE(SYM_OP))
    Dsr_24b=Dsr_24b+DD_temp*CPH*CPH_RL
ENDDO
!===================================================================
Dsr_12b(:,:)=0.
D_12b(3,3)=A6
!D_12b=[0,   0,  0;...
!       0,   0,  0;...
!       0,   0, A6];
RL=RI(10,:)-RI(2,:)
RJ=RI(2,:)-RI(1,:)
DO near_ii=1,1
    SELECT CASE  (near_ii)
        CASE (1) 
            SYM_OP=II_op
    END SELECT
    XI(1)=SYM_OP(1,1)*RJ(1)+SYM_OP(1,2)*RJ(2)+SYM_OP(1,3)*RJ(3)
	XI(2)=SYM_OP(2,1)*RJ(1)+SYM_OP(2,2)*RJ(2)+SYM_OP(2,3)*RJ(3)
	XI(3)=SYM_OP(3,1)*RJ(1)+SYM_OP(3,2)*RJ(2)+SYM_OP(3,3)*RJ(3)
	
    DOT=Qx*XI(1)/aa_wz+Qy*XI(2)/aa_wz+Qz*XI(3)/cc_wz
    DOT_RL=Qx*RL(1)/aa_wz+Qy*RL(2)/aa_wz+Qz*RL(3)/cc_wz
    CPH=EXP(CI*DOT*2*PI)
	CPH_RL=EXP(CI*DOT_RL*2*PI)
    
	DD_temp=MATMUL(MATMUL(SYM_OP,D_12b),TRANSPOSE(SYM_OP))
    Dsr_12b=Dsr_12b+DD_temp*CPH*CPH_RL
ENDDO
!-------------------------------------------------------------------
Dsr_34b(:,:)=0.
D_12b(3,3)=A6
!D_12b=[0,   0,  0;...
!       0,   0,  0;...
!       0,   0, A6];
RL=RI(12,:)-RI(4,:)
RJ=RI(4,:)-RI(3,:)
DO near_ii=1,1
    SELECT CASE  (near_ii)
        CASE (1) 
            SYM_OP=II_op
    END SELECT
    XI(1)=SYM_OP(1,1)*RJ(1)+SYM_OP(1,2)*RJ(2)+SYM_OP(1,3)*RJ(3)
	XI(2)=SYM_OP(2,1)*RJ(1)+SYM_OP(2,2)*RJ(2)+SYM_OP(2,3)*RJ(3)
	XI(3)=SYM_OP(3,1)*RJ(1)+SYM_OP(3,2)*RJ(2)+SYM_OP(3,3)*RJ(3)
	
    DOT=Qx*XI(1)/aa_wz+Qy*XI(2)/aa_wz+Qz*XI(3)/cc_wz
    DOT_RL=Qx*RL(1)/aa_wz+Qy*RL(2)/aa_wz+Qz*RL(3)/cc_wz
    CPH=EXP(CI*DOT*2*PI)
	CPH_RL=EXP(CI*DOT_RL*2*PI)
    
	DD_temp=MATMUL(MATMUL(SYM_OP,D_12b),TRANSPOSE(SYM_OP))
    Dsr_34b=Dsr_34b+DD_temp*CPH*CPH_RL
ENDDO
!===================================================================
Dsr_11p(:,:)=0.
D_11p(1,1)=A7
D_11p(2,2)=B7
D_11p(3,3)=C7
D_11p(1,2)=D7
D_11p(1,3)=E7
D_11p(2,3)=F7
D_11p(2,1)=-D7
D_11p(3,1)=-E7
D_11p(3,2)=F7
!D_11p=[ A7, D7, E7;...
!       -D7, B7, F7;...
!       -E7, F7, C7];
RL=RI(13,:)-RI(1,:)
RJ=RI(1,:)-RI(1,:)
DO near_ii=1,6
    SELECT CASE (near_ii)
        CASE (1) 
            SYM_OP=II_op
        CASE (2) 
            SYM_OP=sigma_v2
        CASE (3) 
            SYM_OP=C3_op
        CASE (4)
            SYM_OP=sigma_v1
        CASE (5) 
            SYM_OP=TRANSPOSE(C3_op)
        CASE (6) 
            SYM_OP=sigma_v3
    END SELECT
    XI(1)=SYM_OP(1,1)*RJ(1)+SYM_OP(1,2)*RJ(2)+SYM_OP(1,3)*RJ(3)
	XI(2)=SYM_OP(2,1)*RJ(1)+SYM_OP(2,2)*RJ(2)+SYM_OP(2,3)*RJ(3)
	XI(3)=SYM_OP(3,1)*RJ(1)+SYM_OP(3,2)*RJ(2)+SYM_OP(3,3)*RJ(3)
	
	XL(1)=SYM_OP(1,1)*RL(1)+SYM_OP(1,2)*RL(2)+SYM_OP(1,3)*RL(3)
	XL(2)=SYM_OP(2,1)*RL(1)+SYM_OP(2,2)*RL(2)+SYM_OP(2,3)*RL(3)
	XL(3)=SYM_OP(3,1)*RL(1)+SYM_OP(3,2)*RL(2)+SYM_OP(3,3)*RL(3)
	
    DOT=Qx*XI(1)/aa_wz+Qy*XI(2)/aa_wz+Qz*XI(3)/cc_wz
    DOT_RL=Qx*XL(1)/aa_wz+Qy*XL(2)/aa_wz+Qz*XL(3)/cc_wz
    CPH=EXP(CI*DOT*2*PI)
	CPH_RL=EXP(CI*DOT_RL*2*PI)
    
 	DD_temp=MATMUL(MATMUL(SYM_OP,D_11p),TRANSPOSE(SYM_OP))
    Dsr_11p=Dsr_11p+DD_temp*CPH*CPH_RL
	!write(6,7001) REAL(TRANSPOSE(Dsr_11p))
ENDDO
!-------------------------------------------------------------------
Dsr_22p(:,:)=0.
D_22p(1,1)=A8
D_22p(2,2)=B8
D_22p(3,3)=C8
D_22p(1,2)=D8
D_22p(1,3)=E8
D_22p(2,3)=F8
D_22p(2,1)=-D8
D_22p(3,1)=-E8
D_22p(3,2)=F8
!D_22p=[ A8, D8, E8;...
!       -D8, B8, F8;...
!       -E8, F8, C8];
RL=RI(16,:)-RI(2,:)
RJ=RI(2,:)-RI(2,:)
DO near_ii=1,6
    SELECT CASE (near_ii)
        CASE (1) 
            SYM_OP=II_op
        CASE (2) 
            SYM_OP=sigma_v2
        CASE (3) 
            SYM_OP=C3_op
        CASE (4) 
            SYM_OP=sigma_v1
        CASE (5) 
            SYM_OP=TRANSPOSE(C3_op)
        CASE (6) 
            SYM_OP=sigma_v3
    END SELECT
    XI(1)=SYM_OP(1,1)*RJ(1)+SYM_OP(1,2)*RJ(2)+SYM_OP(1,3)*RJ(3)
	XI(2)=SYM_OP(2,1)*RJ(1)+SYM_OP(2,2)*RJ(2)+SYM_OP(2,3)*RJ(3)
	XI(3)=SYM_OP(3,1)*RJ(1)+SYM_OP(3,2)*RJ(2)+SYM_OP(3,3)*RJ(3)
	
	XL(1)=SYM_OP(1,1)*RL(1)+SYM_OP(1,2)*RL(2)+SYM_OP(1,3)*RL(3)
	XL(2)=SYM_OP(2,1)*RL(1)+SYM_OP(2,2)*RL(2)+SYM_OP(2,3)*RL(3)
	XL(3)=SYM_OP(3,1)*RL(1)+SYM_OP(3,2)*RL(2)+SYM_OP(3,3)*RL(3)
	
    DOT=Qx*XI(1)/aa_wz+Qy*XI(2)/aa_wz+Qz*XI(3)/cc_wz
    DOT_RL=Qx*XL(1)/aa_wz+Qy*XL(2)/aa_wz+Qz*XL(3)/cc_wz
    CPH=EXP(CI*DOT*2*PI)
	CPH_RL=EXP(CI*DOT_RL*2*PI)
	
    
 	DD_temp=MATMUL(MATMUL(SYM_OP,D_22p),TRANSPOSE(SYM_OP))
    Dsr_22p=Dsr_22p+DD_temp*CPH*CPH_RL
ENDDO
!-------------------------------------------------------------------
Dsr_33p(:,:)=0.
D_33p(1,1)=A7
D_33p(2,2)=B7
D_33p(3,3)=C7
D_33p(1,2)=-D7
D_33p(1,3)=E7
D_33p(2,3)=-F7
D_33p(2,1)=D7
D_33p(3,1)=-E7
D_33p(3,2)=-F7
!D_33p=[ A7, -D7,   E7;...
!        D7,  B7,  -F7;...
!       -E7, -F7,   C7];
RL=RI(19,:)-RI(3,:)
RJ=RI(3,:)-RI(3,:)
DO near_ii=1,6
    SELECT CASE (near_ii)
        CASE (1) 
            SYM_OP=II_op
        CASE (2) 
            SYM_OP=sigma_v2
        CASE (3) 
            SYM_OP=C3_op
        CASE (4)
            SYM_OP=sigma_v1
        CASE (5) 
            SYM_OP=TRANSPOSE(C3_op)
        CASE (6) 
            SYM_OP=sigma_v3
    END SELECT
     XI(1)=SYM_OP(1,1)*RJ(1)+SYM_OP(1,2)*RJ(2)+SYM_OP(1,3)*RJ(3)
	XI(2)=SYM_OP(2,1)*RJ(1)+SYM_OP(2,2)*RJ(2)+SYM_OP(2,3)*RJ(3)
	XI(3)=SYM_OP(3,1)*RJ(1)+SYM_OP(3,2)*RJ(2)+SYM_OP(3,3)*RJ(3)
	
	XL(1)=SYM_OP(1,1)*RL(1)+SYM_OP(1,2)*RL(2)+SYM_OP(1,3)*RL(3)
	XL(2)=SYM_OP(2,1)*RL(1)+SYM_OP(2,2)*RL(2)+SYM_OP(2,3)*RL(3)
	XL(3)=SYM_OP(3,1)*RL(1)+SYM_OP(3,2)*RL(2)+SYM_OP(3,3)*RL(3)
	
    DOT=Qx*XI(1)/aa_wz+Qy*XI(2)/aa_wz+Qz*XI(3)/cc_wz
    DOT_RL=Qx*XL(1)/aa_wz+Qy*XL(2)/aa_wz+Qz*XL(3)/cc_wz
    CPH=EXP(CI*DOT*2*PI)
	CPH_RL=EXP(CI*DOT_RL*2*PI)
	
    
 	DD_temp=MATMUL(MATMUL(SYM_OP,D_33p),TRANSPOSE(SYM_OP))
    Dsr_33p=Dsr_33p+DD_temp*CPH*CPH_RL
ENDDO
!-------------------------------------------------------------------
Dsr_44p(:,:)=0.
D_44p(1,1)=A8
D_44p(2,2)=B8
D_44p(3,3)=C8
D_44p(1,2)=-D8
D_44p(1,3)=E8
D_44p(2,3)=-F8
D_44p(2,1)=D8
D_44p(3,1)=-E8
D_44p(3,2)=-F8
!D_44p=[ A8, -D8,  E8;...
!        D8,  B8, -F8;...
!       -E8, -F8,  C8];
RL=RI(22,:)-RI(4,:)
RJ=RI(4,:)-RI(4,:)
DO near_ii=1,6
    SELECT CASE (near_ii)
        CASE (1) 
            SYM_OP=II_op
        CASE (2) 
            SYM_OP=sigma_v2
        CASE (3) 
            SYM_OP=C3_op
        CASE (4) 
            SYM_OP=sigma_v1
        CASE (5) 
            SYM_OP=TRANSPOSE(C3_op)
        CASE (6) 
            SYM_OP=sigma_v3
    END SELECT
    XI(1)=SYM_OP(1,1)*RJ(1)+SYM_OP(1,2)*RJ(2)+SYM_OP(1,3)*RJ(3)
	XI(2)=SYM_OP(2,1)*RJ(1)+SYM_OP(2,2)*RJ(2)+SYM_OP(2,3)*RJ(3)
	XI(3)=SYM_OP(3,1)*RJ(1)+SYM_OP(3,2)*RJ(2)+SYM_OP(3,3)*RJ(3)
	
	XL(1)=SYM_OP(1,1)*RL(1)+SYM_OP(1,2)*RL(2)+SYM_OP(1,3)*RL(3)
	XL(2)=SYM_OP(2,1)*RL(1)+SYM_OP(2,2)*RL(2)+SYM_OP(2,3)*RL(3)
	XL(3)=SYM_OP(3,1)*RL(1)+SYM_OP(3,2)*RL(2)+SYM_OP(3,3)*RL(3)
	
    DOT=Qx*XI(1)/aa_wz+Qy*XI(2)/aa_wz+Qz*XI(3)/cc_wz
    DOT_RL=Qx*XL(1)/aa_wz+Qy*XL(2)/aa_wz+Qz*XL(3)/cc_wz
    CPH=EXP(CI*DOT*2*PI)
	CPH_RL=EXP(CI*DOT_RL*2*PI)
    
 	DD_temp=MATMUL(MATMUL(SYM_OP,D_44p),TRANSPOSE(SYM_OP))
    Dsr_44p=Dsr_44p+DD_temp*CPH*CPH_RL
ENDDO
!===================================================================
Dsr_14_COP(:,:)=0.
D_14_COP(1,1)=A9
D_14_COP(2,2)=B9
D_14_COP(3,3)=C9
D_14_COP(1,2)=D9
D_14_COP(1,3)=E9
D_14_COP(2,3)=F9
D_14_COP(2,1)=D9
D_14_COP(3,1)=E9
D_14_COP(3,2)=F9
!D_14cop=[ A9,  D9,  E9;...
!          D9,  B9,  F9;...
!          E9,  F9,  C9];
RL=RI(37,:)-RI(4,:)
RJ=RI(4,:)-RI(1,:)
DO near_ii=1,6
    SELECT CASE (near_ii)
        CASE (1) 
            SYM_OP=II_op
        CASE (2) 
            SYM_OP=sigma_v2
        CASE (3) 
            SYM_OP=C3_op
        CASE (4) 
            SYM_OP=sigma_v1
        CASE (5) 
            SYM_OP=TRANSPOSE(C3_op)
        CASE (6) 
            SYM_OP=sigma_v3
    END SELECT
    XI(1)=SYM_OP(1,1)*RJ(1)+SYM_OP(1,2)*RJ(2)+SYM_OP(1,3)*RJ(3)
	XI(2)=SYM_OP(2,1)*RJ(1)+SYM_OP(2,2)*RJ(2)+SYM_OP(2,3)*RJ(3)
	XI(3)=SYM_OP(3,1)*RJ(1)+SYM_OP(3,2)*RJ(2)+SYM_OP(3,3)*RJ(3)
	
	XL(1)=SYM_OP(1,1)*RL(1)+SYM_OP(1,2)*RL(2)+SYM_OP(1,3)*RL(3)
	XL(2)=SYM_OP(2,1)*RL(1)+SYM_OP(2,2)*RL(2)+SYM_OP(2,3)*RL(3)
	XL(3)=SYM_OP(3,1)*RL(1)+SYM_OP(3,2)*RL(2)+SYM_OP(3,3)*RL(3)
	
    DOT=Qx*XI(1)/aa_wz+Qy*XI(2)/aa_wz+Qz*XI(3)/cc_wz
    DOT_RL=Qx*XL(1)/aa_wz+Qy*XL(2)/aa_wz+Qz*XL(3)/cc_wz
    CPH=EXP(CI*DOT*2*PI)
	CPH_RL=EXP(CI*DOT_RL*2*PI)
    
 	DD_temp=MATMUL(MATMUL(SYM_OP,D_14_COP),TRANSPOSE(SYM_OP))
    Dsr_14_COP=Dsr_14_COP+DD_temp*CPH*CPH_RL
ENDDO
!-------------------------------------------------------------------
Dsr_23_COP(:,:)=0.
D_23_COP(1,1)=A9
D_23_COP(2,2)=B9
D_23_COP(3,3)=C9
D_23_COP(1,2)=D9
D_23_COP(1,3)=-E9
D_23_COP(2,3)=-F9
D_23_COP(2,1)=D9
D_23_COP(3,1)=-E9
D_23_COP(3,2)=-F9
!D_23cop=[ A9,  D9,  -E9;...
!          D9,  B9,  -F9;...
!         -E9, -F9,   C9];
RL=RI(38,:)-RI(3,:)
RJ=RI(3,:)-RI(2,:)
DO near_ii=1,6
    SELECT CASE (near_ii)
        CASE (1) 
            SYM_OP=II_op
        CASE (2) 
            SYM_OP=sigma_v2
        CASE (3) 
            SYM_OP=C3_op
        CASE (4) 
            SYM_OP=sigma_v1
        CASE (5) 
            SYM_OP=TRANSPOSE(C3_op)
        CASE (6) 
            SYM_OP=sigma_v3
    END SELECT
    XI(1)=SYM_OP(1,1)*RJ(1)+SYM_OP(1,2)*RJ(2)+SYM_OP(1,3)*RJ(3)
	XI(2)=SYM_OP(2,1)*RJ(1)+SYM_OP(2,2)*RJ(2)+SYM_OP(2,3)*RJ(3)
	XI(3)=SYM_OP(3,1)*RJ(1)+SYM_OP(3,2)*RJ(2)+SYM_OP(3,3)*RJ(3)
	
	XL(1)=SYM_OP(1,1)*RL(1)+SYM_OP(1,2)*RL(2)+SYM_OP(1,3)*RL(3)
	XL(2)=SYM_OP(2,1)*RL(1)+SYM_OP(2,2)*RL(2)+SYM_OP(2,3)*RL(3)
	XL(3)=SYM_OP(3,1)*RL(1)+SYM_OP(3,2)*RL(2)+SYM_OP(3,3)*RL(3)
	
    DOT=Qx*XI(1)/aa_wz+Qy*XI(2)/aa_wz+Qz*XI(3)/cc_wz
    DOT_RL=Qx*XL(1)/aa_wz+Qy*XL(2)/aa_wz+Qz*XL(3)/cc_wz
    CPH=EXP(CI*DOT*2*PI)
	CPH_RL=EXP(CI*DOT_RL*2*PI)
    
 	DD_temp=MATMUL(MATMUL(SYM_OP,D_23_COP),TRANSPOSE(SYM_OP))
    Dsr_23_COP=Dsr_23_COP+DD_temp*CPH*CPH_RL
ENDDO
!===================================================================
DSR(1,1,:,:)=Dsr_11p
DSR(1,2,:,:)=Dsr_12+Dsr_12b
DSR(1,3,:,:)=Dsr_13+Dsr_13b
DSR(1,4,:,:)=Dsr_14+Dsr_14b+Dsr_14_COP

DSR(2,2,:,:)=Dsr_22p
DSR(2,3,:,:)=Dsr_23+Dsr_23b+Dsr_23_COP
DSR(2,4,:,:)=Dsr_24+Dsr_24b

DSR(3,3,:,:)=Dsr_33p
DSR(3,4,:,:)=Dsr_34+Dsr_34b

DSR(4,4,:,:)=Dsr_44p

!===================================================================
! complete Dynamic matrix using translational invariance------------
DO K1=2,4
DO K2=1,K1-1
DO ii=1,3
    DO jj=1,3
        DSR(K1,K2,ii,jj)=CONJG(DSR(K2,K1,jj,ii))    
    ENDDO
ENDDO
ENDDO
ENDDO
!===================================================================
! define self-interaction for Q1=0;
DDR(:,:,:)=0.
IF (Q1.LT.0.00001) THEN
    DO ion_ii=1,4
        DO ii=1,3
            DO jj=1,3
                DD=0.
                DO ion_jj=1,4
                    DD=DD-DSR(ion_ii,ion_jj,ii,jj)
                ENDDO
                DDR(ion_ii,ii,jj)=DD
            ENDDO
        ENDDO
    ENDDO
ENDIF



!WRITE(6,7001) REAL(TRANSPOSE(Dsr_24))
!WRITE(6,7001) REAL(TRANSPOSE(Dsr_24b))


7001 FORMAT(3F10.4)




END SUBROUTINE



















