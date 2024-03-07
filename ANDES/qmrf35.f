C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
C     This file is part of FEDEM - https://openfedem.org
C
      SUBROUTINE QMRF35(EM,X,Y,Z,THK,RHO)
C
C***********************************************************************
C
C     QMRF35 DETERMINES A LUMPED MASS MATRIX  EM  FOR THE QUADRILATERAL
C     MEMBRANE ELEMENT WITH ROTATIONAL DEGREES OF FREEDOM OF P.G.BERGAN
C     AND C.A.FELIPPA.
C
C     PROGRAMMED BY :  KETIL AAMNES
C     DATE/VERSION  :  081085 / 1.0
C
C***********************************************************************
C
      IMPLICIT NONE
C
      DOUBLE PRECISION  EM(24),EMM(2,9),X(4),Y(4),Z(4),THK,RHO
      DOUBLE PRECISION  AUX1,COSG,SING,SL21,SL31
      DOUBLE PRECISION  X1,X2,X3,Y1,Y2,Y3,Z1,Z2,Z3,XL(4),YL(4)
      DOUBLE PRECISION  AREA,XTP,YTP,H,B,B1,B2,S1,S2,SL
      DOUBLE PRECISION  IXX,IYY1,IYY2,IYY,IZZ1,IZZ2,IZZ
      DOUBLE PRECISION  IX,IY,IZ,M,GAMMA,TR(3,3)
      INTEGER           I
C
C ---------------------------- Local coordinates -------------------
C
      DO 10 I=1,2
C
C
      X1 = 0.
      Y1 = 0.
      Z1 = 0.
C
      X2 = X(I+1)-X(1)
      Y2 = Y(I+1)-Y(1)
      Z2 = Z(I+1)-Z(1)
C
      X3 = X(I+2)-X(1)
      Y3 = Y(I+2)-Y(1)
      Z3 = Z(I+2)-Z(1)
C
      SL21 = DSQRT(X2*X2 + Y2*Y2 + Z2*Z2)
      SL31 = DSQRT(X3*X3 + Y3*Y3 + Z3*Z3)
C
      COSG = (X3*X2 + Y3*Y2 + Z3*Z2)/(SL31*SL21)
      AUX1 = 1.-COSG*COSG
C
      IF (AUX1) 40,40,50
   40 SING = 0.
      GO TO 60
   50 SING = DSQRT(AUX1)
C
   60 XL(1) = 0.
      XL(2) = SL21
      XL(3) = SL31*COSG
      YL(1) = 0.
      YL(2) = 0.
      YL(3) = SL31*SING
C
C ---------------------------- Element area ------------------------
C
      AREA = (XL(2)*YL(3))/2.
C
      IF (AREA) 80,80,70
C
C ---------------------------- Koordinates for TP ------------------
C
   70 XTP = ((XL(3)+((XL(2)-XL(3))/2.)))*(2./3.)
      YTP = (YL(3))/3.
C
C ------------------------------------------------------------------
C
      H = YL(3)
      B = XL(2)
C
      B1 = XL(3)
      B2 = XL(2)-XL(3)
C
      S1 = DSQRT( (((2./3.)*B1)**2) + (((1./3.)*H)**2) )
      S2 = DSQRT( (((2./3.)*B2)**2) + (((1./3.)*H)**2) )
C
C ---------------------------- Moment of inertia about the X-aksis -
C
      IXX = (RHO*THK*B*H) * (((H**2)/36.) + ((THK**2)/24.))
C
C ---------------------------- Moment of inertia about the Y-aksis -
C
      IYY1 = (RHO*THK*B1*H) * (((B1**2)/36.) + ((THK**2)/24.))
      IYY1 = IYY1 + (((XTP-((2./3.)*B1))**2)*(((B1*H)/2.)*RHO*THK))
C
      IYY2 = (RHO*THK*B2*H) * (((B2**2)/36.) + ((THK**2)/24.))
      IYY2 = IYY2 + (((B1-XTP+((1./3.)*B2))**2)*(((B2*H)/2.)*RHO*THK))
C
      IYY = IYY1 + IYY2
C
C ---------------------------- Moment of inertia about the Z-aksis -
C
      IZZ1 = (RHO*THK*H*B1) * (((B1**2)/4.)+((H**2)/12.))
      IZZ1 = IZZ1 - ((S1**2)*RHO*THK*(B1*H)/2.)
      IZZ1 = IZZ1 + (((XTP-((2./3.)*B1))**2)*((B1*H)/2.)*RHO*THK)
C
      IZZ2 = (RHO*THK*H*B2) * (((B2**2)/4.)+((H**2)/12.))
      IZZ2 = IZZ2 - ((S2**2)*RHO*THK*(B2*H)/2.)
      IZZ2 = IZZ2 + (((B1+((1./3.)*B2)-XTP)**2)*((B2*H)/2.)*RHO*THK)
C
      IZZ = IZZ1 + IZZ2
C
C ---------------------------- The 3 X 3 transformation matrix -----
C
C ---------------------------- V1 ----------------------------------
      TR(1,1) = X2-X1
      TR(2,1) = Y2-Y1
      TR(3,1) = Z2-Z1
      SL = DSQRT(TR(1,1)**2+TR(2,1)**2+TR(3,1)**2)
      TR(1,1) = TR(1,1)/SL
      TR(2,1) = TR(2,1)/SL
      TR(3,1) = TR(3,1)/SL
C ---------------------------- V3 ----------------------------------
      TR(1,2) = X3-X1
      TR(2,2) = Y3-Y1
      TR(3,2) = Z3-Z1
      TR(1,3) = (TR(2,1)*TR(3,2)-TR(2,2)*TR(3,1))
      TR(2,3) = (TR(1,2)*TR(3,1)-TR(1,1)*TR(3,2))
      TR(3,3) = (TR(1,1)*TR(2,2)-TR(1,2)*TR(2,1))
C
      SL = DSQRT(TR(1,3)**2+TR(2,3)**2+TR(3,3)**2)
      TR(1,3) = TR(1,3)/SL
      TR(2,3) = TR(2,3)/SL
      TR(3,3) = TR(3,3)/SL
C ---------------------------- V2 ----------------------------------
      TR(1,2) = TR(2,3)*TR(3,1)-TR(2,1)*TR(3,3)
      TR(2,2) = TR(1,1)*TR(3,3)-TR(1,3)*TR(3,1)
      TR(3,2) = TR(1,3)*TR(2,1)-TR(1,1)*TR(2,3)
C
C ---------------------------- Transform to global
C                              moments of inertia ------------------
      IX = ABS(TR(1,1)*IXX+TR(1,2)*IYY+TR(1,3)*IZZ)
      IY = ABS(TR(2,1)*IXX+TR(2,2)*IYY+TR(2,3)*IZZ)
      IZ = ABS(TR(3,1)*IXX+TR(3,2)*IYY+TR(3,3)*IZZ)
C
C ---------------------------- The elements of EM ------------------
C
      M = (RHO*AREA*THK)/3.
      GAMMA = 1./20.
C
      EMM(I,1) = M
      EMM(I,2) = M
      EMM(I,4) = M
      EMM(I,5) = M
      EMM(I,7) = M
      EMM(I,8) = M
C
      EMM(I,3) = GAMMA*IZ
      EMM(I,6) = EMM(I,3)
      EMM(I,9) = EMM(I,3)
C
   10 CONTINUE
C
      EM( 1) = EMM(1,1) + EMM(2,1)
      EM( 2) = EMM(1,2) + EMM(2,2)
      EM( 6) = EMM(1,3) + EMM(2,3)
C
      EM( 7) = EMM(1,4)
      EM( 8) = EMM(1,5)
      EM(12) = EMM(1,6)
C
      EM(13) = EMM(1,7) + EMM(2,4)
      EM(14) = EMM(1,8) + EMM(2,5)
      EM(18) = EMM(1,9) + EMM(2,6)
C
      EM(19) = EMM(2,7)
      EM(20) = EMM(2,8)
      EM(24) = EMM(2,9)
C
C---- Fordeler masser og treghetsmomenter likt pa de 4 knutepunktene
C
      EM(1) = (EM(1) + EM(7) + EM(13) + EM(19))/4.0
      EM(2) = (EM(2) + EM(8) + EM(14) + EM(20))/4.0
      EM(6) = (EM(6) + EM(12) + EM(18) + EM(24))/4.0
C
      EM(7)  = EM(1)
      EM(13) = EM(1)
      EM(19) = EM(1)
C
      EM(8)  = EM(2)
      EM(14) = EM(2)
      EM(20) = EM(2)
C
      EM(12) = EM(6)
      EM(18) = EM(6)
      EM(24) = EM(6)
C
C ---------------------------- Return ------------------------------
   80 RETURN
C
      END
