C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
C     This file is part of FEDEM - https://openfedem.org
C
      SUBROUTINE TMRF35(EM,X,Y,Z,THK,RHO)
C
C***********************************************************************
C
C     TMRF35 DETERMINES A LUMPED MASS MATRIX  EM  FOR THE TRIANGULAR
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
      DOUBLE PRECISION  EM(9),X(3),Y(3),Z(3),THK,RHO
      DOUBLE PRECISION  COSG,SING,SL21,SL31
      DOUBLE PRECISION  X21,X31,Y21,Y31,Z21,Z31,XL(3),YL(3)
      DOUBLE PRECISION  AREA,XTP,YTP,H,B,B1,B2,S1,S2
      DOUBLE PRECISION  IXX,IYY1,IYY2,IYY,IZZ1,IZZ2,IZZ
      DOUBLE PRECISION  IX,IY,IZ,M,GAMMA,TR(3,3)
C
C ---------------------------- Local coordinates -------------------
C
      X21 = X(2)-X(1)
      Y21 = Y(2)-Y(1)
      Z21 = Z(2)-Z(1)
      X31 = X(3)-X(1)
      Y31 = Y(3)-Y(1)
      Z31 = Z(3)-Z(1)
      SL21 = DSQRT(X21*X21 + Y21*Y21 + Z21*Z21)
      SL31 = DSQRT(X31*X31 + Y31*Y31 + Z31*Z31)
      COSG = (X31*X21 + Y31*Y21 + Z31*Z21)/(SL31*SL21)
      IF (DABS(COSG) .GE. 1.0D0) THEN
         SING = 0.0D0
      ELSE
         SING = DSQRT(1.0D0 - COSG*COSG)
      END IF
C
      XL(1) = 0.
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
      IF (AREA .LE. 0.0D0) GOTO 80
C
C ---------------------------- Koordinates for TP ------------------
C
      XTP = ((XL(3)+((XL(2)-XL(3))/2.)))*(2./3.)
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
C ---------------------------- V1 ----------------------------------
      TR(1,1) = X(2)-X(1)
      TR(2,1) = Y(2)-Y(1)
      TR(3,1) = Z(2)-Z(1)
      S1 = DSQRT(TR(1,1)**2+TR(2,1)**2+TR(3,1)**2)
      TR(1,1) = TR(1,1)/S1
      TR(2,1) = TR(2,1)/S1
      TR(3,1) = TR(3,1)/S1
C ---------------------------- V3 ----------------------------------
      TR(1,2) = X(3)-X(1)
      TR(2,2) = Y(3)-Y(1)
      TR(3,2) = Z(3)-Z(1)
      TR(1,3) = (TR(2,1)*TR(3,2)-TR(2,2)*TR(3,1))
      TR(2,3) = (TR(1,2)*TR(3,1)-TR(1,1)*TR(3,2))
      TR(3,3) = (TR(1,1)*TR(2,2)-TR(1,2)*TR(2,1))
C
      S1 = DSQRT(TR(1,3)**2+TR(2,3)**2+TR(3,3)**2)
      TR(1,3) = TR(1,3)/S1
      TR(2,3) = TR(2,3)/S1
      TR(3,3) = TR(3,3)/S1
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
      EM(1) = M
      EM(2) = EM(1)
      EM(4) = EM(1)
      EM(5) = EM(1)
      EM(7) = EM(1)
      EM(8) = EM(1)
C
      EM(3) = GAMMA*IZ
      EM(6) = EM(3)
      EM(9) = EM(3)
C ---------------------------- Return ------------------------------
   80 RETURN
      END
