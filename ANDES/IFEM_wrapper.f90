!! $Id$
!!==============================================================================
!> @file IFEM_wrapper.f90
!> @brief Wrapper on the Andes shell elements from FEDEM (https://openfedem.org)
!>
!> @author Knut Morten Okstad / SINTEF
!>
!> $date Feb 20 2024

!> @brief Calculates the stiffness matrix for a 3-noded shell element.
subroutine IFEM_ANDES3 (iEL, X0, Thick, Emod, Rny, Rho, EK, IERR)

  use KindModule        , only : dp
  use Andes3ShellModule , only : Andes3shell_stiffmat
  use ManipMatrixModule , only : trans3p

  implicit none

  integer , parameter   :: nelnod = 3, neldof = 6*nelnod
  integer , intent(in)  :: iEL
  real(dp), intent(in)  :: X0(3,nelnod), Thick, Emod, Rny, Rho
  real(dp), intent(out) :: ek(neldof,neldof)
  integer , intent(out) :: ierr

  !! Local variables
  integer , parameter :: ltype  = 2, lpu = 6
  real(dp), parameter :: alpha  = 1.5_dp
  real(dp), parameter :: alphaH = 0.5_dp
  real(dp), parameter :: thetaMaterial = 0.0_dp

  integer  :: i, j
  real(dp) :: X21, Y21, Z21, X31, Y31, Z31, SL21, SL31, COSG, SING
  real(dp) :: Kmat(neldof,neldof), Cmat(6,6)
  real(dp) :: Xl(nelnod), Yl(nelnod), Trel(3,4), Tel(3,3), TelT(3,3)


  !! --- Logic section ---

  ek = 0.0_dp

  !! Set up the constitutive matrix

  Cmat = 0.0_dp
  Cmat(1,1) = Emod*Thick / (1.0_dp - Rny*Rny)
  Cmat(1,2) = Rny * Cmat(1,1)
  Cmat(2,1) = Cmat(1,2)
  Cmat(2,2) = Cmat(1,1)
  Cmat(3,3) = 0.5_dp * Emod*Thick / (1.0_dp + Rny)
  Cmat(4:6,4:6) = Cmat(1:3,1:3) * (Thick*Thick/12.0_dp)

  !! Set up element transformation matrix

  Trel = trans3P(X0(:,1),X0(:,2),X0(1,3),lpu,ierr)
  if (ierr < 0) then
     write(lpu,600) iEL
     ierr = -ierr
600  format('  ** The 3-noded shell element',I6, &
          & ' has zero area and is ignored.' &
          / '     Check your FE model for consistency.')
     return
  end if

  TelT = Trel(1:3,1:3)
  Tel  = transpose(TelT)

  !! Calculate the element stiffness matrix in local axes

  X21 = X0(1,2) - X0(1,1)
  Y21 = X0(2,2) - X0(2,1)
  Z21 = X0(3,2) - X0(3,1)
  X31 = X0(1,3) - X0(1,1)
  Y31 = X0(2,3) - X0(2,1)
  Z31 = X0(3,3) - X0(3,1)
  SL21 = sqrt(X21*X21 + Y21*Y21 + Z21*Z21)
  SL31 = sqrt(X31*X31 + Y31*Y31 + Z31*Z31)
  COSG = (X31*X21 + Y31*Y21 + Z31*Z21)/(SL31*SL21)
  if (COSG > -1.0_dp .and. COSG < 1.0_dp) then
     SING = sqrt(1.0_dp - COSG*COSG)
  else
     SING = 0.0_dp
  end if

  XL(1) = 0.0_dp
  XL(2) = SL21
  XL(3) = SL31*COSG
  YL(1) = 0.0_dp
  YL(2) = 0.0_dp
  YL(3) = SL31*SING

  call Andes3shell_stiffmat (Xl,Yl,Cmat,alpha,alphaH,lType,Kmat,lpu,ierr)
  if (ierr < 0) then
     write(lpu,690) iEL
690  format('*** Failed to compute stiffness matrix for element',I6)
     return
  end if

  !! Transform to global axes

  do i = 1, neldof, 3
     do j = 1, neldof, 3
        EK(i:i+2,j:j+2) = matmul(TelT,matmul(Kmat(i:i+2,j:j+2),Tel))
     end do
  end do

end subroutine IFEM_ANDES3


!> @brief Calculates the stiffness matrix for a 4-noded shell element.
subroutine IFEM_ANDES4 (iEL, X0, Thick, Emod, Rny, Rho, EK, IERR)

  use KindModule       , only : dp
  use Andes4ShellModule, only : Andes4shell_stiffmat
  use ManipMatrixModule, only : trans3p
  use PmatModule       , only : pMatStiff

  implicit none

  integer , parameter   :: nelnod = 4, neldof = 6*nelnod
  integer , intent(in)  :: iEL
  real(dp), intent(in)  :: X0(3,nelnod), Thick, Emod, Rny, Rho
  real(dp), intent(out) :: ek(neldof,neldof)
  integer , intent(out) :: ierr

  !! Local variables
  integer , parameter :: ltype = 2, lpu = 6
  real(dp), parameter :: alpha = 1.5_dp
  real(dp), parameter :: beta  = 0.9_dp

  integer  :: i, j
  real(dp) :: Tmp(neldof*neldof), Xtmp(3)
  real(dp) :: Kmat(neldof,neldof), Pmat(neldof,neldof), Cmat(6,6)
  real(dp) :: Xl(nelnod), Yl(nelnod), Zl(nelnod), Trel(3,4), Tel(3,3), TelT(3,3)


  !! --- Logic section ---

  ek = 0.0_dp

  !! Set up the constitutive matrix

  Cmat = 0.0_dp
  Cmat(1,1) = Emod*Thick / (1.0_dp - Rny*Rny)
  Cmat(1,2) = Rny * Cmat(1,1)
  Cmat(2,1) = Cmat(1,2)
  Cmat(2,2) = Cmat(1,1)
  Cmat(3,3) = 0.5_dp * Emod*Thick / (1.0_dp + Rny)
  Cmat(4:6,4:6) = Cmat(1:3,1:3) * (Thick*Thick/12.0_dp)

  !! Set up the element transformation matrix

  Xtmp = 0.5_dp*(X0(:,3)+X0(:,4))
  Trel = trans3P(X0(:,1),X0(:,2),Xtmp,lpu,ierr)
  if (ierr < 0) then
     write(lpu,600) iEL
     ierr = -ierr
600  format('  ** The 4-noded shell element',I6, &
          & ' has zero area and is ignored.' &
          / '     Check your FE model for consistency.')
     return
  end if

  TelT = Trel(1:3,1:3)
  Tel  = transpose(TelT)

  !! Calculate the element stiffness matrix
  !! Zl has to be close to zero, i.e., best fit in X-Y plane

  Xl(1) = 0.0_dp
  Yl(1) = 0.0_dp
  Zl(1) = 0.0_dp
  do i = 2, nelnod
     Xtmp  = X0(:,i) - X0(:,1)
     Xl(i) = dot_product(Tel(1,:),Xtmp)
     Yl(i) = dot_product(Tel(2,:),Xtmp)
     Zl(i) = dot_product(Tel(3,:),Xtmp)
  end do
  call Andes4shell_stiffmat (Xl,Yl,Zl,Cmat,alpha,beta,ltype,Kmat,LPU,IERR)
  if (ierr < 0) goto 999

  !! Transform to global axes
  do i = 1, neldof, 3
     do j = 1, neldof, 3
        Kmat(i:i+2,j:j+2) = matmul(TelT,matmul(Kmat(i:i+2,j:j+2),Tel))
     end do
  end do

  !! Compute projection matrix that restores stress free rigid body motions
  call pMatStiff (X0,PMAT,LPU,IERR)
  if (ierr < 0) goto 999

  !! Project the stiffness matrix: K = P'*K*P
  call DGEMM ('T','N', NELDOF, NELDOF, NELDOF, 1.0_dp, &
       &      PMAT(1,1), NELDOF, KMAT(1,1), NELDOF, 0.0_dp, TMP(1), NELDOF)
  call DGEMM ('N','N', NELDOF, NELDOF, NELDOF, 1.0_dp, &
       &      TMP(1), NELDOF, PMAT(1,1), NELDOF, 0.0_dp, EK(1,1), NELDOF)

  return

999 write(lpu,*) '*** Failed to compute stiffness matrix for element',iEL

end subroutine IFEM_ANDES4
