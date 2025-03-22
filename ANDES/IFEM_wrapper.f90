!! $Id$
!!==============================================================================
!> @file IFEM_wrapper.f90
!> @brief Wrapper on the Andes shell elements from FEDEM (https://openfedem.org)
!>
!> @author Knut Morten Okstad / SINTEF
!>
!> $date Feb 20 2024

!> @brief Calculates the stiffness- and mass matrix for a 3-noded shell element.
subroutine IFEM_ANDES3 (iEL, X0, Thick, Emod, Rny, Rho, Press, EK, EM, ES, IERR)

  use KindModule        , only : dp
  use Andes3ShellModule , only : Andes3shell_stiffmat
  use ManipMatrixModule , only : trans3p

  implicit none

  integer , parameter   :: nelnod = 3, neldof = 6*nelnod
  integer , intent(in)  :: iEL
  real(dp), intent(in)  :: X0(3,nelnod), Thick, Emod, Rny, Rho, Press(3,nelnod)
  real(dp), intent(out) :: ek(neldof,neldof), em(neldof,neldof), es(neldof)
  integer , intent(out) :: ierr

  !! Local variables
  integer , parameter :: ltype  = 2, lpu = 6
  real(dp), parameter :: alpha  = 1.5_dp
  real(dp), parameter :: alphaH = 0.5_dp
  real(dp), parameter :: thetaMaterial = 0.0_dp

  integer  :: i, j, k
  real(dp) :: X21, Y21, Z21, X31, Y31, Z31, SL21, SL31, COSG, SING
  real(dp) :: Kmat(neldof,neldof), Cmat(6,6), Tmp(neldof)
  real(dp) :: Xl(nelnod), Yl(nelnod), Zl(nelnod), Trel(3,4), Tel(3,3), TelT(3,3)


  !! --- Logic section ---

  if (Thick > 0.0_dp) then

     call DCOPY (NELDOF*NELDOF,0.0_dp,0,EK(1,1),1)

     !! Set up the constitutive matrix

     Cmat = 0.0_dp
     Cmat(1,1) = Emod*Thick / (1.0_dp - Rny*Rny)
     Cmat(1,2) = Rny * Cmat(1,1)
     Cmat(2,1) = Cmat(1,2)
     Cmat(2,2) = Cmat(1,1)
     Cmat(3,3) = 0.5_dp * Emod*Thick / (1.0_dp + Rny)
     Cmat(4:6,4:6) = Cmat(1:3,1:3) * (Thick*Thick/12.0_dp)

  end if

  !! Set up element transformation matrix

  Trel = trans3P(X0(:,1),X0(:,2),X0(1,3),lpu,ierr)
  if (ierr < 0) then
     write(lpu,600) iEL
     ierr = -ierr
     do i = 1, 3
        write(lpu,601) CHAR(ICHAR('W')+i), X0(i,:)
     end do
600  format('  ** The 3-noded shell element',I8, &
          & ' has zero area and is ignored.' &
          / '     Check your FE model for consistency.')
601  format('     ',A1,' =',1P3E13.5)
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

  if (Thick > 0.0_dp) then

     call Andes3shell_stiffmat (Xl,Yl,Cmat,alpha,alphaH,lType,Kmat,lpu,ierr)
     if (ierr < 0) then
        ierr = 9
        write(lpu,*) '*** Failed to compute stiffness matrix for element',iEL
        do i = 1, 3
           write(lpu,601) CHAR(ICHAR('W')+i), X0(i,:)
        end do
        return
     end if

     !! Transform to global axes

     do i = 1, neldof, 3
        do j = 1, neldof, 3
           EK(i:i+2,j:j+2) = matmul(TelT,matmul(Kmat(i:i+2,j:j+2),Tel))
        end do
     end do

  end if

  if (any(Press /= 0.0_dp)) call T3_load (Press,ES)

  if (Rho <= 0.0_dp .or. Thick <= 0.0_dp) return

  !! Compute lumped mass matrix
  Xl = X0(1,:)
  Yl = X0(2,:)
  Zl = X0(3,:)
  call TMRF35 (Tmp(1),Xl(1),Yl(1),Zl(1),Thick,Rho)
  call TEBA35 (Tmp(10),Xl(1),Yl(1),Zl(1),Thick,Rho)
  Tmp(6)   = Tmp(3)
  Tmp(3:5) = Tmp(10:12)
  do i = 1, 6
     do j = 6, 12, 6
        Tmp(i+j) = Tmp(i)
     end do
  end do
  call DCOPY (NELDOF*NELDOF,0.0_dp,0,EM(1,1),1)
  call DCOPY (NELDOF,Tmp(1),1,EM(1,1),NELDOF+1)

contains

  !> @brief Calculates consistent load vector for a 3-noded shell element.
  subroutine T3_load (Pglob,ES)
    real(dp), intent(in)  :: Pglob(:,:)
    real(dp), intent(out) :: ES(:)
    real(dp)              :: Ploc(3,nelnod), Ao12, b(3)
    !! Transform from global to local load intensities
    Ploc = matmul(TelT,Pglob)
    !! Element area (divided by 12)
    Ao12 = ((XL(1)*YL(2) + XL(2)*YL(3) + XL(3)*YL(1)) &
         -  (XL(1)*YL(3) + XL(2)*YL(1) + XL(3)*YL(2))) / 24.0_dp
    !! Determine non-zero terms of the local load vector,
    !! transform to global and add into the element load vector
    do i = 1, 3
       j = 1 + mod(i,3)
       k = 1 + mod(j,3)
       b = Ao12*(2.0_dp*Ploc(:,i) + Ploc(:,j) + Ploc(:,k))
       ES(6*i-5:6*i-3) = matmul(Tel,b)
    end do
  end subroutine T3_load

end subroutine IFEM_ANDES3


!> @brief Calculates the stiffness- and mass matrix for a 4-noded shell element.
subroutine IFEM_ANDES4 (iEL, X0, Thick, Emod, Rny, Rho, EK, EM, IERR)

  use KindModule       , only : dp
  use Andes4ShellModule, only : Andes4shell_stiffmat
  use ManipMatrixModule, only : trans3p
  use PmatModule       , only : pMatStiff

  implicit none

  integer , parameter   :: nelnod = 4, neldof = 6*nelnod
  integer , intent(in)  :: iEL
  real(dp), intent(in)  :: X0(3,nelnod), Thick, Emod, Rny, Rho
  real(dp), intent(out) :: ek(neldof,neldof), em(neldof,neldof)
  integer , intent(out) :: ierr

  !! Local variables
  integer , parameter :: ltype = 2, lpu = 6
  real(dp), parameter :: alpha = 1.5_dp
  real(dp), parameter :: beta  = 0.9_dp

  integer  :: i, j
  real(dp) :: Tmp(neldof*neldof), Xtmp(3), Xnod(3,nelnod)
  real(dp) :: Kmat(neldof,neldof), Pmat(neldof,neldof), Cmat(6,6)
  real(dp) :: Xl(nelnod), Yl(nelnod), Zl(nelnod), Trel(3,4), Tel(3,3), TelT(3,3)


  !! --- Logic section ---

  ierr = 0
  if (Thick <= 0.0_dp) return

  call DCOPY (NELDOF*NELDOF,0.0_dp,0,EK(1,1),1)

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
     do i = 1, 3
        write(lpu,601) CHAR(ICHAR('W')+i), X0(i,:)
     end do
600  format('  ** The 4-noded shell element',I8, &
          & ' has zero area and is ignored.' &
          / '     Check your FE model for consistency.')
601  format('     ',A1,' =',1P4E13.5)
     return
  end if

  TelT = Trel(1:3,1:3)
  Tel  = transpose(TelT)

  !! Need to swap element nodes 3 and 4
  Xnod(:,1:2) = X0(:,1:2)
  Xnod(:,3)   = X0(:,4)
  Xnod(:,4)   = X0(:,3)

  !! Calculate the element stiffness matrix
  !! Zl has to be close to zero, i.e., best fit in X-Y plane

  Xl(1) = 0.0_dp
  Yl(1) = 0.0_dp
  Zl(1) = 0.0_dp
  do i = 2, nelnod
     Xtmp  = Xnod(:,i) - Xnod(:,1)
     Xl(i) = dot_product(Tel(1,:),Xtmp)
     Yl(i) = dot_product(Tel(2,:),Xtmp)
     Zl(i) = dot_product(Tel(3,:),Xtmp)
  end do
  call Andes4shell_stiffmat (Xl,Yl,Zl,Cmat,alpha,beta,ltype,Kmat,LPU,IERR)
  if (ierr < 0) then
     ierr = 9
     goto 999
  end if

  !! Transform to global axes
  do i = 1, neldof, 3
     do j = 1, neldof, 3
        Kmat(i:i+2,j:j+2) = matmul(TelT,matmul(Kmat(i:i+2,j:j+2),Tel))
     end do
  end do

  !! Compute projection matrix that restores stress free rigid body motions
  call pMatStiff (Xnod,PMAT,LPU,IERR)
  if (ierr < 0) goto 999

  !! Project the stiffness matrix: K = P^t*K*P
  call DGEMM ('T','N', NELDOF, NELDOF, NELDOF, 1.0_dp, &
       &      PMAT(1,1), NELDOF, KMAT(1,1), NELDOF, 0.0_dp, TMP(1), NELDOF)
  call DGEMM ('N','N', NELDOF, NELDOF, NELDOF, 1.0_dp, &
       &      TMP(1), NELDOF, PMAT(1,1), NELDOF, 0.0_dp, EK(1,1), NELDOF)

  !! Swap contributions from element nodes 3 and 4
  do i = 13, 18
     call swapRowCol (EK,i,i+6)
  end do

  if (Rho <= 0.0_dp) return

  !! Compute lumped mass matrix
  Xl = Xnod(1,:)
  Yl = Xnod(2,:)
  Zl = Xnod(3,:)
  call QMRF35 (Tmp(1),Xl(1),Yl(1),Zl(1),Thick,Rho)
  call QBEK35 (Tmp(1),Xl(1),Yl(1),Zl(1),Thick,Rho)
  Tmp(25:30) = Tmp(13:18)
  Tmp(13:24) = Tmp(19:30)
  call DCOPY (NELDOF*NELDOF,0.0_dp,0,EM(1,1),1)
  call DCOPY (NELDOF,Tmp(1),1,EM(1,1),NELDOF+1)

  return

999 write(lpu,*) '*** Failed to compute stiffness matrix for element',iEL
  do i = 1, 3
     write(lpu,601) CHAR(ICHAR('W')+i), X0(i,:)
  end do

contains

  !> @brief Swaps two rows and columns in the stiffness matrix
  subroutine swapRowCol (A,i,j)
    real(dp), intent(inout) :: A(:,:)
    integer , intent(in)    :: i, j
    Tmp(1:neldof) = A(:,i)
    A(:,i) = A(:,j)
    A(:,j) = Tmp(1:neldof)
    Tmp(1:neldof) = A(i,:)
    A(i,:) = A(j,:)
    A(j,:) = Tmp(1:neldof)
  end subroutine swapRowCol

end subroutine IFEM_ANDES4


!> @brief Calculates FE stresses at the center of a 4-noded shell element.
subroutine IFEM_STRS24 (iEL, X0, Thick, Emod, Rny, EV, SR, Sigma, IERR)

  use KindModule                    , only : dp
  use PmatModule                    , only : pMatStiff
  use StrainAndStressUtilitiesModule, only : getShellElementAxes
  use StrainAndStressUtilitiesModule, only : getShellStressTrans
  use StrainAndStressUtilitiesModule, only : strainDispQuad4

  implicit none

  integer , parameter   :: nelnod = 4, neldof = 6*nelnod, LPU = 6
  integer , intent(in)  :: iEL
  real(dp), intent(in)  :: X0(3,nelnod), Thick, Emod, Rny
  real(dp), intent(in)  :: EV(neldof) !< Element nodal displacements
  real(dp), intent(out) :: SR(6)      !< Stress resultants at element centre
  real(dp), intent(out) :: Sigma(6)   !< Stresses at element center
  integer , intent(out) :: IERR

  !! Local variables
  integer  :: i
  real(dp) :: T_el(3,3), T_str(2,2), Cmat(3,3), epsU(3), epsL(3), Xnod(3,nelnod)
  real(dp) :: PMAT(neldof,neldof), vld(neldof), B_L(3,neldof), B_U(3,neldof)

  !! --- Logic section ---

  !! Need to swap element nodes 3 and 4
  Xnod(:,1:2) = X0(:,1:2)
  Xnod(:,3)   = X0(:,4)
  Xnod(:,4)   = X0(:,3)

  !! Compute the global-to-local transformation matrix
  call getShellElementAxes (nelnod, Xnod(1,:), Xnod(2,:), Xnod(3,:), &
       &                    T_el(1,:), T_el(2,:), T_el(3,:), LPU, IERR)
  if (ierr /= 0) goto 999

  !! Compute the 2D stress transformation matrix
  call getShellStressTrans (T_el(1,:), T_el(3,:), T_str, LPU, IERR)
  if (ierr /= 0) goto 999

  !! Evaluate the strain-displacement matrix in the element center
  call StrainDispQuad4 (6, Xnod(1,:), Xnod(2,:), Xnod(3,:), &
       &                T_el, 0.0_dp, 0.0_dp, 0.5_dp*Thick, B_U)
  call StrainDispQuad4 (6, Xnod(1,:), Xnod(2,:), Xnod(3,:), &
       &                T_el, 0.0_dp, 0.0_dp,-0.5_dp*Thick, B_L)

  !! Compute projection matrix
  call pMatStiff (Xnod,PMAT,LPU,IERR)
  if (ierr < 0) goto 999

  !! Compute deformational displacement vector through projection
  vld = matmul(PMAT,EV)

  !! Form local coordinate displacement vector
  do i = 1, neldof, 3
     vld(i:i+2) = matmul(T_el,vld(i:i+2))
  end do

  !! Strains at element center, upper and lower surface
  epsU = matmul(B_U,vld)
  epsL = matmul(B_L,vld)

  !! Transform to stress output coordinate system
  call traStrain (epsU,T_str)
  call traStrain (epsL,T_str)

  !! Set up the constitutive matrix

  Cmat = 0.0_dp
  Cmat(1,1) = Emod*Thick / (1.0_dp - Rny*Rny)
  Cmat(1,2) = Rny * Cmat(1,1)
  Cmat(2,1) = Cmat(1,2)
  Cmat(2,2) = Cmat(1,1)
  Cmat(3,3) = 0.5_dp * Emod*Thick / (1.0_dp + Rny)

  !! Stresses at element centre
  Sigma(1:3) = matmul(Cmat,epsU)
  Sigma(4:6) = matmul(Cmat,epsL)

  !! Stress resultants at element centre
  SR(1:3) = (Sigma(1:3) + Sigma(4:6)) * Thick/2.0_dp
  SR(4:6) = (Sigma(1:3) - Sigma(4:6)) * Thick*Thick/12.0_dp

  return

999 write(lpu,*) '*** Failed to compute stresses for element',iEL

contains

  !> @brief Transforms a 2D strain tensor.
  subroutine traStrain (eps,T_str)
#ifdef HAS_FFLLIB
    use FFaTensorTransformsInterface, only : tratensor
#endif
    real(dp), intent(inout) :: eps(3)
    real(dp), intent(in)    :: T_str(2,2)
#ifdef HAS_FFLLIB
    eps(3) = 0.5_dp * eps(3) ! to tensorial shear strain, epsilon_xy
    call tratensor (2,eps,T_str)
    eps(3) = 2.0_dp * eps(3) ! back to enginering shear strain, gamma_xy
#endif
  end subroutine traStrain

end subroutine IFEM_STRS24
