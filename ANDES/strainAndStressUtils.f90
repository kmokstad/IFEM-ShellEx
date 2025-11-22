!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file strainAndStressUtils.f90
!> @brief Utilities for strain- and stress calculation.

!!==============================================================================
!> @brief Module with utility subroutines for strain- and stress calculation.

module StrainAndStressUtilitiesModule

  implicit none

contains

  !!============================================================================
  !> @brief Computes strain-displacement matrix for constant strain triangle.
  !>
  !> @param[in] nndof Number of DOFs per nodal point
  !> @param[in] xEl X-coordinates of the element nodes
  !> @param[in] yEl Y-coordinates of the element nodes
  !> @param[in] zEl Z-coordinates of the element nodes
  !> @param[in] T_el Global-to-local transformation matrix for the element
  !> @param[in] zPos Local position of top shell surface w.r.t. the mid-surface
  !> @param[out] B_el Strain-displacement matrix for the element
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 13 Apr 2000

  subroutine StrainDispCST (nndof, xEl, yEl, zEl, T_el, zPos, B_el)

    use KindModule, only : dp, epsDiv0_p

    integer , intent(in)  :: nndof
    real(dp), intent(in)  :: xEl(:), yEl(:), zEl(:), T_el(3,3), zPos
    real(dp), intent(out) :: B_el(3,nndof,3)

    !! Local variables
    real(dp) :: xLij(3,3), yLij(3,3), vec(3), areaDobbel
    integer  :: in, jn

    !! --- Logic section ---

    B_el = 0.0_dp

    do in = 1, 3
       do jn = 1, 3
          if (in == jn) cycle
          vec(1) = xEl(in) - xEl(jn)
          vec(2) = yEl(in) - yEl(jn)
          vec(3) = zEl(in) - zEl(jn)
          xLij(in,jn) = dot_product(T_el(1,:),vec)
          yLij(in,jn) = dot_product(T_el(2,:),vec)
       end do
    end do
    areaDobbel = xLij(2,1)*yLij(3,1) - xLij(3,1)*yLij(2,1)

    !! X-displacements
    B_el(1,1,1) =  yLij(2,3)/areaDobbel
    B_el(1,1,2) =  yLij(3,1)/areaDobbel
    B_el(1,1,3) =  yLij(1,2)/areaDobbel

    B_el(3,1,1) = -xLij(2,3)/areaDobbel
    B_el(3,1,2) = -xLij(3,1)/areaDobbel
    B_el(3,1,3) = -xLij(1,2)/areaDobbel

    !! Y-displacements
    B_el(2,2,1) = -xLij(2,3)/areaDobbel
    B_el(2,2,2) = -xLij(3,1)/areaDobbel
    B_el(2,2,3) = -xLij(1,2)/areaDobbel

    B_el(3,2,1) =  yLij(2,3)/areaDobbel
    B_el(3,2,2) =  yLij(3,1)/areaDobbel
    B_el(3,2,3) =  yLij(1,2)/areaDobbel

    !! Rotational degrees of freedom
    if (nndof >= 5 .and. abs(zpos) > epsDiv0_p) then
       do in = 1, 3
          B_el(:,4,in) = -zpos*B_el(:,2,in) ! X-rotation
          B_el(:,5,in) =  zpos*B_el(:,1,in) ! Y-rotation
       end do
    end if

  end subroutine StrainDispCST


  !!============================================================================
  !> @brief Computes strain-displacement matrix for the 4-noded quadrilateral.
  !>
  !> @param[in] nndof Number of DOFs per nodal point
  !> @param[in] xEl X-coordinates of the element nodes
  !> @param[in] yEl Y-coordinates of the element nodes
  !> @param[in] zEl Z-coordinates of the element nodes
  !> @param[in] T_el Global-to-local transformation matrix for the element
  !> @param[in] xi First coordinate in domain [-1,1] of evaluation point
  !> @param[in] eta Second coordinate in domain [-1,1] of evaluation point
  !> @param[in] zPos Local position of top shell surface w.r.t. the mid-surface
  !> @param[out] B_el Strain-displacement matrix for the element
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 13 Apr 2000

  subroutine StrainDispQuad4 (nndof, xEl, yEl, zEl, T_el, xi, eta, zPos, B_el)

    use KindModule       , only : dp, epsDiv0_p
    use manipMatrixModule, only : invert22

    integer , intent(in)  :: nndof
    real(dp), intent(in)  :: xEl(:), yEl(:), zEl(:), T_el(3,3), xi, eta, zPos
    real(dp), intent(out) :: B_el(3,nndof,4)

    !! Local variables
    real(dp) :: sfuncD(2,4), jacobi(2,2), jacobi_inv(2,2), XY(4,2), vec(3)
    integer  :: in

    !! --- Logic section ---

    B_el = 0.0_dp

    !! Local XY-coordinates relative to node 1
    XY(1,:) = 0.0_dp
    do in = 2, 4
       vec(1) = xEl(in) - xEl(1)
       vec(2) = yEl(in) - yEl(1)
       vec(3) = zEl(in) - zEl(1)
       XY(in,:) = matmul(T_el(1:2,:),vec)
    end do

    !! Partial derivatives with respect to xi
    sfuncD(1,1) = -(1.0_dp - eta)*0.25_dp
    sfuncD(1,2) =  (1.0_dp - eta)*0.25_dp
    sfuncD(1,3) =  (1.0_dp + eta)*0.25_dp
    sfuncD(1,4) = -(1.0_dp + eta)*0.25_dp

    !! Partial derivatives with respect to eta
    sfuncD(2,1) = -(1.0_dp - xi)*0.25_dp
    sfuncD(2,2) = -(1.0_dp + xi)*0.25_dp
    sfuncD(2,3) =  (1.0_dp + xi)*0.25_dp
    sfuncD(2,4) =  (1.0_dp - xi)*0.25_dp

    !! Jacobi matrix  | dx/dxi   dy/dxi  |
    !!                | dx/deta  dy/deta |
    jacobi = matmul(sfuncD,XY)

    !! Invert the Jacobian, i.e.  | dxi/dx  deta/dx |
    !!                            | dxi/dy  deta/dy |
    jacobi_inv = invert22(jacobi)

    !! X-displacements
    B_el(1,1,:) = jacobi_inv(1,1)*sfuncD(1,:) + jacobi_inv(1,2)*sfuncD(2,:)
    B_el(3,1,:) = jacobi_inv(2,1)*sfuncD(1,:) + jacobi_inv(2,2)*sfuncD(2,:)

    !! Y-displacements
    B_el(2,2,:) = B_el(3,1,:)
    B_el(3,2,:) = B_el(1,1,:)

    !! Rotational degrees of freedom
    if (nndof >= 5 .and. abs(zpos) > epsDiv0_p) then
       do in = 1, 4
          B_el(:,4,in) = -zpos*B_el(:,2,in) ! X-rotation
          B_el(:,5,in) =  zpos*B_el(:,1,in) ! Y-rotation
       end do
    end if

  end subroutine StrainDispQuad4


  !!============================================================================
  !> @brief Calculates the X-axis direction of the globalized coordinate system.
  !>
  !> @param[in] eZ Local Z-axis (normal vector) of the shell surface
  !> @return Globalized X-axis direction vector for the shell surface
  !>
  !> @details This function computes the vector @a V1 defined by the projection
  !> of the global X-axis onto the plane defined by the normal vector @a eZ.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 30 Nov 2000

  function getGlobalizedX (eZ) result(V1)

    use KindModule       , only : dp, epsDiv0_p
    use manipMatrixModule, only : cross_product

    real(dp), intent(in) :: eZ(3)
    real(dp)             :: V1(3), V2(3), VLength2
    real(dp), parameter  :: somewhatSmall_p = 0.01_dp

    !! --- Logic section ---

    if (abs(eZ(2)) > somewhatSmall_p .or. abs(eZ(3)) > somewhatSmall_p) then
       !! Define V1 by projecting global X-axis onto the plane
       V1(1) =  eZ(2)*eZ(2) + eZ(3)*eZ(3)
       V1(2) = -eZ(1)*eZ(2)
       V1(3) = -eZ(1)*eZ(3)
    else
       !! Define V2 by projecting global Y-axis onto the plane, then V1={V2}x{eZ}
       V2(1) = -eZ(2)*eZ(1)
       V2(2) =  eZ(1)*eZ(1) + eZ(3)*eZ(3)
       V2(3) = -eZ(2)*eZ(3)
       V1    = cross_product(V2,eZ)
    end if

    VLength2 = dot_product(V1,V1)
    if (VLength2 > epsDiv0_p*epsDiv0_p) then
       V1 = V1 / sqrt(VLength2)
    else
       V1 = 0.0_dp
    end if

  end function getGlobalizedX


  !!============================================================================
  !> @brief Computes the local element axes for a thin shell element.
  !>
  !> @param[in] nenod Number of element nodes (3 or 4)
  !> @param[in] X Global X-coordinates for the element
  !> @param[in] Y Global Y-coordinates for the element
  !> @param[in] Z Global Z-coordinates for the element
  !> @param[out] V1 Direction of local X-axis for the element
  !> @param[out] V2 Direction of local Y-axis for the element
  !> @param[out] V3 Direction of local Z-axis for the element
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !> @param[in] doGlobalize If .true., the globalized shell axes are computed
  !>
  !> @details The local-to-global transformation matrix is then T = [V1,V2,V3].
  !> Optionally, the globalized shell element axes can be computed instead of
  !> the local element axes.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 30 Nov 2000

  subroutine getShellElementAxes (nenod, X, Y, Z, V1, V2, V3, &
       &                          lpu, ierr, doGlobalize)

    use KindModule       , only : dp, epsDiv0_p
    use manipMatrixModule, only : cross_product

    integer , intent(in)  :: nenod, lpu
    real(dp), intent(in)  :: X(nenod), Y(nenod) ,Z(nenod)
    real(dp), intent(out) :: V1(3), V2(3), V3(3)
    integer , intent(out) :: ierr
    logical , optional    :: doGlobalize

    !! Local variables
    real(dp) :: VN
    logical  :: useGlobalized

    !! --- Logic section ---

    ierr = 0
    if (present(doGlobalize)) then
       useGlobalized = doGlobalize
    else
       useGlobalized = .false.
    end if

    !! Define the shell normal vector (local Z-axis)
    if (nenod == 3) then
       !! Use X12 and X13
       V1(1) = X(2) - X(1)
       V1(2) = Y(2) - Y(1)
       V1(3) = Z(2) - Z(1)
       V2(1) = X(3) - X(1)
       V2(2) = Y(3) - Y(1)
       V2(3) = Z(3) - Z(1)
    else if (nenod == 4) then
       !! Use the diagonals X13 and X24
       V1(1) = X(3) - X(1)
       V1(2) = Y(3) - Y(1)
       V1(3) = Z(3) - Z(1)
       V2(1) = X(4) - X(2)
       V2(2) = Y(4) - Y(2)
       V2(3) = Z(4) - Z(2)
    else
       write(lpu,*) '*** getShellElementAxes: Invalid element type ',nenod
       ierr = -1
       return
    end if

    V3 = cross_product(V1,V2)
    VN = dot_product(V3,V3)
    if (VN > epsDiv0_p*epsDiv0_p) then
       V3 = V3 / sqrt(VN)
    else
       ierr = 2
       goto 900
    end if

    !! Define the local- or globalized X-axis
    if (useGlobalized) then
       V1 = getGlobalizedX(V3)
    else if (nenod == 4) then
       !! Use the vector between node 1 and 2 (equivalent to Femlib::DIRC_QUAD)
       V1(1) = X(2) - X(1)
       V1(2) = Y(2) - Y(1)
       V1(3) = Z(2) - Z(1)
       !! Project V1 onto the plane defined by the normal vector V3
       V2 = cross_product(V3,V1)
       V1 = cross_product(V2,V3)
    end if
    VN = dot_product(V1,V1)
    if (VN > epsDiv0_p*epsDiv0_p) then
       V1 = V1 / sqrt(VN)
    else
       ierr = 3
       goto 900
    end if

    !! And finally the Y-axis
    V2 = cross_product(V3,V1)

    return

900 continue
    write(lpu,*) '*** getShellElementAxes: Degenerated shell element',ierr
    write(lpu,600) V1,V2,V3,7-2*ierr,7-2*ierr,VN,epsDiv0_p*epsDiv0_p
600 format(26X,'V1 =',1P3E13.5 / 26X,'V2 =',1P3E13.5 / 26X,'V3 =',1P3E13.5, &
         / 26X,'V',I1,'*V',I1,' = ',1PE12.5,' tol = ',E12.5)

  end subroutine getShellElementAxes


  !!============================================================================
  !> @brief Computes the stress transformation matrix for a thin shell element.
  !>
  !> @param[in] eX X-axis of the element coordinate system
  !> @param[in] eZ Z-axis of the element coordinate system
  !> @param[out] T In-plane transformation matrix
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @details This subroutine computes the 2D transformation matrix from the
  !> element coordinate system to the continuous stress output coordinate
  !> system for a shell element.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 30 Nov 2000

  subroutine getShellStressTrans (eX, eZ, T, lpu, ierr)

    use KindModule       , only : dp, epsDiv0_p
    use manipMatrixModule, only : cross_product

    real(dp), intent(in)  :: eX(3), eZ(3)
    real(dp), intent(out) :: T(2,2)
    integer , intent(in)  :: lpu
    integer , intent(out) :: ierr

    !! Local variables
    real(dp) :: V1(3), V2(3), CA, SA

    !! --- Logic section ---

    !! Find the local X-axis (V1) of the stress output coordinate system
    V1 = getGlobalizedX(eZ)
    if (dot_product(V1,V1) <= epsDiv0_p*epsDiv0_p) then
       write(lpu,*) '*** getShellStressTrans: Invalid shell normal vector',eZ
       ierr = -1
       return
    end if

    !! Find the rotation from aX to V1 in terms of cos(angle) and sin(angle)
    V2 = cross_product(eX,V1)
    CA = dot_product(V1,eX)
    SA = sign(sqrt(dot_product(V2,V2)),dot_product(V2,eZ))

    !! Set up the in-plane rotation matrix
    T(1,1) =  CA
    T(1,2) =  SA
    T(2,1) = -SA
    T(2,2) =  CA

    ierr = 0

  end subroutine getShellStressTrans


  !!============================================================================
  !> @brief Transforms 2D strains to a globalized shell coordinate system.
  !>
  !> @param[in] eX X-axis of the element coordinate system
  !> @param[in] eZ Z-axis of the element coordinate system
  !> @param eps 2D strain tensor components
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @callergraph @callgraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 11 Nov 2025

  subroutine transform2Dstrain (eX, eZ, eps, lpu, ierr)

    use KindModule                  , only : dp
    use FFaTensorTransformsInterface, only : tratensor

    real(dp), intent(in)    :: eX(3), eZ(3)
    real(dp), intent(inout) :: eps(:)
    integer , intent(in)    :: lpu
    integer , intent(out)   :: ierr

    !! Local variables
    integer  :: i
    real(dp) :: T_str(2,2)

    !! --- Logic section ---

    call getShellStressTrans (eX, eZ, T_str, lpu, ierr)
    if (ierr < 0) return

    do i = 1, size(eps), 3
       eps(i+2) = 0.5_dp * eps(i+2) ! to tensorial shear strain, epsilon_xy
       call tratensor (2,eps(i:i+2),T_str)
       eps(i+2) = 2.0_dp * eps(i+2) ! back to enginering shear strain, gamma_xy
    end do

  end subroutine transform2Dstrain

end module StrainAndStressUtilitiesModule
