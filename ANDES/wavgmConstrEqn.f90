!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file wavgmConstrEqn.f90
!>
!> @brief Weighted Average Motion constraint handling.

!!==============================================================================
!> @brief Computes constraint equation coefficients for a WAVGM element.
!>
!> @param[in] iel Element index
!> @param[in] lDof Local index of dependent DOF to compute coefficients for
!> @param[in] nM Number of independent nodes in current element
!> @param[in] nW SiNumber of independent nodes in current element
!> @param[in] indC Nodal component indices (common for all nodes)
!> @param[in] tenc Table of nodal coordinates for current element
!> @param[in] weight Independent DOF weights for current element
!> @param[in] epsX Relative geometric tolerance for WAVGM elements
!> @param dX Element nodal coordinates relative to the centre of gravity
!> @param work Work array
!> @param[out] omega Resulting constraint equation coefficients
!> @param[in] ipsw Print switch
!> @param[in] lpu File unit number for res-file output
!>
!> @details This subroutine pre-processes a Weighted Average Motion (WAVGM)
!> elements (also known as RBE3 in Nastran). The coefficients that couples
!> a dependent DOF at the reference node of a Weighted AVerage Motion element
!> to the independent nodal DOFs of the same element are computed.
!>
!> @callergraph
!>
!> @author Knut Morten Okstad
!>
!> @date 24 Sep 2002

subroutine wavgmConstrEqn (iel,lDof,nM,nW,indC,tenc,weight, &
     &                     epsX,dX,work,omega,ipsw,lpu)

  use kindModule       , only : dp, epsDiv0_p
  use manipMatrixModule, only : writeObject

  implicit none

  logical , parameter     :: reComputeCG = .true.
  integer , parameter     :: nndof = 6
  integer , intent(in)    :: iel, lDof, nM, nW, indC(nndof), ipsw, lpu
  real(dp), intent(in)    :: tenc(3,1+nM), weight(nW), epsX
  real(dp), intent(inout) :: dX(3,1+nM), work(nM)
  real(dp), intent(out)   :: omega(nndof*nM)

  !! Local variables
  integer           :: i, j, k, iFrst, iLast
  integer , save    :: lastIC = 0, lastEl = 0
  real(dp), save    :: tolX(3) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
  real(dp)          :: sumWM, tolWM
  character(len=32) :: label

  !! --- Logic section ---

  omega = 0.0_dp
  if (lDof < 1 .or. lDof > nndof) return

  if (iel /= lastEl) then
     lastEl = iel
     lastIC = 0
  end if

  !!                        k
  !! Cyclic permutation:   / \
  !!                      i - j
  i = 1 + mod(lDof-1,3)
  j = 1 + mod(i,3)
  k = 1 + mod(j,3)

  if (indC(lDof) > 0) then

     !! Direct coupling of dependent DOF to the corresponding independent DOFs

     iFrst = indC(lDof)
     iLast = iFrst + nM-1
     if (iLast > nW) goto 9
     sumWM = sum(weight(iFrst:iLast))
     if (sumWM > epsDiv0_p) then
        call DAXPY (nM,1.0_dp/sumWM,weight(iFrst),1,omega(lDof),nndof)
     else if (ipsw > 0) then
        write(lpu,600) iel,lDof,lDof,sumWM,epsDiv0_p
     end if

  else if (indC(lDof) < 0) then

     !! Direct coupling with uniform weights
     call DCOPY (nM,1.0_dp/real(nM,dp),0,omega(lDof),nndof)

  end if

  if (lDof > 3) then

     !! Coupling of rotational dependent DOF to translational DOFs j,k

     if (indC(j) /= 0) then
        iFrst = indC(j)
        iLast = iFrst + nM-1
        if (iFrst > 0 .and. iLast > nW) goto 9
        if (reComputeCG .and. iFrst /= lastIC) then
           lastIC = iFrst
           if (iFrst > 0) then
              call computePointCoords (weight(iFrst:iLast))
           else ! uniform weights
              call computePointCoords ()
           end if
        end if
        if (iFrst > 0) then
           sumWM = sumSqrW(dX(j,2:1+nM),dX(k,2:1+nM),weight(iFrst:iLast))
           work(1:nM) = dX(k,2:1+nM)*weight(iFrst:iLast)
        else ! uniform weights
           sumWM = sumSqr(dX(j,2:1+nM),dX(k,2:1+nM)) * real(nM,dp)
           call DCOPY (nM,dX(k,2),size(dX,1),work(1),1)
        end if
        tolWM = tolX(j)*tolX(j) + tolX(k)*tolX(k)
        if (sumWM > tolWM) then
           call DAXPY (nM,-1.0_dp/sumWM,work(1),1,omega(j),nndof)
        else if (ipsw > 0) then
           write(lpu,600) iel,lDof,j,sumWM,tolWM
           write(lpu,620) work(1:nM) / sumWM
        end if
     end if
     if (indC(k) /= 0) then
        iFrst = indC(k)
        iLast = iFrst + nM-1
        if (iFrst > 0 .and. iLast > nW) goto 9
        if (reComputeCG .and. iFrst /= lastIC) then
           lastIC = iFrst
           if (iFrst > 0) then
              call computePointCoords (weight(iFrst:iLast))
           else ! uniform weights
              call computePointCoords ()
           end if
        end if
        if (iFrst > 0) then
           sumWM = sumSqrW(dX(j,2:1+nM),dX(k,2:1+nM),weight(iFrst:iLast))
           work(1:nM) = dX(j,2:1+nM)*weight(iFrst:iLast)
        else ! uniform weights
           sumWM = sumSqr(dX(j,2:1+nM),dX(k,2:1+nM)) * real(nM,dp)
           call DCOPY (nM,dX(j,2),size(dX,1),work(1),1)
        end if
        tolWM = tolX(j)*tolX(j) + tolX(k)*tolX(k)
        if (sumWM > tolWM) then
           call DAXPY (nM,1.0_dp/sumWM,work(1),1,omega(k),nndof)
        else if (ipsw > 0) then
           write(lpu,600) iel,lDof,k,sumWM,tolWM
           write(lpu,620) work(1:nM) / sumWM
        end if
     end if

  else if (lDof < 4) then

     !! Coupling of translational dependent DOF to translational DOFs i,j
     !! and rotational DOF 3+k through the eccentricity e_j
     !! of the dependent node w.r.t. the independent node centroid

     if (indC(i) /= 0) then
        iFrst = indC(i)
        iLast = iFrst + nM-1
        if (iFrst > 0 .and. iLast > nW) goto 9
        if (reComputeCG .and. iFrst /= lastIC) then
           lastIC = iFrst
           if (iFrst > 0) then
              call computePointCoords (weight(iFrst:iLast))
           else ! uniform weights
              call computePointCoords ()
           end if
        end if
        if (abs(dX(j,1)) > tolX(j)) then
           if (iFrst > 0) then
              sumWM = sumSqrW(dX(i,2:1+nM),dX(j,2:1+nM),weight(iFrst:iLast))
              work(1:nM) = dX(j,2:1+nM)*weight(iFrst:iLast)
           else
              sumWM = sumSqr(dX(i,2:1+nM),dX(j,2:1+nM)) * real(nM,dp)
              call DCOPY (nM,dX(j,2),size(dX,1),work(1),1)
           end if
           tolWM = tolX(i)*tolX(i) + tolX(j)*tolX(j)
           if (sumWM > tolWM) then
              call DAXPY (nM,dX(j,1)/sumWM,work(1),1,omega(i),nndof)
           else if (ipsw > 0) then
              write(lpu,600) iel,lDof,i,sumWM,tolWM
              write(lpu,620) work(1:nM) * dX(j,1)/sumWM
           end if
        else if (ipsw > 0) then
           write(lpu,610) iel,lDof,i,dX(j,1),tolX(j)
        end if
     end if
     if (indC(j) /= 0) then
        iFrst = indC(j)
        iLast = iFrst + nM-1
        if (iFrst > 0 .and. iLast > nW) goto 9
        if (reComputeCG .and. iFrst /= lastIC) then
           lastIC = iFrst
           if (iFrst > 0) then
              call computePointCoords (weight(iFrst:iLast))
           else ! uniform weights
              call computePointCoords ()
           end if
        end if
        if (abs(dX(j,1)) > tolX(j)) then
           if (iFrst > 0) then
              sumWM = sumSqrW(dX(i,2:1+nM),dX(j,2:1+nM),weight(iFrst:iLast))
              work(1:nM) = dX(i,2:1+nM)*weight(iFrst:iLast)
           else
              sumWM = sumSqr(dX(i,2:1+nM),dX(j,2:1+nM)) * real(nM,dp)
              call DCOPY (nM,dX(i,2),size(dX,1),work(1),1)
           end if
           if (sumWM > tolWM) then
              call DAXPY (nM,-dX(j,1)/sumWM,work(1),1,omega(j),nndof)
           else if (ipsw > 0) then
              write(lpu,600) iel,lDof,j,sumWM,tolWM
              write(lpu,620) work(1:nM) * dX(j,1)/sumWM
           end if
        else if (ipsw > 0) then
           write(lpu,610) iel,lDof,j,dX(j,1),tolX(j)
        end if
     end if
     if (3+k <= nndof .and. indC(3+k) > 0) then
        iFrst = indC(3+k)
        iLast = iFrst + nM-1
        if (iLast > nW) goto 9
        if (reComputeCG .and. iFrst /= lastIC) then
           lastIC = iFrst
           call computePointCoords (weight(iFrst:iLast))
        end if
        if (abs(dX(j,1)) > tolX(j)) then
           sumWM = sum(weight(iFrst:iLast))
           if (sumWM > epsDiv0_p) then
              call DAXPY (nM,-dX(j,1)/sumWM,weight(iFrst),1,omega(3+k),nndof)
           else if (ipsw > 0) then
              write(lpu,600) iel,lDof,3+k,sumWM,epsDiv0_p
           end if
        else if (ipsw > 0) then
           write(lpu,610) iel,lDof,3+k,dX(j,1),tolX(j)
        end if
     end if

     !! Coupling of translational dependent DOF to translational DOFs i,k
     !! and rotational DOF 3+j through the eccentricity e_k
     !! of the dependent node w.r.t. the independent node centroid

     if (indC(i) /= 0) then
        iFrst = indC(i)
        iLast = iFrst + nM-1
        if (iFrst > 0 .and. iLast > nW) goto 9
        if (reComputeCG .and. iFrst /= lastIC) then
           lastIC = iFrst
           if (iFrst > 0) then
              call computePointCoords (weight(iFrst:iLast))
           else ! uniform weights
              call computePointCoords ()
           end if
        end if
        if (abs(dX(k,1)) > tolX(k)) then
           if (iFrst > 0) then
              sumWM = sumSqrW(dX(k,2:1+nM),dX(i,2:1+nM),weight(iFrst:iLast))
              work(1:nM) = dX(k,2:1+nM)*weight(iFrst:iLast)
           else
              sumWM = sumSqr(dX(k,2:1+nM),dX(i,2:1+nM)) * real(nM,dp)
              call DCOPY (nM,dX(k,2),size(dX,1),work(1),1)
           end if
           tolWM = tolX(k)*tolX(k) + tolX(i)*tolX(i)
           if (sumWM > tolWM) then
              call DAXPY (nM,dX(k,1)/sumWM,work(1),1,omega(i),nndof)
           else if (ipsw > 0) then
              write(lpu,600) iel,lDof,i,sumWM,tolWM
              write(lpu,620) work(1:nM) * dX(k,1)/sumWM
           end if
        else if (ipsw > 0) then
           write(lpu,610) iel,lDof,i,dX(k,1),tolX(k)
        end if
     end if
     if (indC(k) /= 0) then
        iFrst = indC(k)
        iLast = iFrst + nM-1
        if (iFrst > 0 .and. iLast > nW) goto 9
        if (reComputeCG .and. iFrst /= lastIC) then
           lastIC = iFrst
           if (iFrst > 0) then
              call computePointCoords (weight(iFrst:iLast))
           else ! uniform weights
              call computePointCoords ()
           end if
        end if
        if (abs(dX(k,1)) > tolX(k)) then
           if (iFrst > 0) then
              sumWM = sumSqrW(dX(k,2:1+nM),dX(i,2:1+nM),weight(iFrst:iLast))
              work(1:nM) = dX(i,2:1+nM)*weight(iFrst:iLast)
           else
              sumWM = sumSqr(dX(k,2:1+nM),dX(i,2:1+nM)) * real(nM,dp)
              call DCOPY (nM,dX(i,2),size(dX,1),work(1),1)
           end if
           tolWM = tolX(k)*tolX(k) + tolX(i)*tolX(i)
           if (sumWM > tolWM) then
              call DAXPY (nM,-dX(k,1)/sumWM,work(1),1,omega(k),nndof)
           else if (ipsw > 0) then
              write(lpu,600) iel,lDof,k,sumWM,tolWM
              write(lpu,620) work(1:nM) * dX(k,1)/sumWM
           end if
        else if (ipsw > 0) then
           write(lpu,610) iel,lDof,k,dX(k,1),tolX(k)
        end if
     end if
     if (3+j <= nndof .and. indC(3+j) > 0) then
        iFrst = indC(3+j)
        iLast = iFrst + nM-1
        if (iLast > nW) goto 9
        if (reComputeCG .and. iFrst /= lastIC) then
           lastIC = iFrst
           call computePointCoords (weight(iFrst:iLast))
        end if
        if (abs(dX(k,1)) > tolX(k)) then
           sumWM = sum(weight(iFrst:iLast))
           if (sumWM > epsDiv0_p) then
              call DAXPY (nM,dX(k,1)/sumWM,weight(iFrst),1,omega(3+j),nndof)
           else if (ipsw > 0) then
              write(lpu,600) iel,lDof,3+j,sumWM,epsDiv0_p
           end if
        else if (ipsw > 0) then
           write(lpu,610) iel,lDof,3+j,dX(k,1),tolX(k)
        end if
     end if

  end if

5 if (ipsw > 4) then
     write(lpu,*)
     write(label,"('     Omega for dependent DOF',I2,' :')") lDof
     call writeObject(omega,lpu,label,nndof)
  end if
  if (ipsw > 0) call flush (lpu)
  return

9 continue
  write(lpu,690) iFrst, iLast, nW, iel
  goto 5

600 format('  ** Warning: WAVGM element',I10, &
         & ': Ignored coupling by dependent DOF',I2,' to independent DOFs',I2 &
         / 14X, 'due to small denominator',1PE13.5,'  tolerance =',E12.5 )
610 format('  ** Warning: WAVGM element',I10, &
         & ': Ignored coupling by dependent DOF',I2,' to independent DOFs',I2 &
         / 14X, 'due to small eccentricity',1PE13.5,'  tolerance =',E12.5 )
620 format(14X,1P10E13.5/(14X,1P10E13.5))
690 format(' *** Error: Indices out of range:',3I6,/12X,'For WAVGM element',I10)

contains

  !> @brief Computes the nodal coordinates relative to the centre of gravity.
  !> @details If the weights (W) are present, a weighted CoG is used.
  !> This subroutine also recomputes the geometric tolerances (tolX).
  subroutine computePointCoords (W)
    real(dp), intent(in), optional :: W(:)
    integer  :: n
    real(dp) :: X0(3), Ws
    X0 = 0.0_dp
    Ws = 0.0_dp
    if (present(W)) then
       do n = 2, 1+nM
          X0 = X0 + tenc(:,n)*W(n-1)
       end do
       Ws = sum(W)
    else
       do n = 2, 1+nM
          X0 = X0 + tenc(:,n)
       end do
       Ws = real(nM,dp)
    end if
    X0 = X0 / Ws
    do n = 1, 1+nM
       dX(:,n) = tenc(:,n) - X0
    end do
    do n = 1, 3
       tolX(n) = max(maxval(-dX(n,1:nM))*epsX, &
            &        maxval( dX(n,1:nM))*epsX, epsDiv0_p)
    end do
    if (ipsw > 3) then
       if (present(W)) then
          write(lpu,699) tolX, X0, Ws, dX(:,1), (n,dX(:,1+n),W(n),n=1,nM)
       else
          write(lpu,699) tolX, X0, Ws, dX(:,1), (n,dX(:,1+n),1.0_dp,n=1,nM)
       end if
699    format(/25X,'X-coor',7X,'Y-coor',7X,'Z-coor',7X,'Weight', &
            & / 5X,'Tolerance       ',1P3E13.5, &
            & / 5X,'Centroid        ',1P4E13.5, &
            & / 5X,'Reference node  ',1P3E13.5, &
            & /(5X,'Element node', I4,1P4E13.5))
    end if
  end subroutine computePointCoords

  !> @brief Returns the weighted square sum of two arrays.
  !> @details The following sum is returned: @code
  !> Sum_i=1,nM {W(i)*(dX(i)^2+dY(i)^2)}
  !> @endcode
  function sumSqrW (dX,dY,W)
    real(dp),intent(in) :: dX(:), dY(:), W(:)
    real(dp)            :: sumSqrW
    sumSqrW = sum(W*(dX*dX+dY*dY))
  end function sumSqrW

  !> @brief Returns the square sum of two arrays.
  !> @details The following sum is returned: @code
  !> Sum_i=1,nM {dX(i)^2+dY(i)^2}
  !> @endcode
  function sumSqr (dX,dY)
    real(dp),intent(in) :: dX(:), dY(:)
    real(dp)            :: sumSqr
    sumSqr = sum(dX*dX+dY*dY)
  end function sumSqr

end subroutine wavgmConstrEqn
