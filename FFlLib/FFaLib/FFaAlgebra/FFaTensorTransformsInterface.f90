!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module FFaTensorTransformsInterface

  implicit none

  interface

     function vonmises (N,S)
       integer , parameter  :: dp = kind(1.0D0)
       real(dp)             :: vonmises
       integer , intent(in) :: N
       real(dp), intent(in) :: S(*)
     end function vonmises

     subroutine princval (N,S,Pval)
       integer , parameter   :: dp = kind(1.0D0)
       integer , intent(in)  :: N
       real(dp), intent(in)  :: S(*)
       real(dp), intent(out) :: Pval(*)
     end subroutine princval

     subroutine maxshearvalue (N,Pval,maxShear)
       integer , parameter   :: dp = kind(1.0D0)
       integer , intent(in)  :: N
       real(dp), intent(in)  :: Pval(*)
       real(dp), intent(out) :: maxShear
     end subroutine maxshearvalue

     subroutine maxshear (N,Pval,Pdir,Shear,Sdir)
       integer , parameter   :: dp = kind(1.0D0)
       integer , intent(in)  :: N
       real(dp), intent(in)  :: Pval(*), Pdir(*)
       real(dp), intent(out) :: Shear, Sdir(*)
     end subroutine maxshear

     subroutine tratensor (N,S,Tmat)
       integer , parameter     :: dp = kind(1.0D0)
       integer , intent(in)    :: N
       real(dp), intent(inout) :: S(*)
       real(dp), intent(in)    :: Tmat(*)
     end subroutine tratensor

     subroutine trainertia (N,I,Xvec,Mass)
       integer , parameter     :: dp = kind(1.0D0)
       integer , intent(in)    :: N
       real(dp), intent(inout) :: I(*)
       real(dp), intent(in)    :: XVec(*), Mass
     end subroutine trainertia

  end interface

end module FFaTensorTransformsInterface
