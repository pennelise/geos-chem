! $Id: CO_strat_pl.f,v 1.3 2005/03/29 15:52:38 bmy Exp $
      SUBROUTINE CO_STRAT_PL( COPROD, COLOSS )
!
!******************************************************************************
!  Subroutine CO_STRAT_PL computes net production of CO above the 
!  annual mean tropopause using archived rates for P(CO) and L(CO).
!  (bnd, qli, bmy, 12/9/99, 7/20/04)
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) COPROD : (REAL*4) Zonally averaged P(CO) in [v/v/s]
!  (2 ) COLOSS : (REAL*4) Zonally averaged L(CO) in [1/s]
!
!  NOTES:
!  (1 ) P(CO) and L(CO) rates were provided by Dylan Jones.  Bob Yantosca
!        has regridded these rates to the GEOS-1 and GEOS-STRAT vertical 
!        levels, for both 2 x 2.5 and 4 x 5 degree resolution.
!  (2 ) CO_STRAT_PL was specially adapted from Bryan Duncan's routine 
!        "CO_strat.f" for use with the simple chemistry module "schem.f".
!  (3 ) Now reference AD from "dao_mod.f".  Now reference IDTCO and IDTCH2O
!        from "tracerid_mod.f". (bmy, 11/6/02)
!  (4 ) Now use function GET_TS_CHEM from "time_mod.f".  Updated comments.
!        (bmy, 2/11/03)
!  (5 ) Now references STT from "tracer_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules 
      USE DAO_MOD,      ONLY : AD
      USE TIME_MOD,     ONLY : GET_TS_CHEM
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD, ONLY : IDTCO, IDTCH2O

      IMPLICIT NONE
     
#     include "CMN_SIZE"   ! Size parameters
#     include "CMN"        ! LPAUSE
#     include "CMN_O3"     ! XNUMOLAIR

      ! Arguments
      REAL*4,  INTENT(IN) :: COPROD(JJPAR,LLPAR)
      REAL*4,  INTENT(IN) :: COLOSS(JJPAR,LLPAR)

      ! Local variables
      INTEGER             :: I, J, L, M, N, LMIN

      REAL*8              :: BAIRDENS, DT, GCO, STTTOGCO

      ! External functions
      REAL*8, EXTERNAL    :: BOXVL

      !=================================================================
      ! CO_STRAT_PL begins here!
      !=================================================================

      ! Chemistry timestep [s]
      DT = GET_TS_CHEM() * 60d0

      !=================================================================
      ! Loop over all stratospheric grid boxes ( L >= LPAUSE(I,J) ). 
      !
      ! Compute the net CO from the P(CO) and L(CO) rates that are 
      ! stored in the COPROD and COLOSS arrays.
      !
      ! Unit conversion to/from [kg/box] and [molec/cm3] is required.
      ! The conversion factor is STTTOGCO, which is given below.
      !
      !   kg CO       box     |   mole CO   | 6.022e23 molec CO       
      !  ------- * -----------+-------------+-------------------  
      !    box      BOXVL cm3 | 28e-3 kg CO |     mole CO             
      !
      !  =  molec CO
      !     --------
      !       cm3
      !=================================================================
      DO L = MINVAL( LPAUSE ), LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Skip tropospheric grid boxes
         IF ( L < LPAUSE(I,J) ) CYCLE
         !IF ( ITS_IN_THE_TROP(I,J,L) ) CYCLE

         ! conversion factor from [kg/box] to [molec/cm3]
         STTTOGCO = 6.022d23 / ( 28d-3 * BOXVL(I,J,L) )

         ! Convert STT from [kg/box] to [molec/cm3]
         GCO = STT(I,J,L,IDTCO) * STTTOGCO

         ! Air density in molec/cm3
         BAIRDENS = AD(I,J,L) * XNUMOLAIR / BOXVL(I,J,L)

         ! Apply P(CO) and L(CO) rates to GCO
         GCO = GCO * ( 1d0 - COLOSS(J,L) * DT ) +
     &               ( COPROD(J,L) * DT * BAIRDENS )

         ! Compute production of CH2O (qli, 12/9/99)
         STT(I,J,L,IDTCH2O) = STT(I,J,L,IDTCH2O) + 
     &                        XNUMOL(IDTCO) / XNUMOL(IDTCH2O) * 
     &                        COPROD(J,L)   * BAIRDENS        / 
     &                        STTTOGCO

         ! Convert STT from [molec/cm3] to [kg/box]
         STT(I,J,L,IDTCO) = GCO / STTTOGCO
      ENDDO
      ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE CO_STRAT_PL
