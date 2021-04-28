	SUBROUTINE ALLOCA_VARS_SED(IOP)

!---------------------------------------------------------------------------------
! ALLOCA_VARS_SED.f90
!---------------------------------------------------------------------------------
! Discussion:
!
! This routine allocate or deallocate the MGB-SED memory for the main variables
!
! Usage:
!
! CALL ALLOCA_VARS_SED(IOP)
!
! uses modules, functions, and subroutines
!
! *Module VARS_MAIN !module with the model main variables
! *Module SED_VARS
!
! opens
!
! * no files are opened in this routine
!
! reads
!
! * no files are read in this routine
!
! creates
!
! * no files are created in this routine
!
!---------------------------------------------------------------------------------
! Licensing:
!
! This code is distributed under the...
!
! Created at 02 July 2019

! By: Hugo Fagundes
! Based on ALLOCA_VARS (Fernando Fan)
!
!---------------------------------------------------------------------------------
! Variables and Parameters:
! *Variables declarations and routines calls are all commented below.
!---------------------------------------------------------------------------------
! End of header
!---------------------------------------------------------------------------------

! Variables and Modules

	USE VARS_MAIN
    USE SED_VARS !Hugo Fagundes 16/07/2019
	IMPLICIT NONE	
	INTEGER IOP
    SAVE

! Select the case how the subroutine is working: for allocation, or for deallocation.

	ALLOCV_CASE: SELECT CASE (IOP) !Checks for allocation or deallocation. IOP=0 means allocation, and IOP=1 means deallocation.
    CASE (0) ! Allocation         
        ALLOCATE (SSOBS(NMAXHID,NT,AUXSS)) !OBSERVATIONS RELATED TO SEDIMENT DATA
        ALLOCATE (AUXOBS(1000*NMAXHID),AUXCALC(1000*NMAXHID)) !STORAGE AVERAGE FOR EACH STATION TO COMPUTE SPATIAL CORRELATION 
        ALLOCATE (SSCORRELtemp(NMAXHID), SSCORRELesp(NMAXHID), SSCORRELtudo(NMAXHID)) !Temporal, spatial and global correlation
        ALLOCATE (RMSE(NMAXHID),CoefB(NMAXHID),DCPERM (NMAXHID), R2SS(NMAXHID),ERRVS(NMAXHID),SSCORRELAC(NMAXHID),SSKGE(NMAXHID),KGEa(NMAXHID))  !RMSE, Beta Coefficient, Slope Sed. duration curve, NASH SS, VOL. SS, CORREL TEMPORAL, KLING GUPTA, alfa Kling Gupta !Hugo Fagundes
        ALLOCATE (ALFsed(NB),BETsed(NB),AuxTKS(NB)) ! MUSLE EQUATION COEFF, AND AUXILIAR COEF. RELATED TO DELAY TIME OF LINEAR RESERVOIR
        ALLOCATE (CSSbase(NB,NB)) ! BACKGROUND SUSPENDED SEDIMENT CONCENTRATION
        ALLOCATE (CARGAaux(NC)) !Carga auxiliar HUGO FAGUNDES 14/08/2019
        !ALLOCATE (RESULTSS(1000,200)) !RECORD FINAL RESULTS OF VALUES FROM OBJETCTIVE FUNCTIONS FOR EACH DATASET (CALIBRATION MODULE)      
    CASE (1) ! Deallocation
	    DEALLOCATE (SSOBS) !OBSERVATIONS RELATED TO SEDIMENT DATA
        DEALLOCATE (AUXOBS,AUXCALC) !STORAGE AVERAGE FOR EACH STATION TO COMPUTE SPATIAL CORRELATION 
        DEALLOCATE (R2SS,ERRVS,SSCORRELAC,SSKGE) !NASH SS, VOL. SS, SEDIMENT r Correl, Sediment KLING GUPTA
        DEALLOCATE (SSCORRELtemp, SSCORRELesp, SSCORRELtudo) !Temporal, spatial and global correlation
        DEALLOCATE (ALFsed,BETsed,AuxTKS) ! MUSLE EQUATION COEFF, AND AUXILIAR COEF. RELATED TO DELAY TIME OF LINEAR RESERVOIR
        DEALLOCATE (CSSbase) ! BACKGROUND SUSPENDED SEDIMENT CONCENTRATION
        DEALLOCATE (NOBSS,NUMSS) !NUMBER SEDIMENT GAUGING STATIONS
        DEALLOCATE (ISSOBS,ISS) !STORAGE ID CODES OF CATCHMENTS FOR EACH GAUGE
        DEALLOCATE (Kusle, Cusle, Pusle, Rgros, Ksdr) !MUSLE FACTORS    (SEE SED_PARAM ROUTINE)
        DEALLOCATE (Mareia, Msilte, Margila, Morg, Mrocha) !SOIL TEXTURE VARIABLES  (SEE SED_PARAM ROUTINE)
        DEALLOCATE (LSAcu, NPXU, SDR) !VARIABLES TO COMPUTE MUSLE EQUATION  (SEE SED_PARAM ROUTINE)
        DEALLOCATE (PSC, SLC, QSC, QSAUX) !VARIABLES RELATED TO MUSLE EQUATION  (SEE SED_PARAM ROUTINE)
        DEALLOCATE (SLU)!VARIABLES RELATED TO MUSLE EQUATION    (SEE SED_PARAM ROUTINE)
        DEALLOCATE (QSU, PSU, SLXU)!VARIABLES RELATED TO MUSLE EQUATION (SEE SED_PARAM ROUTINE)
        DEALLOCATE (CAREIA, CSILTE, CARGIL) !CONCENTRATIONS THAT WILL BE WRITTEN
        DEALLOCATE (CSSCAL,SEDstore)  !SUSPENDED SEDIMENT CONCENTRATION IN OUTLET OF UNIT CATCHMENTS
        !DEALLOCATE (RESULTSS) !RECORD FINAL RESULTS OF VALUES FROM OBJETCTIVE FUNCTIONS FOR EACH DATASET (CALIBRATION MODULE)       
        
        
	CASE DEFAULT
		STOP ' ERROR: IOP UNKNOWN AT ROUTINE ALLOCA_CALIB!!!'
	END SELECT ALLOCV_CASE

	! End of the routine.
	END
