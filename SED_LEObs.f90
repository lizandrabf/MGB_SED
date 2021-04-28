    !*********************************************************************************
    !
    !  SUBROUTINE LEOBSed reads the files with observed sediment data (File with extension .txt)
    !
    !---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This routine reads the files with observed sediment data (File with extension .txt)

	!
	!	 LEOBSed is called inside 1main.
	!
	!	 Saves sediment time series: 
	!     	QOBS(.,.)
	!
    !
    !  	Usage:
    !
    !    * no subroutine is called in this subroutine
    !
    !    where
    !
    !    * no arguments are passed in this subroutine
    !
    !    uses modules and functions
    !
    !    * module     VARS_MAIN   in      VARS_MAIN.f90
	!
    !	 opens
    !
    !    * Files with extension .txt containing time series of observed sediment data. 
    !
    !    reads
    !
    !    * File with extension .txt containing time series of observed SSC, turbidity, TSS and surface reflectance
    !
    !    creates
    !
    !    * Does not create files
    !    
    !
    !---------------------------------------------------------------------------------
    !  Licensing:
    !
    !    This code is distributed under the ...
    !
    !  Version/Modified: 
    !
    !    22 September 2017 (By: Hugo Fagundes)    
    !
    !  Authors:
    !
    !    Original fortran version by Walter Collischonn
    !    Present fortran version by:
    !    * Walter Collischonn
    !    * Rodrigo Cauduro Dias de Paiva
    !    * Diogo da Costa Buarque
    !    * Paulo Pontes Rógenes
    !    * Mino  Viana Sorribas
    !    * Fernando Mainardi Fan
    !    * Juan Martin Bravo 
    !    * Hugo de Oliveira Fagundes
    !
    !  Main Reference:
    !
    !    Walter Collischonn,
    !    Modelo de Grandes Bacias - Thesis
    !    Porto Alegre, 2001
    !    ISBN: XXXXXXXXXXX,
    !
    !---------------------------------------------------------------------------------
    !  Variables and Parameters:
    !
    !   *Variables declarations and routines calls are all commented below.
	!	* All variables are global!?
    !
    !---------------------------------------------------------------------------------		
    
    
    SUBROUTINE LEOBSed

	USE VARS_MAIN
    USE SED_VARS
	IMPLICIT NONE
	INTEGER I,J,K,L,KB,KSS,AUX

	DO KSS=1,AUXSS
        AUX=NOBSS(KSS)
        OPEN(FILSS,FILE=INPUT_DIRECTORY //ARQOBSSS(KSS),STATUS='OLD')
	    READ(FILSS,704)(CABE(K),K=1,AUX)
	    DO IT=1,NT
           READ(FILSS,*)I,J,L,(SSOBS(KB,IT,KSS),KB=1,AUX)    
        ENDDO
        CLOSE (FILSS)
        704 FORMAT(<AUX>A10)
    ENDDO

	RETURN
	END