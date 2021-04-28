    !*********************************************************************************
    !
    !  SUBROUTINE LESUBST reads  the file with time series of stream flow to substitute simulation at selected points.  (File: QSUBST.QSB)
    !
    !---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This routine reads the file with time series of stream flow to substitute simulation at selected points.  (File: QSUBST.QSB)

	!
	!	 LESUBST is called inside 1main.
	!
	!	 Saves discharge time series: 
	!     	QLIDO(.,.)
	!
	!
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
    !    * File with extension .QOB containing time series of observed stream flow. 
    !
    !    reads
    !
    !    * File with extension .QOB containing time series of observed stream flow
    !
    !    creates
    !
    !    * Does not create files
    !    
    !
    !---------------------------------------------------------------------------------
    !  Licensing:
    !
    !   This program is free software: you can redistribute it and/or modify
    !   it under the terms of the GNU General Public License as published by
    !   the Free Software Foundation, either version 3 of the License, or
    !   (at your option) any later version.
    !
    !   This program is distributed in the hope that it will be useful,
    !   but WITHOUT ANY WARRANTY; without even the implied warranty of
    !   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    !   GNU General Public License for more details.
    !
    !   You should have received a copy of the GNU General Public License
    !   along with this program.  If not, see <http://www.gnu.org/licenses/>
    !
    !  Version/Modified: 
    !
    !    2014.26.11 - 25 November 2014 (By: Rodrigo Paiva)    
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

	SUBROUTINE LESUBST
	!Reads discharge time series to be substituted by simulation at selected points. 
	USE VARS_MAIN
	IMPLICIT NONE
	INTEGER I,J,K,L
	OPEN(FILSUBS,FILE=INPUT_DIRECTORY // ''//ARQSUBST,STATUS='OLD')
	READ(FILSUBS,701)(CABE(K),K=1,NUMSUBST)
	DO IT=1,NT
	! READ(FILSUBS,702)CABE(1),(QLIDO(K,IT),K=1,NUMSUBST)
	READ(FILSUBS,*)I,J,L,(QLIDO(K,IT),K=1,NUMSUBST)
	ENDDO
	CLOSE (FILSUBS)
701	FORMAT(20A10)
702	FORMAT(A20,20F10.2)
	RETURN
	END
