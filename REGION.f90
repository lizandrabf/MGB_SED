    !---------------------------------------------------------------------------------
    !  Discussion:
    !
    !    This sub-routine calculates the geomorphologic relationships of depth (Inercial model) and width (Muskingum-Cunge and Inercial model)
    !    as Drainage Area function.
    !
    ! Usage:
    !
    !    * no modules, functions, or subroutines are used in this funcation
    !
    ! uses modules and functions
    !
    !    * module VARS_MAIN in VARSMAIN.f90
    !    * module VARS_INERC (Inertial Version)
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
    ! 2015.06.21 - 21 June 2015 (By: Fernando Mainardi Fan) 
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
    ! Variables and Parameters:
    ! *Variables declarations and routines calls are all commented below.
    !---------------------------------------------------------------------------------
    ! End of header
    !---------------------------------------------------------------------------------
    
	SUBROUTINE REGION

    ! Variables and Parameters
	USE VARS_MAIN
	USE VARS_INERC
	IMPLICIT NONE

	REAL QMESP !Specific flow L/S/KM2

!	QMESP=17. !Specific flow (L/S/KM2) avaliable to Taquari-Antas river 
	
	QMESP=30. !Specific flow (L/S/KM2) avaliable to Amazon river
	
	QRX=QMESP*ARX/1000. 
!    WRITE(*,*) 'QMESP, ARX, QRX = ', QMESP, ARX, QRX
!    PAUSE
	!geomorphologic relationships (width) avaliable to Taquari-Antas river 
	!BRX=5.2466*ARX**0.4106

	!geomorphologic relationships (width) avaliable to Amazon river 
	BRX=0.8054*ARX**0.5289

    !geomorphologic relationships (depth) avaliable to Amazon river (only used in Inertial model)
    HRX=0.9*ARX**0.19 
    !HRX=1.4351*ARX**0.1901 !FOR PURUS RIVER
    
    !geomorphologic relationships (width) avaliable to Paraná river 
	!BRX=1.0391*ARX**0.48
    !BRX=0.8054*ARX**0.5289 !for PURUS RIVER

	
	

	RETURN
	END
