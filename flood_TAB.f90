	!---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This subroutine generates a level x volume measured from COTA_AREA.flp data  (Inertial version).
    !
    !
    ! Usage:
    !
    ! *
    !
    ! uses modules, functions, and subroutines
    !
    ! * USE VARS_MAIN
    ! * USE VARS_INERC (only to Inertial version)
    !
    ! opens
    !
    ! * no files are created in this routine
    !
    ! reads
    !
    ! * no files are created in this routine
    !
    ! creates
    !
    ! * no files are created in this routine
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
    ! 2015.07.06 - 07 July 2015 (by Paulo Pontes)
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
    
    SUBROUTINE flood_TAB

    !Variables and parameters
    USE VARS_MAIN
    USE VARS_INERC
    IMPLICIT NONE
    INTEGER K
    REAL VOLINCR

    DO IC=1,NC
        
        !From 1 to 3 (points), it is calculated the area and volume of river cross-section and updates the units to meters
        ZTAB(1,IC)=ZFUNDOFL(IC)-HRIO(IC) 
        VTAB(1,IC)=0.0 
        ATAB(1,IC)=0.0
        ZTAB(2,IC)=ZFUNDOFL(IC) 
        VTAB(2,IC)=DBLE(BRIO(IC))*HRIO(IC)*DBLE(SRIO(IC))*1000.0 
        ATAB(2,IC)=0.0
        ZTAB(3,IC)=ZFL(1,IC) 
        VTAB(3,IC)=VTAB(2,IC)+0.5*(AFL(1,IC)+0.0)*(ZTAB(3,IC)-ZTAB(2,IC))*1000000.0
        ATAB(3,IC)=AFL(1,IC)
        
        !From next points the volume are calculated using area and level from trapezium approach
        DO K=1,NPFL(IC)-1
            ZTAB(K+3,IC)=ZFL(K+1,IC)
            VOLINCR=(ZFL(K+1,IC)-ZFL(K,IC))*1000000.0*((AFL(K,IC)+AFL(K+1,IC))*0.5) 
            VTAB(K+3,IC)=VTAB(K+2,IC)+VOLINCR
            ATAB(K+3,IC)=AFL(K+1,IC)
        ENDDO
        
    ENDDO  
    
    RETURN
    END
