	!---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This subroutine is the main subroutine of Inertial model  (Inertial version).
    !
    !
    ! Usage: flood_timestep, flood_continuity and flood_discharge
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

    SUBROUTINE flood_inercial
    
    ! Variables and parameters
    USE VARS_MAIN
    USE VARS_INERC
    use IFPORT
    
    implicit none
    REAL(8) elapsed_time 
    integer K !auxiliary counter
    
    
    !-------------------------------------------------------------------------------------
    tflood=0.0
    AFLTUDO(1,IT)=0.0 !inicializa a variável de areas inundadas da bacia
 
    do while (tflood<dtflood0)
        call flood_timestep
        dtflood=min(dtflood,dtflood0-tflood)
        tflood=tflood+dtflood
        ! Inertial equation:
        call flood_discharge
        ! Continuity equation:
        call flood_continuity 
          
    enddo
    QJ2=Q2FL
    QTUDO(:,IT)=Q2fl
    YTUDO(IT,:)=Yfl
    
	return
	endsubroutine
    
    
    