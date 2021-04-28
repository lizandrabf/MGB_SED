	!---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This subroutine calculates the time-step of inertial model  (Inertial version).
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
	
	subroutine flood_timestep
	!--------------------------------------------------------------------------------------
	! Variables and parameters:
	use VARS_INERC
	use VARS_MAIN
	implicit none
	
	real*8 dtminIC
	real*8 dtteste
	!-------------------------------------------------------------------------------------
    
    dtteste=DBLE(DTP)
    


    DO IC=1,NC  !Catchment loop
        ! MAX H to calculates the time-step of Inertial model:
        hmaxfl=(Hfl(IC))
        ! Correct if equal to zero:
        hmaxfl=max(hmaxfl,0.001)
        ! Compute time interval:
        dtflood=alpha*SRIO(IC)*1000./(g*hmaxfl)**0.5
        if(dtflood<dtteste) dtteste=dtflood
    ENDDO

	
    ! MAX time-step of MGB input data (e.g. dialy):
    dtfloodmax=DBLE(DTP)
    dtminIC=dtteste
    dtflood=min(dtfloodmax,dtminIC)   !dt for all catchments calculated as function of H. 
 	
 	endsubroutine