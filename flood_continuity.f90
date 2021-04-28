	!---------------------------------------------------------------------------------
    !  Discussion:
    ! 
    !    This subroutine calculates the level and depth for each catchment from Continuity equation   (Inertial version).
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

	subroutine flood_continuity

	!--------------------------------------------------------------------------------------
	! Variables and parameters:
	use VARS_INERC
	use VARS_MAIN
	implicit none
	
	real*8 areajtab, yjtab

    !-------------------------------------------------------------------------------------
    do iC=1,nC
    
        IB=IBAC(IC) 
		!IF((IB<SUBini).OR.(IB>SUBfim))CYCLE ! Manual controls for especialized calibration
    
            !Number of upstream catchments of IC
            Nentradas = MINIMONT(iC,1)
            
            !Sum of downstream IC flows
            if(Nentradas==0)then
                SumQup=0.0
            else
                SumQup = SUM(Q2fl(MINIMONT(iC,2:1+Nentradas)))
            endif
            !to use lateral connections
            !Vol2(iC) = Vol1(iC)+dtflood*(SumQup+QCEL2(iC)-Q2fl(iC)+Q2viz(IC))-(EVQ(IC)*dtflood*Area2(iC)*1000000.0/(DTP*1000.))+(P(IC)*dtflood*Area2(iC)*1000000.0/(DTP*1000.))
            
            !without lateral connections in the continuity equation (just remove Q2viz(IC) from above equation
            Vol2(iC) = Vol1(iC)+dtflood*(SumQup+QCEL2(iC)-Q2fl(iC))-(E0agua(IC)*dtflood*Area2(iC)*1000000.0/(DTP*1000.))+(P(IC)*dtflood*Area2(iC)*1000000.0/(DTP*1000.))
            Vol2(iC) = max(Vol2(iC),0.0)
            
            !Interpolates the Area and Level from Volume
            CALL hunt(VTAB(:,iC),Vol2(iC),jtab(iC),ATAB(:,iC),Areajtab,ZTAB(:,IC),yjtab,NPFL(IC)+2)
            Area2(iC) = min(Areajtab,ACEL(iC))
            
            !Updates variables:
            y2_fl=max(yjtab,ZTAB(1,IC))
            y2_fl=y2_fl+0.001
            
            !Calculates depth
            Hfl(iC)=y2_fl-ZTAB(1,IC)
            Yfl(iC)=y2_fl       
            
            !Updates the volume
            Vol1(iC)=Vol2(iC)  
            
            !Saves flooded area for the whole basin for this catchment
            
            if (tflood/dtflood0==1) AFLTUDO(1,IT)=AFLTUDO(1,IT)+Area2(IC) 

    enddo
	
	endsubroutine