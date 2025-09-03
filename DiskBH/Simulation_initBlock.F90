!!****if* source/Simulation/SimulationMain/DiskNSRefRot/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(IN) :: blockID) 
!!                       
!!
!!
!!
!! DESCRIPTION
!!
!!  Disco de acreção com velocidade angular Kepleriana
!!
!!
!! 
!! ARGUMENTS
!!
!! PARAMETERS
!!  sim_mas              Massa da Estrela
!!  sim_temp             Temperatura do disco
!!  sim_dotmas           variação de massa
!!  sim_k		Constante de Bolztmann
!!  sim_mash             Massa do hidrogenio
!!  sim_g		Constante gravitaconal 
!!  sim_mu		Peso molecular padrão
!!  ptmass	 	Point mass if external field
!!
!!
!!***

subroutine Simulation_initBlock(blockID)

  use Simulation_data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getCellCoords, Grid_putPointData


  implicit none

#include "constants.h"
#include "Flash.h"

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)
  
  integer, intent(in) :: blockID
  

  integer :: i, j, k, n
  integer :: iMax, jMax, kMax
  


  real :: xx, yy,  zz, xxL, xxR
  
  real :: lPosn0, lPosn
  

  real,allocatable, dimension(:) ::xCenter,xLeft,xRight,yCoord,zCoord

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  integer, dimension(MDIM) :: axis

  
  real :: rhoZone, velxZone, velyZone, velzZone, presZone, & 
       eintZone, enerZone, ekinZone
  
  logical :: gcell = .true.
  integer :: istat 

  
  ! dump some output to stdout listing the paramters
  if (sim_meshMe == MASTER_PE) then
     
     
1    format (1X, 1P, 4(A7, E13.7, :, 1X))
2    format (1X, 1P, 2(A7, E13.7, 1X), A7, I13)
     
  endif
  
  
  ! get the integer index information for the current block
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  
  sizeX = blkLimitsGC(HIGH,IAXIS)
  sizeY = blkLimitsGC(HIGH,JAXIS)
  sizeZ = blkLimitsGC(HIGH,KAXIS)
  allocate(xLeft(sizeX),stat=istat)
  allocate(xRight(sizeX),stat=istat)
  allocate(xCenter(sizeX),stat=istat)
  allocate(yCoord(sizeY),stat=istat)
  allocate(zCoord(sizeZ),stat=istat)
  xCenter = 0.0
  xLeft = 0.0
  xRight = 0.0
  yCoord = 0.0
  zCoord = 0.0

  if (NDIM == 3) call Grid_getCellCoords&
                      (KAXIS, blockId, CENTER,gcell, zCoord, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords&
                      (JAXIS, blockId, CENTER,gcell, yCoord, sizeY)

  call Grid_getCellCoords(IAXIS, blockId, LEFT_EDGE, gcell, xLeft, sizeX)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCenter, sizeX)
  call Grid_getCellCoords(IAXIS, blockId, RIGHT_EDGE, gcell, xRight, sizeX)

!------------------------------------------------------------------------------

! Loop over cells in the block.  For each, compute the physical position of 
! its left and right edge and its center as well as its physical width.  
! Then decide which side of the initial discontinuity it is on and initialize 
! the hydro variables appropriately.


  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     
     ! get the coordinates of the cell center in the z-direction
     zz = zCoord(k)     
     
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        
        ! get the coordinates of the cell center in the y-direction
        yy = yCoord(j)
        
        ! The position of the shock in the current yz-row.
        
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           
           ! get the cell center, left, and right positions in x
           xx  = xCenter(i)
           xxL = xLeft(i)
	   xxR = xRight(i)

	      velxZone = -0.01*sqrt(sim_mas*sim_g/xx)
              velyZone = sqrt(sim_mas*sim_g/xx)
              velzZone = 0.0

	      rhoZone  = sim_dotmas / (4.0*3.1416*1.0E-2*xx*xx*0.01*sqrt(sim_mas*sim_g/xx))

	      presZone = (rhoZone * sim_k * sim_temp) / (sim_mu * sim_mash)

           axis(IAXIS) = i
           axis(JAXIS) = j
           axis(KAXIS) = k 

           !put in value of default species
           if (NSPECIES > 0) then
              call Grid_putPointData(blockID, CENTER, SPECIES_BEGIN, EXTERIOR, &
                   axis, 1.0e0-(NSPECIES-1)*sim_smallX)


              !if there is only 1 species, this loop will not execute
              do n = SPECIES_BEGIN+1,SPECIES_END
                 call Grid_putPointData(blockID, CENTER, n, EXTERIOR, &
                      axis, sim_smallX)
              enddo
           end if

           ! Compute the gas energy and set the gamma-values needed for the equation of 
           ! state.
           ekinZone = 0.5 * (velxZone**2 + & 
                velyZone**2 + & 
                velzZone**2)
           
           eintZone = presZone / (sim_gamma-1.)
           eintZone = eintZone / rhoZone
           enerZone = eintZone + ekinZone
           enerZone = max(enerZone, sim_smallP)
           
           ! store the variables in the current zone via Grid put methods
           ! data is put stored one cell at a time with these calls to Grid_putData           


           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rhoZone)
           call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, presZone)
           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, velxZone)
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, velyZone)
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, velzZone)

#ifdef ENER_VAR
           call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, axis, enerZone)   
#endif
#ifdef EINT_VAR
           call Grid_putPointData(blockId, CENTER, EINT_VAR, EXTERIOR, axis, eintZone)   
#endif
#ifdef GAME_VAR          
           call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, axis, sim_gamma)
#endif
#ifdef GAMC_VAR
           call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, axis, sim_gamma)
#endif


        enddo
     enddo
  enddo

!! Cleanup!  Must deallocate arrays

  deallocate(xLeft)
  deallocate(xRight)
  deallocate(xCenter)
  deallocate(yCoord)
  deallocate(zCoord)


 
  return
end subroutine Simulation_initBlock










