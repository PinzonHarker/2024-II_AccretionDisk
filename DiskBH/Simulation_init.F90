!!****if* source/Simulation/SimulationMain/DiscoNSRefRot/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!!
!! DESCRIPTION
!!
!!  Disco de acreção SB
!!
!! ARGUMENTS
!!
!!  
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
!!***

subroutine Simulation_init()
  
  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stamp
  use Driver_interface, ONLY : Driver_getMype
  implicit none
#include "Flash.h"
#include "constants.h"
  

  call RuntimeParameters_get('smallp', sim_smallP)
  call RuntimeParameters_get('smallx', sim_smallX) 
  
  call RuntimeParameters_get('gamma', sim_gamma)
  
  call RuntimeParameters_get('sim_temp', sim_temp)
  
   call RuntimeParameters_get('sim_mas', sim_mas)

  call RuntimeParameters_get('sim_dotmas', sim_dotmas)
  
  call RuntimeParameters_get('sim_k', sim_k)

  call RuntimeParameters_get('sim_mash', sim_mash)

  call RuntimeParameters_get('sim_g', sim_g)
  
  call RuntimeParameters_get('sim_mu', sim_mu)

  call RuntimeParameters_get('xmin',      sim_xmin)
  call RuntimeParameters_get('xmax',      sim_xmax)
  call RuntimeParameters_get('ymin',      sim_ymin)
  call RuntimeParameters_get('ymax',      sim_ymax)
  call Driver_getMype(MESH_COMM, sim_meshMe)


  call Logfile_stamp( "Disk Reflect + Rotation",  &
       "[Simulation_init]")


end subroutine Simulation_init







