!!****if* source/Simulation/SimulationMain/DiscoNSRefRot/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data
!!
!! DESCRIPTION
!!
!!  Henrique Hirsch Disco SB
!!
!! ARGUMENTS
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
!!
!!***

module Simulation_data

  implicit none

  !! *** Runtime Parameters *** !!

  real, save :: sim_mas, sim_dotmas, sim_mash
  real, save :: sim_temp
  real, save :: sim_k
  real, save :: sim_g
  real, save :: sim_mu
  real, save :: sim_gamma, sim_smallP, sim_smallX
  real, save :: sim_xmin,sim_xmax,sim_ymin,sim_ymax

  !! *** Variables pertaining to Simulation Setup  *** !!
  logical, save :: sim_gCell

integer, save :: sim_meshMe
end module Simulation_data


