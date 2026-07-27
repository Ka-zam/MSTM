module input_execution
   use simulation_averaging, only: run_configuration_average, run_incidence_average, &
                                   run_random_orientation_configuration_average
   use simulation_execution, only: execute_simulation
   implicit none(type, external)
   private

   public :: execute_simulation, run_configuration_average, run_incidence_average, &
             run_random_orientation_configuration_average
end module input_execution
