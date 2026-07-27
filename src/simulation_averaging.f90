module simulation_averaging
   use configuration_averaging, only: run_configuration_average
   use incidence_averaging, only: run_incidence_average
   use random_orientation_averaging, only: run_random_orientation_configuration_average
   implicit none(type, external)
   private

   public :: run_configuration_average, run_incidence_average, &
             run_random_orientation_configuration_average
end module simulation_averaging
