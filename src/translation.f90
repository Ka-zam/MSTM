module translation
   use translation_expansions, only: clear_stored_trans_mat, external_to_external_expansion, &
                                     external_to_internal_expansion, general_interaction_matrix, interaction_radius
   use translation_operator, only: shiftcoefficient, translation_data
   use translation_surface_interactions, only: periodic_lattice_sphere_interaction, spheresurfaceinteraction

   implicit none(type, external)
   private

   public :: translation_data
   public :: interaction_radius
   public :: clear_stored_trans_mat
   public :: general_interaction_matrix
   public :: periodic_lattice_sphere_interaction
   public :: spheresurfaceinteraction
   public :: external_to_external_expansion
   public :: external_to_internal_expansion
   public :: shiftcoefficient

end module translation
