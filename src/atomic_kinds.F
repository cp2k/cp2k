MODULE atomic_kinds

  USE basis_set_types,                 ONLY: gto_basis_set_type
  USE external_potential_types,        ONLY: all_potential_type,&
                                             gth_potential_type
  USE kinds,                           ONLY: dp

  IMPLICIT NONE

  PRIVATE

! *** Parameters ***

  INTEGER, PARAMETER :: max_kind_label_length = 8,&
                        max_basis_set_name_length = 60

! *** Define the atomic kind info type ***

  TYPE kind_info_type
    TYPE(all_potential_type), POINTER        :: all_potential
    TYPE(gth_potential_type), POINTER        :: gth_potential
    TYPE(gto_basis_set_type), POINTER        :: aux_basis_set
    TYPE(gto_basis_set_type), POINTER        :: orb_basis_set
    CHARACTER(LEN=max_kind_label_length)     :: label
    CHARACTER(LEN=max_basis_set_name_length) :: aux_basis_set_name,&
                                                orb_basis_set_name,&
                                                potential_name
    CHARACTER(LEN=2)                         :: element_symbol
    REAL(KIND = dp)                                 :: alpha_core_charge,&
                                                ccore_charge,&
                                                core_charge_radius,mass,&
                                                prj_paw_radius
    INTEGER                                  :: natom,z,zeff
    INTEGER                                  :: number_of_grid_points
    INTEGER, DIMENSION(:), POINTER           :: atom_list,elec_conf
  END TYPE kind_info_type

  TYPE(kind_info_type), DIMENSION(:), POINTER :: kind_info

! *** Public data types ***

  PUBLIC :: kind_info_type,kind_info

! *****************************************************************************

END MODULE atomic_kinds
