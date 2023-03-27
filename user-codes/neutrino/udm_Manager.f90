module udm_Manager
use udm_Parameter

!=======================================================================
! [udm_int]
use udm_int_neutrino,   caller_udm_int_neutrino  => caller
! [udm_part]
! use udm_part_MuonPlus,  caller_udm_part_MuonPlus => caller
!=======================================================================

implicit none
contains




subroutine user_defined_interaction(action,index)
integer action,index
!=======================================================================
call caller_udm_int_neutrino(action,index)
!=======================================================================
end subroutine user_defined_interaction




subroutine user_defined_particle(action)
integer action,index
do index=1,udm_part_nMax
!=======================================================================
! call caller_udm_part_MuonPlus(action,index)
!=======================================================================
enddo
end subroutine user_defined_particle

end module udm_Manager
