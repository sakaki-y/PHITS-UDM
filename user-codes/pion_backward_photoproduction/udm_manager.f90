module udm_Manager
use udm_Parameter

!=======================================================================
! [udm_int]
use udm_int_pion_backward_photoproduction,   caller_udm_int_PBP  => caller
!=======================================================================

implicit none
contains




subroutine user_defined_interaction(action,index)
integer action,index
!=======================================================================
call caller_udm_int_PBP(action,index)
!=======================================================================
end subroutine user_defined_interaction




subroutine user_defined_particle(action)
integer action,index
do index=1,udm_part_nMax
enddo
end subroutine user_defined_particle

end module udm_Manager
