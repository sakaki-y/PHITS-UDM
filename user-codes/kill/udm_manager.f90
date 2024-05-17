module udm_Manager
use udm_Parameter

!=======================================================================
use udm_part_kill,   caller_udm_part_kill  => caller
!=======================================================================

implicit none
contains




subroutine user_defined_interaction(action,index)
integer action,index
end subroutine user_defined_interaction




subroutine user_defined_particle(action)
integer action,index
do index=1,udm_part_nMax
!=======================================================================
call caller_udm_part_kill(action,index)
!=======================================================================
enddo
end subroutine user_defined_particle

end module udm_Manager
