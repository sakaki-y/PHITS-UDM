module udm_Manager
use udm_Parameter

!=======================================================================
! [udm_int]
use udm_int_sample,   caller_udm_int_sample  => caller
! [udm_part]
use udm_part_sample,   caller_udm_part_sample  => caller
!=======================================================================

implicit none
contains




subroutine user_defined_interaction(action,index)
integer action,index
!=======================================================================
call caller_udm_int_sample(action,index)
!=======================================================================
end subroutine user_defined_interaction




subroutine user_defined_particle(action)
integer action,index
do index=1,udm_part_nMax
!=======================================================================
call caller_udm_part_sample(action,index)
!=======================================================================
enddo
end subroutine user_defined_particle

end module udm_Manager
