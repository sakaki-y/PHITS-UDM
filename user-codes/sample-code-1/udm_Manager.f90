module udm_Manager
use udm_Parameter

!=======================================================================
! [udm_int]
use udm_int_sample_1,   caller_udm_int_sample_1  => caller
use udm_int_sample_2,   caller_udm_int_sample_2  => caller
! [udm_part]
use udm_part_sample_1,  caller_udm_part_sample_1 => caller
use udm_part_sample_2,  caller_udm_part_sample_2 => caller
!=======================================================================

implicit none
contains




subroutine user_defined_interaction(action,index)
integer action,index
!=======================================================================
call caller_udm_int_sample_1(action,index)
call caller_udm_int_sample_2(action,index)
!=======================================================================
end subroutine user_defined_interaction




subroutine user_defined_particle(action)
integer action,index
do index=1,udm_part_nMax
!=======================================================================
call caller_udm_part_sample_1(action,index)
call caller_udm_part_sample_2(action,index)
!=======================================================================
enddo
end subroutine user_defined_particle

end module udm_Manager
