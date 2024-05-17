module udm_Manager
use udm_Parameter

!=======================================================================
! [udm_int]
use udm_int_new_boson_radiation,  caller_udm_int_new_boson_radiation  => caller
use udm_int_sample,  caller_udm_int_sample  => caller
use udm_int_ALPs,  caller_udm_int_ALPs  => caller

! [udm_part]
use udm_part_new_boson, caller_udm_part_new_boson => caller
use udm_part_MuonMinus_NoDecay, caller_udm_part_MuonMinus_NoDecay => caller
use udm_part_MuonPlus_NoDecay , caller_udm_part_MuonPlus_NoDecay  => caller
use udm_part_sample , caller_udm_part_sample  => caller
use udm_part_ALPs, caller_udm_part_ALPs => caller
!=======================================================================

implicit none
contains




subroutine user_defined_interaction(action,index)
integer action,index
!=======================================================================
call caller_udm_int_new_boson_radiation(action,index)
call caller_udm_int_sample(action,index)
call caller_udm_int_ALPs(action,index)
!=======================================================================
end subroutine user_defined_interaction




subroutine user_defined_particle(action)
integer action,index
do index=1,udm_part_nMax
!=======================================================================
call caller_udm_part_new_boson(action,index)
call caller_udm_part_MuonMinus_NoDecay(action,index)
call caller_udm_part_MuonPlus_NoDecay (action,index)
call caller_udm_part_sample(action,index)
call caller_udm_part_ALPs(action,index)
!=======================================================================
enddo
end subroutine user_defined_particle

end module udm_Manager
