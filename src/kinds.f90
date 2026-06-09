module kinds_mod
    use, intrinsic :: iso_fortran_env, only: &
        sp => real32, & ! single precision, range: -3.40282347E+38 to 3.40282347E+38, smallest positive: 1.17549435E-38
        dp => real64, & ! double precision, range: -1.7976931348623157E+308 to 1.7976931348623157E+308, smallest positive: 2.2250738585072014E-308
        qp => real128,& ! quad precision
        i1 => int8, &   ! 1 byte integer, range: -128 to 127
        i2 => int16, &  ! 2 byte integer, range: -32768 to 32767
        i4 => int32, &  ! 4 byte integer, range: -2147483648 to 2147483647
        i8 => int64     ! 8 byte integer, range: -9223372036854775808 to 9223372036854775807
    use :: stdlib_sorting, only: &
        radix_sort, &
        ord_sort
    implicit none
    
    type :: random_data
        class(*) , allocatable :: arr( : ) !> Polymorphic variable
        logical :: sorted_data
        integer(i4) :: length
        real(dp) :: real_length
        

    contains
        
        !> Procedure for receiving data
        procedure :: receive_data
        
        !> Procedure for sorting data
        procedure :: sort_data => sorting_random_data
    end type

contains

subroutine receive_data( this , r_data )
    class(random_data) , intent(inout) :: this
    class(*) , intent(in) :: r_data(:)

    if (allocated(this%arr)) deallocate(this%arr)
    
    select type( r_data )
        type is ( real(sp) )
            allocate( this%arr , source=r_data )
        type is ( real(dp) )
            allocate( this%arr , source=r_data )
        type is ( integer(i1) )
            allocate( this%arr , source=r_data )
        type is ( integer(i2) )
            allocate( this%arr , source=r_data )
        type is ( integer(i4) )
            allocate( this%arr , source=r_data )
        type is ( integer(i8) )
            allocate( this%arr , source=r_data )
        class default
            !> If generic_array is a string or an unsupported type, print an error message
            print *, "Error: data type not supported."
        end select
end subroutine

subroutine sorting_random_data( this , pre_ordering , reverse )
    class(random_data) , intent(inout) :: this
    logical, intent(in), optional :: pre_ordering, reverse
    logical :: is_arr_pre_ordering, ordering_in_reverse

    !> Set default values for optional arguments
    is_arr_pre_ordering = .false.
    if (present(pre_ordering)) is_arr_pre_ordering = pre_ordering
    ordering_in_reverse = .false.
    if (present(reverse)) ordering_in_reverse = reverse

    !> Check if the array is allocated
    if (.not. allocated(this%arr)) return

    select type( generic_array => this%arr )
        type is ( real(sp) )
            if (is_arr_pre_ordering) then
                call ord_sort( generic_array , reverse=ordering_in_reverse)
            else
                call radix_sort( generic_array , reverse=ordering_in_reverse)
           endif
        type is ( real(dp) )
            if (is_arr_pre_ordering) then
                call ord_sort( generic_array , reverse=ordering_in_reverse)
            else
                call radix_sort( generic_array , reverse=ordering_in_reverse)
           endif
        type is ( integer(i1) )
            if (is_arr_pre_ordering) then
                call ord_sort( generic_array , reverse=ordering_in_reverse)
            else
                call radix_sort( generic_array , reverse=ordering_in_reverse)
           endif
        type is ( integer(i2) )
            if (is_arr_pre_ordering) then
                call ord_sort( generic_array , reverse=ordering_in_reverse)
            else
                call radix_sort( generic_array , reverse=ordering_in_reverse)
           endif
        type is ( integer(i4) )
            if (is_arr_pre_ordering) then
                call ord_sort( generic_array , reverse=ordering_in_reverse)
            else
                call radix_sort( generic_array , reverse=ordering_in_reverse)
           endif
        type is ( integer(i8) )
            if (is_arr_pre_ordering) then
                call ord_sort( generic_array , reverse=ordering_in_reverse)
            else
                call radix_sort( generic_array , reverse=ordering_in_reverse)
           endif
        class default
            !> If generic_array is a string or an unsupported type, print an error message
            print *, "Error: data type not supported."
        end select
        this%sorted_data = .true.
        this%length = size( this%arr )
        this%real_length = real( size( this%arr ) , dp )
end subroutine

end module kinds_mod    