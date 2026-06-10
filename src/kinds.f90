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
        !> arr express, in real double precision, any type of numerical data
        real(dp) , allocatable :: arr( : )
        integer(i4) :: len      !> Integer type data length for loops
        real(dp) :: real_len    !> Real type data length for calculations
        !> Flag to diferentiate between discrete and continuous methods for MLE
        logical :: data_is_discrete
        !> Flag to save in the memory if data is already sorted
        logical :: sorted_data  
    contains
        !> Procedure for receiving data
        procedure :: receive_data        
        !> Procedure for sorting data
        procedure :: sort_data => sorting_random_data
    end type

contains

subroutine receive_data( this , r_data )
    !> This subroutine can recive any numerical data type (except qp), and
    !   converts the data to real double precision, to simply the treatment
    !   However, as a good pratice, we've manteined the "wraps" to the power law distribution, if one wants to evaluete pdf/cdf/ccdf at any numerical type value
    class(random_data) , intent(inout) :: this
    class(*) , intent(in) :: r_data(:) !> Polymorphic variable

    if (allocated(this%arr)) deallocate(this%arr)
    allocate( this%arr(size(r_data)) )
    select type( r_data )
        type is ( real(sp) )
            this%data_is_discrete = .FALSE.
            this%arr = real(r_data, dp)
        type is ( real(dp) )
            this%data_is_discrete = .FALSE.
            this%arr = r_data
        type is ( integer(i1) )
            this%data_is_discrete = .TRUE.
            this%arr = real(r_data, dp)
        type is ( integer(i2) )
            this%data_is_discrete = .TRUE.
            this%arr = real(r_data, dp)
        type is ( integer(i4) )
            this%data_is_discrete = .TRUE.
            this%arr = real(r_data, dp)
        type is ( integer(i8) )
            this%data_is_discrete = .TRUE.
            this%arr = real(r_data, dp)
        class default
            !> If generic_array is a string or an unsupported type, print an error message
            print *, "Error: data type not supported."
        end select

        this%sorted_data = .FALSE.
end subroutine

subroutine sorting_random_data( this , pre_ordering , reverse )
    !> This subroutines plays the role of an interface to ordering the data
    class(random_data) , intent(inout) :: this
    logical, intent(in), optional :: pre_ordering, reverse
    logical :: is_arr_pre_ordering, ordering_in_reverse

    !> Set default values for optional arguments
    is_arr_pre_ordering = .FALSE.
    if (present(pre_ordering)) is_arr_pre_ordering = pre_ordering
    ordering_in_reverse = .FALSE.
    if (present(reverse)) ordering_in_reverse = reverse

    !> Check if the array is allocated
    if (.not. allocated(this%arr)) return

    if (is_arr_pre_ordering) then
        call ord_sort( this%arr , reverse=ordering_in_reverse)
    else
        call radix_sort( this%arr , reverse=ordering_in_reverse)
    endif
    
    this%sorted_data = .TRUE.
    this%len = size( this%arr )
    this%real_len = real( size( this%arr ) , dp )
end subroutine

end module kinds_mod    