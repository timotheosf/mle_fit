module mle_kinds_mod
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
        !> Procedures for receiving data
        procedure :: receive_data     
        procedure :: from_file => receive_data_from_file   
        !> Procedure for sorting data
        procedure :: sort_data => sorting_random_data
    end type

    type :: clock_time !> A benchmark time class
        integer(i4) :: time_at_start
        integer(i4) :: time_at_stop
        real(dp) :: elapsed
        real(dp) :: time_rate
    
    contains
        procedure :: start => start_time_count
        procedure ::  stop =>  stop_time_count
        procedure :: total => elapsed_time
    endtype

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
            error stop "Error: data type not supported."
        end select

        this%sorted_data = .FALSE.
end subroutine

subroutine receive_data_from_file( this , file_unit , skip_title )
    class(random_data) , intent(inout) :: this
    integer(i4) , intent(in) :: file_unit
    logical , intent(in) , optional :: skip_title
    character(len=50) :: string
    integer(i4) :: io_status , file_len , i
    logical :: x_is_integer , last_x_is_integer
    
    rewind( file_unit ) !> Rewind file passed

    if (present(skip_title)) then
        if (skip_title) read(file_unit,*)
    endif

    file_len = 0
    read_file: do
        read(file_unit,*, iostat=io_status) string
        if ( io_status/=0 ) exit read_file
        if (index(string, '.')==0) then
            x_is_integer = .TRUE.
        else
            x_is_integer = .FALSE.
        endif
        if (file_len/=0) then
            if ( x_is_integer .neqv. last_x_is_integer ) error stop "Error: data type must be unique"
        endif
        last_x_is_integer = x_is_integer
        file_len = file_len + 1
    enddo read_file
    if (file_len == 0) error stop "Error: file doens't have any data"
    this%data_is_discrete = x_is_integer
    rewind( file_unit ) !> Rweind file again
    if (present(skip_title)) then
        if (skip_title) read(file_unit,*)
    endif
    allocate( this%arr(file_len) )
    do i = 1 , file_len
        read(file_unit,*) this%arr(i)
    enddo

    call this%sort_data()
end subroutine
    

subroutine sorting_random_data( this , pre_ordering , reverse , work_buffer )
    !> This subroutines plays the role of an interface to ordering the data
    class(random_data) , intent(inout) :: this
    logical, intent(in), optional :: pre_ordering, reverse
    logical :: is_arr_pre_ordering, ordering_in_reverse
    real(dp), intent(inout), optional :: work_buffer(:)

    !> Set default values for optional arguments
    is_arr_pre_ordering = .FALSE.
    if (present(pre_ordering)) is_arr_pre_ordering = pre_ordering
    ordering_in_reverse = .FALSE.
    if (present(reverse)) ordering_in_reverse = reverse
    if (present(work_buffer)) then
        if (size(work_buffer)/=size(this%arr)) error stop "Sorting random data: Passed work buffer has to be the same size data array"
    endif

    !> Check if the array is allocated
    if (.not. allocated(this%arr)) return

        if (is_arr_pre_ordering) then
        if (present(work_buffer)) then
            call ord_sort( this%arr , work=work_buffer, reverse=ordering_in_reverse)
        else
            call ord_sort( this%arr , reverse=ordering_in_reverse)
        endif
    else
        if (present(work_buffer)) then
            call radix_sort( this%arr , work=work_buffer, reverse=ordering_in_reverse)
        else
            call radix_sort( this%arr , reverse=ordering_in_reverse)
        endif
    endif
    
    this%sorted_data = .TRUE.
    this%len = size( this%arr )
    this%real_len = real( size( this%arr ) , dp )
end subroutine

subroutine start_time_count( this )
    class(clock_time) , intent(inout) :: this
    integer(i4) :: rate
    call system_clock(count_rate=rate)
    this%time_rate = real(rate,dp)
    call system_clock(count=this%time_at_start)
end subroutine

subroutine stop_time_count( this )
    class(clock_time) , intent(inout) :: this    
    call system_clock(count=this%time_at_stop)
    this%elapsed = real( this%time_at_stop - this%time_at_start , dp ) / this%time_rate
end subroutine

function elapsed_time( this ) result( elapsed_seconds )
    class(clock_time) , intent(in) :: this
    real(dp) :: elapsed_seconds
    elapsed_seconds = real( this%time_at_stop - this%time_at_start , dp ) / this%time_rate
end function

end module mle_kinds_mod    