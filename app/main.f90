program main
  use kinds_mod
  use rndgen_mod
  implicit none
  type(rndgen) :: generator
  integer :: seed = 294727492
  
  call generator%init(seed)
  print*, generator%rnd()
end program main
