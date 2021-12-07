program test
implicit none

character(100) :: filename, line
integer :: i, j

call get_command_argument(1,filename)

open(unit = 10, file = filename, action = 'read')

read(10,*) i

do j = 1, i
	read(10,*) line
	write(*,*) line
end do

call get_command_argument(0,filename)

write(*,*) filename

end program test
