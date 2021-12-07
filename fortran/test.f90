program test
use points_mod
use helix_func
implicit none

real(wp) :: a, b, n, m, p
type(point) :: cp
type(vector) :: v, garbage

a = 1._wp
b = .5_wp
n = 100._wp
m = 50._wp
p = 1.5_wp
cp%x = 0.9_wp
cp%y = 0._wp
cp%z = 0._wp

call velocity_dim(v, garbage, 1._wp, a, b, n, cp, m, p, 0._wp, 'r')

write(*,*) v

end program test

