program compare
use dislin
implicit none

integer :: i
integer, parameter :: res = 50
real, dimension(res) :: J, CL, CD, mul, mud
real :: garbage
character(len = 2*6) :: leg

open(unit = 10, file = 'comparison.txt')

do i = 1, res
	read(10,*) J(i), CL(i), CD(i)
	mul(i) = .8665
	mud(i) = .0286
end do

! call metafl('xwin')
call metafl('pdf')
call window(2000,10,1333,1000)
call page(4000,3000)
!call scrmod('reverse')
call disini
call psfont('Times-Roman')
call texmod('on')

call axslen(1000,1000)
call axspos(500,1250)
call name('Coefficient of Lift','y')
call name('Advance Ratio, J','x')
call axsscl('log','x')
call labels('log','x')
! call labtyp('vert','x')
call labdig(1,'y')
call graf(1., 3., 1., 1., 0.7, 1.5, 0.7, .1)
call color('half')
! call grid(0,1)
call color('fore')
call curve(J, CL, res)
call dash
call curve(J, mul, res)
call solid
call legini(leg,2,6)
call legtit('')
call legbgd(0)
!call legpos(2200,1600)
call legpat(5,1,-1,-1,-1,1)
call leglin(leg,'MachUp',1)
call legpat(0,1,-1,-1,-1,2)
call leglin(leg,"PLL",2)
call legend(leg,7)
call endgrf

call axslen(1000,1000)
call axspos(2000,1250)
call name('Coefficient of Drag','y')
call name('Advance Ratio, J','x')
call axsscl('log','x')
call labels('log','x')
! call labtyp('vert','x')
call labdig(3,'y')
call graf(1., 3., 1., 1., -.175, .05, -.175, .025)
call color('half')
! call grid(0,1)
call color('fore')
call curve(J, CD, res)
call dash
call curve(J, mud, res)
call solid
call legini(leg,2,6)
call legtit('')
call legbgd(0)
!call legpos(2200,1600)
call legpat(5,1,-1,-1,-1,1)
call leglin(leg,'MachUp',1)
call legpat(0,1,-1,-1,-1,2)
call leglin(leg,"PLL",2)
call legend(leg,6)
call endgrf


call disfin


end program compare
