program compare
use dislin
implicit none

integer :: i
integer, parameter :: res = 101
real, dimension(res) :: J, ct_pll, cp_pll, cl_pll, ct_bet, cp_bet, cl_bet
real :: garbage
character(len = 2*3) :: leg

open(unit = 10, file = 'pll_ct_vs_J.txt')
open(unit = 20, file = 'bet_Ct_vs_J.txt')

do i = 1, res
	read(10,*) J(i), ct_pll(i), cp_pll(i), cl_pll(i)
	read(20,*) garbage, ct_bet(i), cp_bet(i), cl_bet(i)
end do

!call metafl('xwin')
call metafl('pdf')
call window(2000,10,1333,1000)
call page(4000,3000)
!call scrmod('reverse')
call disini
call psfont('Times-Roman')
call texmod('on')

call axslen(1000,1000)
call axspos(500,1250)
call name('Coefficient of Thrust','y')
call name('Advance Ratio, J','x')
call labtyp('vert','x')
call labdig(2,'y')
call graf(0., 1., 0., .1, 0., .101, 0., .01)
call color('half')
call grid(1,1)
call color('fore')
call curve(J, ct_pll, res)
call dash
call curve(J, ct_bet, res)
call solid
call legini(leg,2,3)
call legtit('Method')
call legbgd(0)
!call legpos(2200,1600)
call legpat(5,1,-1,-1,-1,1)
call leglin(leg,'BET',1)
call legpat(0,1,-1,-1,-1,2)
call leglin(leg,"PLL",2)
call legend(leg,5)
call endgrf

call axslen(1000,1000)
call axspos(2000,1250)
call name('Coefficient of Power','y')
call name('Advance Ratio, J','x')
call labtyp('vert','x')
call labdig(3,'y')
call graf(0., 1., 0., .1, 0., .05, 0., .005)
call color('half')
call grid(1,1)
call color('fore')
call curve(J, cp_pll, res)
call dash
call curve(J, cp_bet, res)
call solid
call legini(leg,2,15)
call legtit('Method')
call legbgd(0)
!call legpos(2200,1600)
call legpat(5,1,-1,-1,-1,1)
call leglin(leg,'BET',1)
call legpat(0,1,-1,-1,-1,2)
call leglin(leg,"PLL",2)
call legend(leg,5)
call endgrf

call axslen(1000,1000)
call axspos(500,2700)
call name('Coefficient of Torque','y')
call name('Advance Ratio, J','x')
call labtyp('vert','x')
call labdig(4,'y')
call graf(0., 1., 0., .1, 0., .011, 0., .0011)
call color('half')
call grid(1,1)
call color('fore')
call curve(J, cl_pll, res)
call dash
call curve(J, cl_bet, res)
call solid
call legini(leg,2,15)
call legtit('Method')
call legbgd(0)
!call legpos(2200,1600)
call legpat(5,1,-1,-1,-1,1)
call leglin(leg,'BET',1)
call legpat(0,1,-1,-1,-1,2)
call leglin(leg,"PLL",2)
call legend(leg,5)
call endgrf


call disfin


end program compare
