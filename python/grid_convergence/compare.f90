program compare
use dislin
implicit none

integer :: i
integer, parameter :: res = 100
real, dimension(res) :: zeta, t_pll, d_pll, t_bet, d_bet, tw_pll, tw_bet
real :: garbage, tmax, dmax, temp
character(len = 2*3) :: leg

open(unit = 10, file = 'pll_blade_loading.txt')
open(unit = 20, file = 'bet_blade_loading.txt')

do i = 1, res
	read(10,*) zeta(i), d_pll(i), t_pll(i), tw_pll(i)
	read(20,*) garbage, d_bet(i), t_bet(i), tw_bet(i)
end do

temp = maxval(-t_pll)
tmax = maxval(-t_bet)
if (temp > tmax) tmax = temp
tmax = tmax / .9

temp = maxval(d_pll)
dmax = maxval(d_bet)
if (temp > dmax) dmax = temp
dmax = dmax / .9

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
call name('Thrust (lbf)','y')
call name('Blade Span (ft)','x')
call labtyp('vert','x')
call labdig(3,'y')
call graf(0., 1., 0., .1, 0., tmax, 0., .1 * tmax)
call color('half')
call grid(1,1)
call color('fore')
call curve(zeta, -t_pll, res)
call dash
call curve(zeta, -t_bet, res)
call solid
call legini(leg,2,3)
call legtit('')
call legbgd(0)
!call legpos(2200,1600)
call legpat(5,1,-1,-1,-1,1)
call leglin(leg,'BET',1)
call legpat(0,1,-1,-1,-1,2)
call leglin(leg,"PLL",2)
call legend(leg,8)
call endgrf

call axslen(1000,1000)
call axspos(2000,1250)
call name('Drag (lbf)','y')
call name('Blade Span (ft)','x')
call labtyp('vert','x')
call labdig(3,'y')
call graf(0., 1., 0., .1, 0., dmax, 0., .1 * dmax)
call color('half')
call grid(1,1)
call color('fore')
call curve(zeta, d_pll, res)
call dash
call curve(zeta, d_bet, res)
call solid
call legini(leg,2,3)
call legtit('')
call legbgd(0)
!call legpos(2200,1600)
call legpat(5,1,-1,-1,-1,1)
call leglin(leg,'BET',1)
call legpat(0,1,-1,-1,-1,2)
call leglin(leg,"PLL",2)
call legend(leg,7)
call endgrf


call disfin

end program compare
