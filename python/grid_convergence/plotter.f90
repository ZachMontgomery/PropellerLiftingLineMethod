program plotter
use dislin
implicit none

integer :: i, j

integer, parameter :: res_grid = 5, res_J = 11

real, dimension(res_J) :: adv_ratio
real, dimension(res_J,res_grid) :: CT, CP, CL, C
real, dimension(res_grid) :: nodes
real :: ymax, ymin, y_dist
character(len = 3) :: str
character(len = res_grid * 3) :: leg
character(len = 21) :: ylabel

open(unit = 10, file = 'C_vs_J_bet.txt')

do i = 1, res_grid
	do j = 1, res_J
		read(10,*) adv_ratio(j), CT(j,i), CP(j,i), CL(j,i), nodes(i)
	end do
end do

! ylabel = 'Coefficient of Thrust'
! C = CT
ylabel = 'Coefficient of Power'
C = CP

do j = 1, res_J
	
	call metafl('xwin')
! 	call metafl('pdf')
	call window(2000,10,1333,1000)
	call page(4000,3000)
	
	call disini
	call psfont('Times-Roman')
	call texmod('on')
	
	call messag('$J = $',900,100)
	call number(adv_ratio(j),1,999,999)
	
	call axslen(1000,1000)
	call axspos(500,1250)
	call name(ylabel,'y')
	call name('Nodes per Blade','x')
	call axsscl('log','x')
	call labels('log','x')
	
	ymax = maxval(C(j,:))
	ymin = minval(C(j,:))
	y_dist = ymax - ymin
	ymax = ymax + .1 * y_dist
	ymin = ymin - .1 * y_dist
	
	i = int(ceiling(abs(log10(.1*y_dist))))
	call labdig(i,'y')
	
	call graf(.95,2.3,1.,1.,ymin,ymax,ymin,.1*y_dist)
	call incmrk(-1)
	call marker(2)
	call curve(nodes, C(j,:), res_grid)
	call disfin
	
end do

call metafl('xwin')
! call metafl('pdf')
call scrmod('reverse')
call window(2000,10,1333,1000)
call page(4000,3000)

call disini
call psfont('Times-Roman')
call texmod('on')

! call messag('$J = $',900,100)
! call number(adv_ratio(j),1,999,999)

call axslen(1000,1000)
call axspos(500,1250)
call name(ylabel,'y')
call name('Advance Ratio','x')
call labtyp('vert','x')
! call axsscl('log','x')
! call labels('log','x')

ymax = maxval(C(:,:))
ymin = minval(C(:,:))
ymax = .045
y_dist = ymax - ymin
ymax = ymax + .1 * y_dist
ymin = ymin - .1 * y_dist

j = int(ceiling(abs(log10(.1*y_dist))))
call labdig(j,'y')

call graf(-.1,1.1,-.1,.1,ymin,ymax,ymin,.1*y_dist)
! call chncrv('color')
call setvlt('grey')

call legini(leg,res_grid,3)
call legtit('Nodes per Blade')
! call legbgd(0)
j = 10

10 format(I3)

call incmrk(-1)


do i = 1, res_grid
	
	call marker(15)
	call hsymbl(71-i*10)
	call setclr(211-i*30)
	call curve(adv_ratio, C(:,i), res_J)
! 	call sendbf
! 	read(*,*)
! 	str = char(j)
	write(str,10) j
	call color('fore')
! 	call legval(real(71-i*10)/10.,'symbol')
	call leglin(leg,str,i)
	
	j = j * 2
end do

call legend(leg,9)

! call rlsymb(0,.5,.05)

call disfin

end program plotter
