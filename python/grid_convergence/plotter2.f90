program plotter
use dislin
implicit none

integer :: i, j

real, dimension(1) :: adv_ratio
real, dimension(1,7) :: CT, CP, CL
real, dimension(7) :: nodes
real :: ymax, ymin, y_dist
character(len = 3) :: str
character(len = 7 * 3) :: leg

open(unit = 10, file = 'C_vs_J_take_2.txt')

do i = 1, 7
	do j = 1, 1
		read(10,*) adv_ratio(j), CT(j,i), CP(j,i), CL(j,i), nodes(i)
	end do
end do

do j = 1, 1
	
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
	call name('Coefficient of Power','y')
	call name('Nodes per Blade','x')
	call axsscl('log','x')
	call labels('log','x')
	
	ymax = maxval(CT(j,:))
	ymin = minval(CT(j,:))
	y_dist = ymax - ymin
	ymax = ymax + .1 * y_dist
	ymin = ymin - .1 * y_dist
	
	i = int(ceiling(abs(log10(.1*y_dist))))
	call labdig(i,'y')
	
	call graf(.95,2.9,1.,1.,ymin,ymax,ymin,.1*y_dist)
	call incmrk(-1)
	call marker(2)
	call curve(nodes, CT(j,:), 7)
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
call name('Coefficient of Power','y')
call name('Advance Ratio','x')
call labtyp('vert','x')
! call axsscl('log','x')
! call labels('log','x')

ymax = maxval(CT(:,:))
ymin = minval(CT(:,:))
y_dist = ymax - ymin
ymax = ymax + .1 * y_dist
ymin = ymin - .1 * y_dist

j = int(ceiling(abs(log10(.1*y_dist))))
call labdig(j,'y')

call graf(-.1,1.1,-.1,.1,ymin,ymax,ymin,.1*y_dist)
! call chncrv('color')
call setvlt('grey')

call legini(leg,7,3)
call legtit('Nodes per Blade')
! call legbgd(0)
j = 10

10 format(I3)

call incmrk(-1)


do i = 1, 7
	
	call marker(15)
	call hsymbl(71-i*10)
	call setclr(211-i*30)
	call curve(adv_ratio, CT(:,i), 1)
! 	call sendbf
! 	read(*,*)
! 	str = char(j)
	write(str,10) j
	call color('fore')
! 	call legval(real(71-i*10)/10.,'symbol')
	call leglin(leg,str,i)
	
	j = j * 2
end do

call legend(leg,7)

call disfin

end program plotter
