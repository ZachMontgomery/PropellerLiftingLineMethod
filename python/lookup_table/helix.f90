
!***************************************************************************
! dimensional
!***************************************************************************
subroutine calc_nu(nu, a1, a2, b, rl, control, n, m, p, phi, int_type, ij)
	implicit none
	integer, parameter :: wp = selected_real_kind(p=8)
	real(wp), intent(in) :: a1, a2, b, n, m, p, phi, rl
	real(wp), dimension(0:2), intent(in) :: control
	real(wp), dimension(0:2) :: temp
	real(wp), dimension(0:2), intent(inout) :: nu
!f2py intent(in,out) :: nu
	character, intent(in) :: int_type
	logical, intent(in) :: ij
	real(wp), dimension(0:2) :: v, w, r1, r2
	real(wp), dimension(0:2) :: j1, j2
	integer :: i, steps
	real(wp), allocatable, dimension(:) :: zeta
	real(wp), parameter :: pi = 3.14159265358979323846_wp
	
	call integrate(v, a1, b, rl, control, n, m, p, phi, int_type)
	!*******************************************************************************************
	call integrate(w, a2, b, rl, control, n, m, p, phi, int_type)
	!************************************************************************************
	if (.not. ij) then
		j1(0) = a1 * cos(phi)
		j1(1) = a1 * sin(phi)
		j1(2) = 0._wp
		
		j2(0) = a2 * cos(phi)
		j2(1) = a2 * sin(phi)
		j2(2) = 0._wp
		
		r1 = control - j1
		r2 = control - j2
		
		call straight_segment(temp, r1,r2)
		temp = -rl * temp
	end if
	! set the output variable
	nu = w - v + temp
	
	
end subroutine calc_nu
!***************************************************************************
! my functions
!***************************************************************************
subroutine dvx(dx,zeta,denom,cx1,cx2,cx3,cx4,s,c)
	implicit none
	integer, parameter :: wp = selected_real_kind(p=8)
	real(wp), intent(in) :: zeta,denom,cx1,cx2,cx3,cx4,s,c
	real(wp), intent(out) :: dx
	dx = (cx1 + cx2 * s + (cx3 + cx4 * zeta) * c)/denom
end subroutine dvx
!***************************************************************************
subroutine dvy(dy,zeta,denom,cy1,cy2,cy3,cy4,s,c)
	implicit none
	integer, parameter :: wp = selected_real_kind(p=8)
	real(wp), intent(in) :: zeta,denom,cy1,cy2,cy3,cy4,s,c
	real(wp), intent(out) :: dy
	dy = (cy1 + cy2 * c + (cy3 + cy4 * zeta) * s)/denom
end subroutine dvy
!***************************************************************************
subroutine dvz(dz,denom,cz1,cz2,cz3,s,c)
	implicit none
	integer, parameter :: wp = selected_real_kind(p=8)
	real(wp), intent(in) :: denom,cz1,cz2,cz3,s,c
	real(wp), intent(out) :: dz
	dz = (cz1 + cz2 * c + cz3 * s)/denom
end subroutine dvz
!***************************************************************************
subroutine integrate(v, a, b, rl, control, n, m, p, phi, int_type)
	implicit none
	integer, parameter :: wp = selected_real_kind(p=8)
	real(wp), dimension(0:2), intent(inout) :: v
!f2py intent(in,out) :: v
	real(wp), intent(in) :: a, b, rl, n, m, p, phi
	real(wp), dimension(0:2), intent(in) :: control
	character, intent(in) :: int_type
	
	real(wp) :: c0, c1, c2, cx1, cx2, cx3, cx4, cy1, cy2, cy3, cy4, cz1, cz2, cz3, d1, d2
	real(wp) :: dzeta, theta, denom, c, s, zeta_temp
	integer :: steps, i
	real(wp), allocatable, dimension(:) :: zeta
	real(wp), dimension(0:2) :: k1, k2, k3, k4, temp, dv_o, dv_n
	real(wp), parameter :: pi = 3.14159265358979323846_wp
	
	steps = nint(n*m)
	allocate(zeta(0:steps))
	! calculate constants
	c0 = 2._wp * pi * n
	c1 = n / 4._wp / pi
	c2 = a * n / 2._wp
	! x function contstants
	cx1 = -b * control(1)
	cx2 = -rl * a * b
	cx3 = -rl * 2._wp * pi * a * control(2)
	cx4 = rl * 2._wp * pi * a * b * n
	! y function contstants
	cy1 = b * control(0)
	cy2 = -a * b
	cy3 = 2._wp * pi * a * control(2)
	cy4 = -2._wp * pi * a * b * n
	! z function constants
	cz1 = -rl * a
	cz2 = rl * control(0)
	cz3 = -control(1)
	d1 = b * n
	d2 = rl * a
	! calculate size of integration step and number of steps
	dzeta = 1._wp / n / m
	do i = 0, steps
		zeta(i) = (real(i,wp)*dzeta)**p
	end do
	! initialize vectors
	v = 0._wp
	k1 = v
	k2 = v
	k3 = v
	k4 = v
	temp = v
	! initialize values for loop
	theta = phi
	s = sin(theta)
	c = cos(theta)
	denom = ((control(0)-a*c)**2._wp+(control(1)+d2*s)**2._wp+control(2)**2._wp)**(1.5_wp)
	call dvx(dv_o(0),0._wp,denom,cx1,cx2,cx3,cx4,s,c)
	call dvy(dv_o(1), 0._wp,denom,cy1,cy2,cy3,cy4,s,c)
	call dvz(dv_o(2), denom, cz1, cz2, cz3, s, c)
	! check for trapezoidal rule or runge-kutta 4
	if (int_type == 't') then
		! loop for all segments
		do i = 1, steps
			! calculate common items between the x, y, and z velocities
			dzeta = zeta(i) - zeta(i-1)
			theta = c0 * zeta(i) + phi
			s = sin(theta)
			c = cos(theta)
			denom = ((control(0)-a*c)**2._wp+(control(1)+d2*s)**2._wp+(control(2)-d1*zeta(i))**2._wp)**(1.5_wp)
			call dvx(dv_n(0), zeta(i),denom,cx1,cx2,cx3,cx4,s,c)
			call dvy(dv_n(1), zeta(i),denom,cy1,cy2,cy3,cy4,s,c)
			call dvz(dv_n(2), denom,cz1,cz2,cz3,s,c)
			! add on the trapezoidal area
			v = v + dzeta / 2._wp * (dv_o + dv_n)
			! save values
			dv_o = dv_n
		end do
		v(0) = -rl * v(0) * c1
		v(1) = -rl * v(1) * c1
		v(2) = -rl * v(2) * c2
	elseif (int_type == 'r') then
		! loop for all segments
		do i = 1, steps
			! calculate dzeta
			dzeta = zeta(i) - zeta(i-1)
			! calculate k1
			k1 = dv_o
			! calculate k2
			zeta_temp = zeta(i-1) + dzeta / 2._wp
			theta = c0 * zeta_temp + phi
			s = sin(theta)
			c = cos(theta)
			denom = ((control(0)-a*c)**2._wp+(control(1)+d2*s)**2._wp+(control(2)-d1*zeta_temp)**2._wp)**(1.5_wp)
			call dvx(k2(0), zeta_temp,denom,cx1,cx2,cx3,cx4,s,c)
			call dvy(k2(1), zeta_temp,denom,cy1,cy2,cy3,cy4,s,c)
			call dvz(k2(2), denom,cz1,cz2,cz3,s,c)
			! calculate k3
			!k3 = k2
			! calculate k4
			theta = c0 * zeta(i) + phi
			s = sin(theta)
			c = cos(theta)
			denom = ((control(0)-a*c)**2._wp+(control(1)+d2*s)**2._wp+(control(2)-d1*zeta(i))**2._wp)**(1.5_wp)
			call dvx(k4(0), zeta(i),denom,cx1,cx2,cx3,cx4,s,c)
			call dvy(k4(1), zeta(i),denom,cy1,cy2,cy3,cy4,s,c)
			call dvz(k4(2), denom,cz1,cz2,cz3,s,c)
			dv_n = k4
			! add on to the final values
			v = v + dzeta / 6._wp * (k1 + 4._wp * k2 + k4)
			! save values
			dv_o = dv_n
		end do
		v(0) = -rl * v(0) * c1
		v(1) = -rl * v(1) * c1
		v(2) = -rl * v(2) * c2
	else
		write(*,*)
		write(*,*) '**************************************'
		write(*,*) "Error int_type must be 't' or 'r'"
		write(*,*) '**************************************'
		write(*,*)
	end if
	!*******************************************************************************************
	deallocate(zeta)
end subroutine integrate
!***************************************************************************
subroutine integrate_nondim(v, lambda, rl, cyl, n, m, p, phi, int_type)
	implicit none
	integer, parameter :: wp = selected_real_kind(p=8)
	real(wp), dimension(0:2), intent(inout) :: v
!f2py intent(in,out) :: v
	real(wp), intent(in) :: lambda, rl, n, m, p, phi
	real(wp), dimension(0:2), intent(in) :: cyl
	character, intent(in) :: int_type
	
	real(wp) :: c0, c1, c2, cx1, cx2, cx3, cx4, cy1, cy2, cy3, cy4, cz1, cz2, cz3, d1, d2
	real(wp) :: dzeta, theta, denom, c, s, zeta_temp
	integer :: steps, i
	real(wp), allocatable, dimension(:) :: zeta
	real(wp), dimension(0:2) :: k1, k2, k3, k4, temp, dv_o, dv_n, control
	real(wp), parameter :: pi = 3.14159265358979323846_wp
	
	control(0) = cyl(0) * cos(cyl(1))
	control(1) = cyl(0) * sin(cyl(1))
	control(2) = cyl(2)
	
	steps = nint(n*m)
	allocate(zeta(0:steps))
	! calculate constants
	c0 = 2._wp * pi * n
	c1 = n / 4._wp / pi
	c2 = n / 2._wp
	! x function contstants
	cx1 = -lambda * control(1)
	cx2 = -rl * lambda
	cx3 = -rl * 2._wp * pi * control(2)
	cx4 = rl * 2._wp * pi * lambda * n
	! y function contstants
	cy1 = lambda * control(0)
	cy2 = -lambda
	cy3 = 2._wp * pi * control(2)
	cy4 = -2._wp * pi * lambda * n
	! z function constants
	cz1 = -rl
	cz2 = rl * control(0)
	cz3 = -control(1)
	d1 = lambda * n
	d2 = rl
	! calculate size of integration step and number of steps
	dzeta = 1._wp / n / m
	do i = 0, steps
		zeta(i) = (real(i,wp)*dzeta)**p
	end do
	! initialize vectors
	v = 0._wp
	k1 = v
	k2 = v
	k3 = v
	k4 = v
	temp = v
	! initialize values for loop
	theta = phi
	s = sin(theta)
	c = cos(theta)
	denom = ((control(0)-c)**2._wp+(control(1)+d2*s)**2._wp+control(2)**2._wp)**(1.5_wp)
	call dvx(dv_o(0),0._wp,denom,cx1,cx2,cx3,cx4,s,c)
	call dvy(dv_o(1), 0._wp,denom,cy1,cy2,cy3,cy4,s,c)
	call dvz(dv_o(2), denom, cz1, cz2, cz3, s, c)
	! check for trapezoidal rule or runge-kutta 4
	if (int_type == 't') then
		! loop for all segments
		do i = 1, steps
			! calculate common items between the x, y, and z velocities
			dzeta = zeta(i) - zeta(i-1)
			theta = c0 * zeta(i) + phi
			s = sin(theta)
			c = cos(theta)
			denom = ((control(0)-c)**2._wp+(control(1)+d2*s)**2._wp+(control(2)-d1*zeta(i))**2._wp)**(1.5_wp)
			call dvx(dv_n(0), zeta(i),denom,cx1,cx2,cx3,cx4,s,c)
			call dvy(dv_n(1), zeta(i),denom,cy1,cy2,cy3,cy4,s,c)
			call dvz(dv_n(2), denom,cz1,cz2,cz3,s,c)
			! add on the trapezoidal area
			v = v + dzeta / 2._wp * (dv_o + dv_n)
			! save values
			dv_o = dv_n
		end do
		v(0) = -rl * v(0) * c1
		v(1) = -rl * v(1) * c1
		v(2) = -rl * v(2) * c2
	elseif (int_type == 'r') then
		! loop for all segments
		do i = 1, steps
			! calculate dzeta
			dzeta = zeta(i) - zeta(i-1)
			! calculate k1
			k1 = dv_o
			! calculate k2
			zeta_temp = zeta(i-1) + dzeta / 2._wp
			theta = c0 * zeta_temp + phi
			s = sin(theta)
			c = cos(theta)
			denom = ((control(0)-c)**2._wp+(control(1)+d2*s)**2._wp+(control(2)-d1*zeta_temp)**2._wp)**(1.5_wp)
			call dvx(k2(0), zeta_temp,denom,cx1,cx2,cx3,cx4,s,c)
			call dvy(k2(1), zeta_temp,denom,cy1,cy2,cy3,cy4,s,c)
			call dvz(k2(2), denom,cz1,cz2,cz3,s,c)
			! calculate k3
			!k3 = k2
			! calculate k4
			theta = c0 * zeta(i) + phi
			s = sin(theta)
			c = cos(theta)
			denom = ((control(0)-c)**2._wp+(control(1)+d2*s)**2._wp+(control(2)-d1*zeta(i))**2._wp)**(1.5_wp)
			call dvx(k4(0), zeta(i),denom,cx1,cx2,cx3,cx4,s,c)
			call dvy(k4(1), zeta(i),denom,cy1,cy2,cy3,cy4,s,c)
			call dvz(k4(2), denom,cz1,cz2,cz3,s,c)
			dv_n = k4
			! add on to the final values
			v = v + dzeta / 6._wp * (k1 + 4._wp * k2 + k4)
			! save values
			dv_o = dv_n
		end do
		v(0) = -rl * v(0) * c1
		v(1) = -rl * v(1) * c1
		v(2) = -rl * v(2) * c2
	else
		write(*,*)
		write(*,*) '**************************************'
		write(*,*) "Error int_type must be 't' or 'r'"
		write(*,*) '**************************************'
		write(*,*)
	end if
	!*******************************************************************************************
	deallocate(zeta)
end subroutine integrate_nondim
!***************************************************************************
subroutine interpolate(x, y1, y2, y, x1, x2)
	implicit none
	integer, parameter :: wp = selected_real_kind(p=8)
	real(wp), intent(in) :: y1, y2, y, x1, x2
	real(wp), intent(out) :: x
	x = (y-y1)*(x2-x1)/(y2-y1)+x1
end subroutine interpolate










!***************************************************************************
! HAM functions
!***************************************************************************
subroutine straight_segment(r,r1,r2)
	implicit none
	integer, parameter :: wp = selected_real_kind(p=8)
	real(wp), dimension(0:2), intent(in) :: r1, r2
	real(wp), dimension(0:2), intent(inout) :: r
!f2py intent(in,out) :: r
	real(wp) :: r1mag, r2mag, r1dotr2
	real(wp), dimension(0:2) :: r1crossr2
	real(wp), parameter :: pi = 3.14159265358979323846_wp
	
	r1crossr2(0) = r1(1) * r2(2) - r1(2) * r2(1)
	r1crossr2(1) = r1(2) * r2(0) - r1(0) * r2(2)
	r1crossr2(2) = r1(0) * r2(1) - r1(1) * r2(0)
	
	r1dotr2 = r1(0) * r2(0) + r1(1) * r2(1) + r1(2) * r2(2)
	
	r1mag = sqrt(r1(0)**2. + r1(1)**2. + r1(2)**2.)
	r2mag = sqrt(r2(0)**2. + r2(1)**2. + r2(2)**2.)
	
	r = .25_wp / pi * (r1mag + r2mag) / r1mag / r2mag / (r1mag * r2mag + r1dotr2) * r1crossr2
	
end subroutine straight_segment
!***************************************************************************
