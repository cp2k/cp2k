  SELECT CASE(lp)
 CASE(0)
 CALL integrate_core_0(grid(1,1,1),coef_xyz(1),pol_x(0,-cmax), pol_y(1,0,-cmax),pol_z(1,0,-cmax), &
                        map(-cmax,1),sphere_bounds(1),cmax,gridbounds(1,1))
 CASE(1)
 CALL integrate_core_1(grid(1,1,1),coef_xyz(1),pol_x(0,-cmax), pol_y(1,0,-cmax),pol_z(1,0,-cmax), &
                        map(-cmax,1),sphere_bounds(1),cmax,gridbounds(1,1))
 CASE(2)
 CALL integrate_core_2(grid(1,1,1),coef_xyz(1),pol_x(0,-cmax), pol_y(1,0,-cmax),pol_z(1,0,-cmax), &
                        map(-cmax,1),sphere_bounds(1),cmax,gridbounds(1,1))
 CASE(3)
 CALL integrate_core_3(grid(1,1,1),coef_xyz(1),pol_x(0,-cmax), pol_y(1,0,-cmax),pol_z(1,0,-cmax), &
                        map(-cmax,1),sphere_bounds(1),cmax,gridbounds(1,1))
 CASE(4)
 CALL integrate_core_4(grid(1,1,1),coef_xyz(1),pol_x(0,-cmax), pol_y(1,0,-cmax),pol_z(1,0,-cmax), &
                        map(-cmax,1),sphere_bounds(1),cmax,gridbounds(1,1))
 CASE(5)
 CALL integrate_core_5(grid(1,1,1),coef_xyz(1),pol_x(0,-cmax), pol_y(1,0,-cmax),pol_z(1,0,-cmax), &
                        map(-cmax,1),sphere_bounds(1),cmax,gridbounds(1,1))
 CASE(6)
 CALL integrate_core_6(grid(1,1,1),coef_xyz(1),pol_x(0,-cmax), pol_y(1,0,-cmax),pol_z(1,0,-cmax), &
                        map(-cmax,1),sphere_bounds(1),cmax,gridbounds(1,1))
 CASE(7)
 CALL integrate_core_7(grid(1,1,1),coef_xyz(1),pol_x(0,-cmax), pol_y(1,0,-cmax),pol_z(1,0,-cmax), &
                        map(-cmax,1),sphere_bounds(1),cmax,gridbounds(1,1))
 CASE(8)
 CALL integrate_core_8(grid(1,1,1),coef_xyz(1),pol_x(0,-cmax), pol_y(1,0,-cmax),pol_z(1,0,-cmax), &
                        map(-cmax,1),sphere_bounds(1),cmax,gridbounds(1,1))
 CASE(9)
 CALL integrate_core_9(grid(1,1,1),coef_xyz(1),pol_x(0,-cmax), pol_y(1,0,-cmax),pol_z(1,0,-cmax), &
                        map(-cmax,1),sphere_bounds(1),cmax,gridbounds(1,1))
  CASE DEFAULT
 CALL integrate_core_default(grid(1,1,1),coef_xyz(1),pol_x(0,-cmax),pol_y(1,0,-cmax),pol_z(1,0,-cmax), &
                        map(-cmax,1),sphere_bounds(1),lp,cmax,gridbounds(1,1))
  END SELECT
