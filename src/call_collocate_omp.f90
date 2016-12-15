  SELECT CASE(lp)
 CASE(0)
 CALL collocate_core_0(lgrid%r(1,ithread_l),coef_xyz(1),pol_x(0,-cmax), pol_y(1,0,-cmax),pol_z(1,0,-cmax), &
                        map(-cmax,1),sphere_bounds(1),cmax,gridbounds(1,1))
 CASE(1)
 CALL collocate_core_1(lgrid%r(1,ithread_l),coef_xyz(1),pol_x(0,-cmax), pol_y(1,0,-cmax),pol_z(1,0,-cmax), &
                        map(-cmax,1),sphere_bounds(1),cmax,gridbounds(1,1))
 CASE(2)
 CALL collocate_core_2(lgrid%r(1,ithread_l),coef_xyz(1),pol_x(0,-cmax), pol_y(1,0,-cmax),pol_z(1,0,-cmax), &
                        map(-cmax,1),sphere_bounds(1),cmax,gridbounds(1,1))
 CASE(3)
 CALL collocate_core_3(lgrid%r(1,ithread_l),coef_xyz(1),pol_x(0,-cmax), pol_y(1,0,-cmax),pol_z(1,0,-cmax), &
                        map(-cmax,1),sphere_bounds(1),cmax,gridbounds(1,1))
 CASE(4)
 CALL collocate_core_4(lgrid%r(1,ithread_l),coef_xyz(1),pol_x(0,-cmax), pol_y(1,0,-cmax),pol_z(1,0,-cmax), &
                        map(-cmax,1),sphere_bounds(1),cmax,gridbounds(1,1))
 CASE(5)
 CALL collocate_core_5(lgrid%r(1,ithread_l),coef_xyz(1),pol_x(0,-cmax), pol_y(1,0,-cmax),pol_z(1,0,-cmax), &
                        map(-cmax,1),sphere_bounds(1),cmax,gridbounds(1,1))
 CASE(6)
 CALL collocate_core_6(lgrid%r(1,ithread_l),coef_xyz(1),pol_x(0,-cmax), pol_y(1,0,-cmax),pol_z(1,0,-cmax), &
                        map(-cmax,1),sphere_bounds(1),cmax,gridbounds(1,1))
 CASE(7)
 CALL collocate_core_7(lgrid%r(1,ithread_l),coef_xyz(1),pol_x(0,-cmax), pol_y(1,0,-cmax),pol_z(1,0,-cmax), &
                        map(-cmax,1),sphere_bounds(1),cmax,gridbounds(1,1))
 CASE(8)
 CALL collocate_core_8(lgrid%r(1,ithread_l),coef_xyz(1),pol_x(0,-cmax), pol_y(1,0,-cmax),pol_z(1,0,-cmax), &
                        map(-cmax,1),sphere_bounds(1),cmax,gridbounds(1,1))
 CASE(9)
 CALL collocate_core_9(lgrid%r(1,ithread_l),coef_xyz(1),pol_x(0,-cmax), pol_y(1,0,-cmax),pol_z(1,0,-cmax), &
                        map(-cmax,1),sphere_bounds(1),cmax,gridbounds(1,1))
  CASE DEFAULT
 CALL collocate_core_default(lgrid%r(1,ithread_l),coef_xyz(1),pol_x(0,-cmax),pol_y(1,0,-cmax),pol_z(1,0,-cmax), &
                        map(-cmax,1),sphere_bounds(1),lp,cmax,gridbounds(1,1))
  END SELECT
