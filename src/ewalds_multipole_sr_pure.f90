                ! Compute the Short Range constribution according the task
                IF (debug_this_module) THEN
                   f         = HUGE(0.0_dp)
                   tij       = HUGE(0.0_dp)
                   tij_a     = HUGE(0.0_dp)
                   tij_ab    = HUGE(0.0_dp)
                   tij_abc   = HUGE(0.0_dp)
                   tij_abcd  = HUGE(0.0_dp)
                   tij_abcde = HUGE(0.0_dp)
                END IF
                r     = SQRT(rab2)
                irab2 = 1.0_dp/rab2
                ir    = 1.0_dp/r

                ! Compute the radial function
#if 0
                ! code for point multipole with screening
                IF (debug_this_module.AND.debug_r_space.AND.(.NOT.debug_g_space)) THEN
                   f(0)  = ir
                   tmp   = 0.0_dp
                ELSE
                   f(0)  = erfc(alpha*r)*ir
                   tmp   = EXP(-alpha**2*rab2)*oorootpi
                END IF
                fac = 1.0_dp
                DO i = 1, 5
                   fac  = fac*REAL(2*i-1,KIND=dp)
                   f(i) = irab2*(f(i-1)+ tmp*((2.0_dp*alpha**2)**i)/(fac*alpha))
                END DO
#endif
#if 0
                ! code for gaussian multipole with screening
                IF (debug_this_module.AND.debug_r_space.AND.(.NOT.debug_g_space)) THEN
                   f(0)  = ir
                   tmp1   = 0.0_dp
                   tmp2   = 0.0_dp
                ELSE
                   f(0)  = erf(beta*r)*ir - erf(alpha*r)*ir
                   tmp1   = EXP(-alpha**2*rab2)*oorootpi
                   tmp2   = EXP(-beta**2*rab2)*oorootpi
                END IF
                fac = 1.0_dp
                DO i = 1, 5
                   fac  = fac*REAL(2*i-1,KIND=dp)
                   f(i) = irab2*(f(i-1) + tmp1*((2.0_dp*alpha**2)**i)/(fac*alpha) - tmp2*((2.0_dp*beta**2)**i)/(fac*beta))
                END DO
#endif
#if 0
                IF (debug_this_module.AND.debug_r_space.AND.(.NOT.debug_g_space)) THEN
                   f(0)  = ir
                   tmp   = 0.0_dp
                ELSE
                   f(0)  = erf(alpha*r)*ir
                   tmp   = EXP(-alpha**2*rab2)*oorootpi
                END IF
                fac = 1.0_dp
                DO i = 1, 5
                   fac  = fac*REAL(2*i-1,KIND=dp)
                   f(i) = irab2*(f(i-1) - tmp*((2.0_dp*alpha**2)**i)/(fac*alpha))
                END DO
#endif
#if 1
                ! code for point multipole without screening
                f(0)  = ir
                DO i = 1, 5
                   f(i) = irab2*f(i-1)
                END DO
#endif
#if 0
                ! code for gaussian multipole without screening
                IF (debug_this_module.AND.debug_r_space.AND.(.NOT.debug_g_space)) THEN
                   f(0)  = ir
                   tmp   = 0.0_dp
                ELSE
                   f(0)  = erf(beta*r)*ir
                   tmp   = EXP(-beta**2*rab2)*oorootpi
                END IF
                fac = 1.0_dp
                PRINT *, "CHECK", f(0)
                DO i = 1, 5
                   fac  = fac*REAL(2*i-1,KIND=dp)
                   f(i) = irab2*(f(i-1) - tmp*((2.0_dp*beta**2)**i)/(fac*beta))
                END DO
#endif
                ! Compute the Tensor components
                force_eval = do_stress
                IF (task(1,1)) THEN
                   tij         = f(0)*fac_ij
                                                 force_eval = do_forces .OR.do_efield1
                END IF
                IF (task(2,2))                   force_eval = force_eval.OR.do_efield0
                IF (task(1,2).OR.force_eval) THEN
                   force_eval = do_stress
                   tij_a    = - rab*f(1)*fac_ij
                   IF (task(1,2))                force_eval = force_eval.OR. do_forces
                END IF
                IF (task(1,1))                   force_eval = force_eval.OR.do_efield2
                IF (task(3,3))                   force_eval = force_eval.OR.do_efield0
                IF (task(2,2).OR.task(3,1).OR.force_eval) THEN
                   force_eval = do_stress
                   DO b = 1,3
                      DO a = 1,3
                         tmp = rab(a)*rab(b)*fac_ij
                         tij_ab(a,b) = 3.0_dp*tmp*f(2)
                         IF (a==b) tij_ab(a,b) = tij_ab(a,b) - f(1)*fac_ij
                      END DO
                   END DO
                   IF (task(2,2).OR.task(3,1))   force_eval = force_eval.OR. do_forces
                END IF
                IF (task(2,2))                   force_eval = force_eval.OR.do_efield2
                IF (task(3,3))                   force_eval = force_eval.OR.do_efield1
                IF (task(3,2).OR.force_eval) THEN
                   force_eval = do_stress
                   DO c = 1, 3
                      DO b = 1, 3
                         DO a = 1, 3
                            tmp = rab(a)*rab(b)*rab(c)*fac_ij
                            tij_abc(a,b,c) = - 15.0_dp*tmp*f(3)
                            tmp = 3.0_dp*f(2)*fac_ij
                            IF (a==b) tij_abc(a,b,c) = tij_abc(a,b,c) + tmp*rab(c)
                            IF (a==c) tij_abc(a,b,c) = tij_abc(a,b,c) + tmp*rab(b)
                            IF (b==c) tij_abc(a,b,c) = tij_abc(a,b,c) + tmp*rab(a)
                         END DO
                      END DO
                   END DO
                   IF (task(3,2))                force_eval = force_eval.OR. do_forces
                END IF
                IF (task(3,3).OR.force_eval) THEN
                   force_eval = do_stress
                   DO d = 1, 3
                      DO c = 1, 3
                         DO b = 1, 3
                            DO a = 1, 3
                               tmp = rab(a)*rab(b)*rab(c)*rab(d)*fac_ij
                               tij_abcd(a,b,c,d) = 105.0_dp*tmp*f(4)
                               tmp1 = 15.0_dp*f(3)*fac_ij
                               tmp2 =  3.0_dp*f(2)*fac_ij
                               IF (a==b) THEN
                                            tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) - tmp1*rab(c)*rab(d)
                                  IF (c==d) tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) + tmp2
                               END IF
                               IF (a==c) THEN
                                            tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) - tmp1*rab(b)*rab(d)
                                  IF (b==d) tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) + tmp2
                               END IF
                               IF (a==d)    tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) - tmp1*rab(b)*rab(c)
                               IF (b==c) THEN
                                            tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) - tmp1*rab(a)*rab(d)
                                  IF (a==d) tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) + tmp2
                               END IF
                               IF (b==d)    tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) - tmp1*rab(a)*rab(c)
                               IF (c==d)    tij_abcd(a,b,c,d) = tij_abcd(a,b,c,d) - tmp1*rab(a)*rab(b)
                            END DO
                         END DO
                      END DO
                   END DO
                   IF (task(3,3))                force_eval = force_eval.OR. do_forces
                END IF
                IF (force_eval) THEN
                   force_eval = do_stress
                   DO e = 1, 3
                      DO d = 1, 3
                         DO c = 1, 3
                            DO b = 1, 3
                               DO a = 1, 3
                                  tmp = rab(a)*rab(b)*rab(c)*rab(d)*rab(e)*fac_ij
                                  tij_abcde(a,b,c,d,e) = -945.0_dp*tmp*f(5)
                                  tmp1 = 105.0_dp*f(4)*fac_ij
                                  tmp2 =  15.0_dp*f(3)*fac_ij
                                  IF (a==b) THEN
                                               tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(c)*rab(d)*rab(e)
                                     IF (c==d) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(e)
                                     IF (c==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(d)
                                     IF (d==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(c)
                                  END IF
                                  IF (a==c) THEN
                                               tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(b)*rab(d)*rab(e)
                                     IF (b==d) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(e)
                                     IF (b==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(d)
                                     IF (d==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(b)
                                  END IF
                                  IF (a==d) THEN
                                               tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(b)*rab(c)*rab(e)
                                     IF (b==c) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(e)
                                     IF (b==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(c)
                                     IF (c==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(b)
                                  END IF
                                  IF (a==e) THEN
                                               tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(b)*rab(c)*rab(d)
                                     IF (b==c) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(d)
                                     IF (b==d) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(c)
                                     IF (c==d) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(b)
                                  END IF
                                  IF (b==c) THEN
                                               tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(a)*rab(d)*rab(e)
                                     IF (d==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(a)
                                  END IF
                                  IF (b==d) THEN
                                               tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(a)*rab(c)*rab(e)
                                     IF (c==e) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(a)
                                  END IF
                                  IF (b==e) THEN
                                               tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(a)*rab(c)*rab(d)
                                     IF (c==d) tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) - tmp2*rab(a)
                                  END IF
                                  IF (c==d)    tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(a)*rab(b)*rab(e)
                                  IF (c==e)    tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(a)*rab(b)*rab(d)
                                  IF (d==e)    tij_abcde(a,b,c,d,e) = tij_abcde(a,b,c,d,e) + tmp1*rab(a)*rab(b)*rab(c)
                               END DO
                            END DO
                         END DO
                      END DO
                   END DO
                END IF
                eloc  = 0.0_dp
                fr    = 0.0_dp
                ef0_i = 0.0_dp
                ef0_j = 0.0_dp
                ef1_j = 0.0_dp
                ef1_i = 0.0_dp
                ef2_j = 0.0_dp
                ef2_i = 0.0_dp

#if 0

                ! Initialize the damping function.
                IF (kind_a==ikind) THEN
                   ! for atom i
                   SELECT CASE (itype_ij)
                   CASE (tang_toennies)
                      dampsumfi = 1.0_dp
                      xf = 1.0_dp
                      factorial = 1.0_dp
                      DO kk = 1, nkdamp_ij
                         xf = xf*dampa_ij*r
                         factorial = factorial * float(kk)
                         dampsumfi = dampsumfi + (xf/factorial)
                      END DO
                      dampaexpi = dexp(-dampa_ij * r)
                      dampfunci = dampsumfi * dampaexpi * dampfac_ij
                      dampfuncdiffi = -dampa_ij * dampaexpi * &
                                      dampfac_ij * (((dampa_ij * r) ** nkdamp_ij) / &
                                      factorial)
                   CASE DEFAULT
                      dampfunci=0.0_dp
                      dampfuncdiffi=0.0_dp
                   END SELECT

                   ! for atom j
                   SELECT CASE (itype_ji)
                   CASE (tang_toennies)
                      dampsumfj = 1.0_dp
                      xf = 1.0_dp
                      factorial = 1.0_dp
                      DO kk = 1, nkdamp_ji
                         xf = xf*dampa_ji*r
                         factorial = factorial * float(kk)
                         dampsumfj = dampsumfj + (xf/factorial)
                      END DO
                      dampaexpj = dexp(-dampa_ji * r)
                      dampfuncj = dampsumfj * dampaexpj * dampfac_ji
                      dampfuncdiffj = -dampa_ji * dampaexpj * &
                                      dampfac_ji * (((dampa_ji * r) ** nkdamp_ji) / &
                                      factorial)
                   CASE DEFAULT
                      dampfuncj = 0.0_dp
                      dampfuncdiffj = 0.0_dp
                   END SELECT
                ELSE
                   SELECT CASE (itype_ij)
                   CASE(tang_toennies)
                      dampsumfj = 1.0_dp
                      xf = 1.0_dp
                      factorial = 1.0_dp
                      DO kk = 1, nkdamp_ij
                         xf = xf*dampa_ij*r
                         factorial = factorial * float(kk)
                         dampsumfj = dampsumfj + (xf/factorial)
                      END DO
                      dampaexpj = dexp(-dampa_ij * r)
                      dampfuncj = dampsumfj * dampaexpj * dampfac_ij
                      dampfuncdiffj = -dampa_ij * dampaexpj * &
                                      dampfac_ij * (((dampa_ij * r) ** nkdamp_ij) / &
                                      factorial)
                   CASE DEFAULT
                      dampfuncj=0.0_dp
                      dampfuncdiffj=0.0_dp
                   END SELECT

                   !for j
                   SELECT CASE (itype_ji)
                   CASE (tang_toennies)
                      dampsumfi = 1.0_dp
                      xf = 1.0_dp
                      factorial = 1.0_dp
                      DO kk = 1, nkdamp_ji
                         xf = xf*dampa_ji*r
                         factorial = factorial * float(kk)
                         dampsumfi = dampsumfi + (xf/factorial)
                      END DO
                      dampaexpi = dexp(-dampa_ji * r)
                      dampfunci = dampsumfi * dampaexpi * dampfac_ji
                      dampfuncdiffi = -dampa_ji * dampaexpi * &
                                      dampfac_ji * (((dampa_ji * r) ** nkdamp_ji) / &
                                      factorial)
                   CASE DEFAULT
                      dampfunci = 0.0_dp
                      dampfuncdiffi = 0.0_dp
                   END SELECT
                END IF

                damptij_a = -rab*dampfunci*fac_ij*irab2*ir
                damptji_a = -rab*dampfuncj*fac_ij*irab2*ir
                DO b = 1,3
                   DO a = 1,3
                      tmp = rab(a)*rab(b)*fac_ij
                      damptij_ab(a,b) = tmp*(-dampfuncdiffi*irab2*irab2+3.0_dp*dampfunci*irab2*irab2*ir)
                      damptji_ab(a,b) = tmp*(-dampfuncdiffj*irab2*irab2+3.0_dp*dampfuncj*irab2*irab2*ir)
                      IF (a==b) damptij_ab(a,b) = damptij_ab(a,b) - dampfunci*fac_ij*irab2*ir
                      IF (a==b) damptji_ab(a,b) = damptji_ab(a,b) - dampfuncj*fac_ij*irab2*ir
                   END DO
                END DO

#endif

                ! Initialize the charge, dipole and quadrupole for atom A and B
                IF (debug_this_module) THEN
                   ch_j  = HUGE(0.0_dp)
                   ch_i  = HUGE(0.0_dp)
                   dp_j  = HUGE(0.0_dp)
                   dp_i  = HUGE(0.0_dp)
                   qp_j  = HUGE(0.0_dp)
                   qp_i  = HUGE(0.0_dp)
                END IF
                IF (ANY(task(1,:))) THEN
                   ch_j  = charges(atom_a)
                   ch_i  = charges(atom_b)
                END IF
                IF (ANY(task(2,:))) THEN
                   dp_j  = dipoles(:,atom_a)
                   dp_i  = dipoles(:,atom_b)
                END IF
                IF (ANY(task(3,:))) THEN
                   qp_j  = quadrupoles(:,:,atom_a)
                   qp_i  = quadrupoles(:,:,atom_b)
                END IF
                IF (task(1,1)) THEN
                   ! Charge - Charge
                   eloc = eloc + ch_i*tij*ch_j
                   ! Forces on particle i (locally b)
                   IF (do_forces.OR.do_stress) THEN
                      fr(1) = fr(1) - ch_j * tij_a(1) * ch_i
                      fr(2) = fr(2) - ch_j * tij_a(2) * ch_i
                      fr(3) = fr(3) - ch_j * tij_a(3) * ch_i
                   END IF
                   ! Electric fields
                   IF (do_efield) THEN
                      ! Potential
                      IF (do_efield0) THEN
                         ef0_i = ef0_i + tij * ch_j

                         ef0_j = ef0_j + tij * ch_i
                      END IF
                      ! Electric field
                      IF (do_efield1) THEN
                         ef1_i(1) = ef1_i(1) - tij_a(1) * ch_j
                         ef1_i(2) = ef1_i(2) - tij_a(2) * ch_j
                         ef1_i(3) = ef1_i(3) - tij_a(3) * ch_j

                         ef1_j(1) = ef1_j(1) + tij_a(1) * ch_i
                         ef1_j(2) = ef1_j(2) + tij_a(2) * ch_i
                         ef1_j(3) = ef1_j(3) + tij_a(3) * ch_i

#if 0
                         ef1_i(1) = ef1_i(1) + damptij_a(1) * ch_j
                         ef1_i(2) = ef1_i(2) + damptij_a(2) * ch_j
                         ef1_i(3) = ef1_i(3) + damptij_a(3) * ch_j

                         ef1_j(1) = ef1_j(1) - damptji_a(1) * ch_i
                         ef1_j(2) = ef1_j(2) - damptji_a(2) * ch_i
                         ef1_j(3) = ef1_j(3) - damptji_a(3) * ch_i
#endif

                      END IF
                      ! Electric field gradient
                      IF (do_efield2) THEN
                         ef2_i(1,1) = ef2_i(1,1) - tij_ab(1,1) * ch_j
                         ef2_i(2,1) = ef2_i(2,1) - tij_ab(2,1) * ch_j
                         ef2_i(3,1) = ef2_i(3,1) - tij_ab(3,1) * ch_j
                         ef2_i(1,2) = ef2_i(1,2) - tij_ab(1,2) * ch_j
                         ef2_i(2,2) = ef2_i(2,2) - tij_ab(2,2) * ch_j
                         ef2_i(3,2) = ef2_i(3,2) - tij_ab(3,2) * ch_j
                         ef2_i(1,3) = ef2_i(1,3) - tij_ab(1,3) * ch_j
                         ef2_i(2,3) = ef2_i(2,3) - tij_ab(2,3) * ch_j
                         ef2_i(3,3) = ef2_i(3,3) - tij_ab(3,3) * ch_j

                         ef2_j(1,1) = ef2_j(1,1) - tij_ab(1,1) * ch_i
                         ef2_j(2,1) = ef2_j(2,1) - tij_ab(2,1) * ch_i
                         ef2_j(3,1) = ef2_j(3,1) - tij_ab(3,1) * ch_i
                         ef2_j(1,2) = ef2_j(1,2) - tij_ab(1,2) * ch_i
                         ef2_j(2,2) = ef2_j(2,2) - tij_ab(2,2) * ch_i
                         ef2_j(3,2) = ef2_j(3,2) - tij_ab(3,2) * ch_i
                         ef2_j(1,3) = ef2_j(1,3) - tij_ab(1,3) * ch_i
                         ef2_j(2,3) = ef2_j(2,3) - tij_ab(2,3) * ch_i
                         ef2_j(3,3) = ef2_j(3,3) - tij_ab(3,3) * ch_i
                      END IF
                   END IF
                END IF
                IF (task(2,2)) THEN
                   ! Dipole - Dipole
                   tmp= - (dp_i(1)*(tij_ab(1,1)*dp_j(1)+&
                                    tij_ab(2,1)*dp_j(2)+&
                                    tij_ab(3,1)*dp_j(3))+&
                           dp_i(2)*(tij_ab(1,2)*dp_j(1)+&
                                    tij_ab(2,2)*dp_j(2)+&
                                    tij_ab(3,2)*dp_j(3))+&
                           dp_i(3)*(tij_ab(1,3)*dp_j(1)+&
                                    tij_ab(2,3)*dp_j(2)+&
                                    tij_ab(3,3)*dp_j(3)))
                   eloc = eloc + tmp
                   ! Forces on particle i (locally b)
                   IF (do_forces.OR.do_stress) THEN
                      DO k = 1, 3
                         fr(k) = fr(k) +  dp_i(1)*(tij_abc(1,1,k)*dp_j(1)+&
                                                   tij_abc(2,1,k)*dp_j(2)+&
                                                   tij_abc(3,1,k)*dp_j(3))&
                                       +  dp_i(2)*(tij_abc(1,2,k)*dp_j(1)+&
                                                   tij_abc(2,2,k)*dp_j(2)+&
                                                   tij_abc(3,2,k)*dp_j(3))&
                                       +  dp_i(3)*(tij_abc(1,3,k)*dp_j(1)+&
                                                   tij_abc(2,3,k)*dp_j(2)+&
                                                   tij_abc(3,3,k)*dp_j(3))
                      END DO
                   END IF
                   ! Electric fields
                   IF (do_efield) THEN
                      ! Potential
                      IF (do_efield0) THEN
                         ef0_i = ef0_i - (tij_a(1)*dp_j(1)+&
                                          tij_a(2)*dp_j(2)+&
                                          tij_a(3)*dp_j(3))

                         ef0_j = ef0_j + (tij_a(1)*dp_i(1)+&
                                          tij_a(2)*dp_i(2)+&
                                          tij_a(3)*dp_i(3))
                      END IF
                      ! Electric field
                      IF (do_efield1) THEN
                         ef1_i(1) = ef1_i(1) + (tij_ab(1,1)*dp_j(1)+&
                                                tij_ab(2,1)*dp_j(2)+&
                                                tij_ab(3,1)*dp_j(3))
                         ef1_i(2) = ef1_i(2) + (tij_ab(1,2)*dp_j(1)+&
                                                tij_ab(2,2)*dp_j(2)+&
                                                tij_ab(3,2)*dp_j(3))
                         ef1_i(3) = ef1_i(3) + (tij_ab(1,3)*dp_j(1)+&
                                                tij_ab(2,3)*dp_j(2)+&
                                                tij_ab(3,3)*dp_j(3))

                         ef1_j(1) = ef1_j(1) + (tij_ab(1,1)*dp_i(1)+&
                                                tij_ab(2,1)*dp_i(2)+&
                                                tij_ab(3,1)*dp_i(3))
                         ef1_j(2) = ef1_j(2) + (tij_ab(1,2)*dp_i(1)+&
                                                tij_ab(2,2)*dp_i(2)+&
                                                tij_ab(3,2)*dp_i(3))
                         ef1_j(3) = ef1_j(3) + (tij_ab(1,3)*dp_i(1)+&
                                                tij_ab(2,3)*dp_i(2)+&
                                                tij_ab(3,3)*dp_i(3))
                      END IF
                      ! Electric field gradient
                      IF (do_efield2) THEN
                         ef2_i(1,1) = ef2_i(1,1) + (tij_abc(1,1,1)*dp_j(1)+&
                                                    tij_abc(2,1,1)*dp_j(2)+&
                                                    tij_abc(3,1,1)*dp_j(3))
                         ef2_i(1,2) = ef2_i(1,2) + (tij_abc(1,1,2)*dp_j(1)+&
                                                    tij_abc(2,1,2)*dp_j(2)+&
                                                    tij_abc(3,1,2)*dp_j(3))
                         ef2_i(1,3) = ef2_i(1,3) + (tij_abc(1,1,3)*dp_j(1)+&
                                                    tij_abc(2,1,3)*dp_j(2)+&
                                                    tij_abc(3,1,3)*dp_j(3))
                         ef2_i(2,1) = ef2_i(2,1) + (tij_abc(1,2,1)*dp_j(1)+&
                                                    tij_abc(2,2,1)*dp_j(2)+&
                                                    tij_abc(3,2,1)*dp_j(3))
                         ef2_i(2,2) = ef2_i(2,2) + (tij_abc(1,2,2)*dp_j(1)+&
                                                    tij_abc(2,2,2)*dp_j(2)+&
                                                    tij_abc(3,2,2)*dp_j(3))
                         ef2_i(2,3) = ef2_i(2,3) + (tij_abc(1,2,3)*dp_j(1)+&
                                                    tij_abc(2,2,3)*dp_j(2)+&
                                                    tij_abc(3,2,3)*dp_j(3))
                         ef2_i(3,1) = ef2_i(3,1) + (tij_abc(1,3,1)*dp_j(1)+&
                                                    tij_abc(2,3,1)*dp_j(2)+&
                                                    tij_abc(3,3,1)*dp_j(3))
                         ef2_i(3,2) = ef2_i(3,2) + (tij_abc(1,3,2)*dp_j(1)+&
                                                    tij_abc(2,3,2)*dp_j(2)+&
                                                    tij_abc(3,3,2)*dp_j(3))
                         ef2_i(3,3) = ef2_i(3,3) + (tij_abc(1,3,3)*dp_j(1)+&
                                                    tij_abc(2,3,3)*dp_j(2)+&
                                                    tij_abc(3,3,3)*dp_j(3))

                         ef2_j(1,1) = ef2_j(1,1) - (tij_abc(1,1,1)*dp_i(1)+&
                                                    tij_abc(2,1,1)*dp_i(2)+&
                                                    tij_abc(3,1,1)*dp_i(3))
                         ef2_j(1,2) = ef2_j(1,2) - (tij_abc(1,1,2)*dp_i(1)+&
                                                    tij_abc(2,1,2)*dp_i(2)+&
                                                    tij_abc(3,1,2)*dp_i(3))
                         ef2_j(1,3) = ef2_j(1,3) - (tij_abc(1,1,3)*dp_i(1)+&
                                                    tij_abc(2,1,3)*dp_i(2)+&
                                                    tij_abc(3,1,3)*dp_i(3))
                         ef2_j(2,1) = ef2_j(2,1) - (tij_abc(1,2,1)*dp_i(1)+&
                                                    tij_abc(2,2,1)*dp_i(2)+&
                                                    tij_abc(3,2,1)*dp_i(3))
                         ef2_j(2,2) = ef2_j(2,2) - (tij_abc(1,2,2)*dp_i(1)+&
                                                    tij_abc(2,2,2)*dp_i(2)+&
                                                    tij_abc(3,2,2)*dp_i(3))
                         ef2_j(2,3) = ef2_j(2,3) - (tij_abc(1,2,3)*dp_i(1)+&
                                                    tij_abc(2,2,3)*dp_i(2)+&
                                                    tij_abc(3,2,3)*dp_i(3))
                         ef2_j(3,1) = ef2_j(3,1) - (tij_abc(1,3,1)*dp_i(1)+&
                                                    tij_abc(2,3,1)*dp_i(2)+&
                                                    tij_abc(3,3,1)*dp_i(3))
                         ef2_j(3,2) = ef2_j(3,2) - (tij_abc(1,3,2)*dp_i(1)+&
                                                    tij_abc(2,3,2)*dp_i(2)+&
                                                    tij_abc(3,3,2)*dp_i(3))
                         ef2_j(3,3) = ef2_j(3,3) - (tij_abc(1,3,3)*dp_i(1)+&
                                                    tij_abc(2,3,3)*dp_i(2)+&
                                                    tij_abc(3,3,3)*dp_i(3))
                      END IF
                   END IF
                END IF
                IF (task(2,1)) THEN
                   ! Dipole - Charge
                   tmp=   ch_j*(tij_a(1)*dp_i(1)+&
                                tij_a(2)*dp_i(2)+&
                                tij_a(3)*dp_i(3))&
                        - ch_i*(tij_a(1)*dp_j(1)+&
                                tij_a(2)*dp_j(2)+&
                                tij_a(3)*dp_j(3))
#if 0
                   tmp=  tmp- ch_j*(damptij_a(1)*dp_i(1)+&
                                damptij_a(2)*dp_i(2)+&
                                damptij_a(3)*dp_i(3))&
                        + ch_i*(damptji_a(1)*dp_j(1)+&
                                damptji_a(2)*dp_j(2)+&
                                damptji_a(3)*dp_j(3))
#endif
                   eloc = eloc + tmp
                   ! Forces on particle i (locally b)
                   IF (do_forces.OR.do_stress) THEN
                      DO k = 1, 3
                         fr(k) = fr(k) -  ch_j *(tij_ab(1,k)*dp_i(1)+&
                                                 tij_ab(2,k)*dp_i(2)+&
                                                 tij_ab(3,k)*dp_i(3))&
                                       +  ch_i *(tij_ab(1,k)*dp_j(1)+&
                                                 tij_ab(2,k)*dp_j(2)+&
                                                 tij_ab(3,k)*dp_j(3))
#if 0
                         fr(k) = fr(k) +  ch_j *(damptij_ab(1,k)*dp_i(1)+&
                                                 damptij_ab(2,k)*dp_i(2)+&
                                                 damptij_ab(3,k)*dp_i(3))&
                                       -  ch_i *(damptji_ab(1,k)*dp_j(1)+&
                                                 damptji_ab(2,k)*dp_j(2)+&
                                                 damptji_ab(3,k)*dp_j(3))
#endif
                      END DO
                   END IF
                END IF
                IF (task(3,3)) THEN
                   ! Quadrupole - Quadrupole
                   fac  = 1.0_dp/9.0_dp
                   tmp11 = qp_i(1,1)*(tij_abcd(1,1,1,1)*qp_j(1,1)+&
                                      tij_abcd(2,1,1,1)*qp_j(2,1)+&
                                      tij_abcd(3,1,1,1)*qp_j(3,1)+&
                                      tij_abcd(1,2,1,1)*qp_j(1,2)+&
                                      tij_abcd(2,2,1,1)*qp_j(2,2)+&
                                      tij_abcd(3,2,1,1)*qp_j(3,2)+&
                                      tij_abcd(1,3,1,1)*qp_j(1,3)+&
                                      tij_abcd(2,3,1,1)*qp_j(2,3)+&
                                      tij_abcd(3,3,1,1)*qp_j(3,3))
                   tmp21 = qp_i(2,1)*(tij_abcd(1,1,1,2)*qp_j(1,1)+&
                                      tij_abcd(2,1,1,2)*qp_j(2,1)+&
                                      tij_abcd(3,1,1,2)*qp_j(3,1)+&
                                      tij_abcd(1,2,1,2)*qp_j(1,2)+&
                                      tij_abcd(2,2,1,2)*qp_j(2,2)+&
                                      tij_abcd(3,2,1,2)*qp_j(3,2)+&
                                      tij_abcd(1,3,1,2)*qp_j(1,3)+&
                                      tij_abcd(2,3,1,2)*qp_j(2,3)+&
                                      tij_abcd(3,3,1,2)*qp_j(3,3))
                   tmp31 = qp_i(3,1)*(tij_abcd(1,1,1,3)*qp_j(1,1)+&
                                      tij_abcd(2,1,1,3)*qp_j(2,1)+&
                                      tij_abcd(3,1,1,3)*qp_j(3,1)+&
                                      tij_abcd(1,2,1,3)*qp_j(1,2)+&
                                      tij_abcd(2,2,1,3)*qp_j(2,2)+&
                                      tij_abcd(3,2,1,3)*qp_j(3,2)+&
                                      tij_abcd(1,3,1,3)*qp_j(1,3)+&
                                      tij_abcd(2,3,1,3)*qp_j(2,3)+&
                                      tij_abcd(3,3,1,3)*qp_j(3,3))
                   tmp22 = qp_i(2,2)*(tij_abcd(1,1,2,2)*qp_j(1,1)+&
                                      tij_abcd(2,1,2,2)*qp_j(2,1)+&
                                      tij_abcd(3,1,2,2)*qp_j(3,1)+&
                                      tij_abcd(1,2,2,2)*qp_j(1,2)+&
                                      tij_abcd(2,2,2,2)*qp_j(2,2)+&
                                      tij_abcd(3,2,2,2)*qp_j(3,2)+&
                                      tij_abcd(1,3,2,2)*qp_j(1,3)+&
                                      tij_abcd(2,3,2,2)*qp_j(2,3)+&
                                      tij_abcd(3,3,2,2)*qp_j(3,3))
                   tmp32 = qp_i(3,2)*(tij_abcd(1,1,2,3)*qp_j(1,1)+&
                                      tij_abcd(2,1,2,3)*qp_j(2,1)+&
                                      tij_abcd(3,1,2,3)*qp_j(3,1)+&
                                      tij_abcd(1,2,2,3)*qp_j(1,2)+&
                                      tij_abcd(2,2,2,3)*qp_j(2,2)+&
                                      tij_abcd(3,2,2,3)*qp_j(3,2)+&
                                      tij_abcd(1,3,2,3)*qp_j(1,3)+&
                                      tij_abcd(2,3,2,3)*qp_j(2,3)+&
                                      tij_abcd(3,3,2,3)*qp_j(3,3))
                   tmp33 = qp_i(3,3)*(tij_abcd(1,1,3,3)*qp_j(1,1)+&
                                      tij_abcd(2,1,3,3)*qp_j(2,1)+&
                                      tij_abcd(3,1,3,3)*qp_j(3,1)+&
                                      tij_abcd(1,2,3,3)*qp_j(1,2)+&
                                      tij_abcd(2,2,3,3)*qp_j(2,2)+&
                                      tij_abcd(3,2,3,3)*qp_j(3,2)+&
                                      tij_abcd(1,3,3,3)*qp_j(1,3)+&
                                      tij_abcd(2,3,3,3)*qp_j(2,3)+&
                                      tij_abcd(3,3,3,3)*qp_j(3,3))
                   tmp12 = tmp21
                   tmp13 = tmp31
                   tmp23 = tmp32
                   tmp   = tmp11 + tmp12 + tmp13 + &
                           tmp21 + tmp22 + tmp23 + &
                           tmp31 + tmp32 + tmp33

                   eloc = eloc + fac*tmp
                   ! Forces on particle i (locally b)
                   IF (do_forces.OR.do_stress) THEN
                      DO k = 1, 3
                         tmp11 = qp_i(1,1)*(tij_abcde(1,1,1,1,k)*qp_j(1,1)+&
                                            tij_abcde(2,1,1,1,k)*qp_j(2,1)+&
                                            tij_abcde(3,1,1,1,k)*qp_j(3,1)+&
                                            tij_abcde(1,2,1,1,k)*qp_j(1,2)+&
                                            tij_abcde(2,2,1,1,k)*qp_j(2,2)+&
                                            tij_abcde(3,2,1,1,k)*qp_j(3,2)+&
                                            tij_abcde(1,3,1,1,k)*qp_j(1,3)+&
                                            tij_abcde(2,3,1,1,k)*qp_j(2,3)+&
                                            tij_abcde(3,3,1,1,k)*qp_j(3,3))
                         tmp21 = qp_i(2,1)*(tij_abcde(1,1,2,1,k)*qp_j(1,1)+&
                                            tij_abcde(2,1,2,1,k)*qp_j(2,1)+&
                                            tij_abcde(3,1,2,1,k)*qp_j(3,1)+&
                                            tij_abcde(1,2,2,1,k)*qp_j(1,2)+&
                                            tij_abcde(2,2,2,1,k)*qp_j(2,2)+&
                                            tij_abcde(3,2,2,1,k)*qp_j(3,2)+&
                                            tij_abcde(1,3,2,1,k)*qp_j(1,3)+&
                                            tij_abcde(2,3,2,1,k)*qp_j(2,3)+&
                                            tij_abcde(3,3,2,1,k)*qp_j(3,3))
                         tmp31 = qp_i(3,1)*(tij_abcde(1,1,3,1,k)*qp_j(1,1)+&
                                            tij_abcde(2,1,3,1,k)*qp_j(2,1)+&
                                            tij_abcde(3,1,3,1,k)*qp_j(3,1)+&
                                            tij_abcde(1,2,3,1,k)*qp_j(1,2)+&
                                            tij_abcde(2,2,3,1,k)*qp_j(2,2)+&
                                            tij_abcde(3,2,3,1,k)*qp_j(3,2)+&
                                            tij_abcde(1,3,3,1,k)*qp_j(1,3)+&
                                            tij_abcde(2,3,3,1,k)*qp_j(2,3)+&
                                            tij_abcde(3,3,3,1,k)*qp_j(3,3))
                         tmp22 = qp_i(2,2)*(tij_abcde(1,1,2,2,k)*qp_j(1,1)+&
                                            tij_abcde(2,1,2,2,k)*qp_j(2,1)+&
                                            tij_abcde(3,1,2,2,k)*qp_j(3,1)+&
                                            tij_abcde(1,2,2,2,k)*qp_j(1,2)+&
                                            tij_abcde(2,2,2,2,k)*qp_j(2,2)+&
                                            tij_abcde(3,2,2,2,k)*qp_j(3,2)+&
                                            tij_abcde(1,3,2,2,k)*qp_j(1,3)+&
                                            tij_abcde(2,3,2,2,k)*qp_j(2,3)+&
                                            tij_abcde(3,3,2,2,k)*qp_j(3,3))
                         tmp32 = qp_i(3,2)*(tij_abcde(1,1,3,2,k)*qp_j(1,1)+&
                                            tij_abcde(2,1,3,2,k)*qp_j(2,1)+&
                                            tij_abcde(3,1,3,2,k)*qp_j(3,1)+&
                                            tij_abcde(1,2,3,2,k)*qp_j(1,2)+&
                                            tij_abcde(2,2,3,2,k)*qp_j(2,2)+&
                                            tij_abcde(3,2,3,2,k)*qp_j(3,2)+&
                                            tij_abcde(1,3,3,2,k)*qp_j(1,3)+&
                                            tij_abcde(2,3,3,2,k)*qp_j(2,3)+&
                                            tij_abcde(3,3,3,2,k)*qp_j(3,3))
                         tmp33 = qp_i(3,3)*(tij_abcde(1,1,3,3,k)*qp_j(1,1)+&
                                            tij_abcde(2,1,3,3,k)*qp_j(2,1)+&
                                            tij_abcde(3,1,3,3,k)*qp_j(3,1)+&
                                            tij_abcde(1,2,3,3,k)*qp_j(1,2)+&
                                            tij_abcde(2,2,3,3,k)*qp_j(2,2)+&
                                            tij_abcde(3,2,3,3,k)*qp_j(3,2)+&
                                            tij_abcde(1,3,3,3,k)*qp_j(1,3)+&
                                            tij_abcde(2,3,3,3,k)*qp_j(2,3)+&
                                            tij_abcde(3,3,3,3,k)*qp_j(3,3))
                         tmp12 = tmp21
                         tmp13 = tmp31
                         tmp23 = tmp32
                         fr(k) = fr(k) - fac * ( tmp11 + tmp12 + tmp13 +&
                                                 tmp21 + tmp22 + tmp23 +&
                                                 tmp31 + tmp32 + tmp33  )
                      END DO
                   END IF
                   ! Electric field
                   IF (do_efield) THEN
                      fac = 1.0_dp/3.0_dp
                      ! Potential
                      IF (do_efield0) THEN
                         ef0_i = ef0_i + fac*(tij_ab(1,1)*qp_j(1,1)+&
                                              tij_ab(2,1)*qp_j(2,1)+&
                                              tij_ab(3,1)*qp_j(3,1)+&
                                              tij_ab(1,2)*qp_j(1,2)+&
                                              tij_ab(2,2)*qp_j(2,2)+&
                                              tij_ab(3,2)*qp_j(3,2)+&
                                              tij_ab(1,3)*qp_j(1,3)+&
                                              tij_ab(2,3)*qp_j(2,3)+&
                                              tij_ab(3,3)*qp_j(3,3))

                         ef0_j = ef0_j + fac*(tij_ab(1,1)*qp_i(1,1)+&
                                              tij_ab(2,1)*qp_i(2,1)+&
                                              tij_ab(3,1)*qp_i(3,1)+&
                                              tij_ab(1,2)*qp_i(1,2)+&
                                              tij_ab(2,2)*qp_i(2,2)+&
                                              tij_ab(3,2)*qp_i(3,2)+&
                                              tij_ab(1,3)*qp_i(1,3)+&
                                              tij_ab(2,3)*qp_i(2,3)+&
                                              tij_ab(3,3)*qp_i(3,3))
                      END IF
                      ! Electric field
                      IF (do_efield1) THEN
                         ef1_i(1) = ef1_i(1) - fac*(tij_abc(1,1,1)*qp_j(1,1)+&
                                                    tij_abc(2,1,1)*qp_j(2,1)+&
                                                    tij_abc(3,1,1)*qp_j(3,1)+&
                                                    tij_abc(1,2,1)*qp_j(1,2)+&
                                                    tij_abc(2,2,1)*qp_j(2,2)+&
                                                    tij_abc(3,2,1)*qp_j(3,2)+&
                                                    tij_abc(1,3,1)*qp_j(1,3)+&
                                                    tij_abc(2,3,1)*qp_j(2,3)+&
                                                    tij_abc(3,3,1)*qp_j(3,3))
                         ef1_i(2) = ef1_i(2) - fac*(tij_abc(1,1,2)*qp_j(1,1)+&
                                                    tij_abc(2,1,2)*qp_j(2,1)+&
                                                    tij_abc(3,1,2)*qp_j(3,1)+&
                                                    tij_abc(1,2,2)*qp_j(1,2)+&
                                                    tij_abc(2,2,2)*qp_j(2,2)+&
                                                    tij_abc(3,2,2)*qp_j(3,2)+&
                                                    tij_abc(1,3,2)*qp_j(1,3)+&
                                                    tij_abc(2,3,2)*qp_j(2,3)+&
                                                    tij_abc(3,3,2)*qp_j(3,3))
                         ef1_i(3) = ef1_i(3) - fac*(tij_abc(1,1,3)*qp_j(1,1)+&
                                                    tij_abc(2,1,3)*qp_j(2,1)+&
                                                    tij_abc(3,1,3)*qp_j(3,1)+&
                                                    tij_abc(1,2,3)*qp_j(1,2)+&
                                                    tij_abc(2,2,3)*qp_j(2,2)+&
                                                    tij_abc(3,2,3)*qp_j(3,2)+&
                                                    tij_abc(1,3,3)*qp_j(1,3)+&
                                                    tij_abc(2,3,3)*qp_j(2,3)+&
                                                    tij_abc(3,3,3)*qp_j(3,3))

                         ef1_j(1) = ef1_j(1) + fac*(tij_abc(1,1,1)*qp_i(1,1)+&
                                                    tij_abc(2,1,1)*qp_i(2,1)+&
                                                    tij_abc(3,1,1)*qp_i(3,1)+&
                                                    tij_abc(1,2,1)*qp_i(1,2)+&
                                                    tij_abc(2,2,1)*qp_i(2,2)+&
                                                    tij_abc(3,2,1)*qp_i(3,2)+&
                                                    tij_abc(1,3,1)*qp_i(1,3)+&
                                                    tij_abc(2,3,1)*qp_i(2,3)+&
                                                    tij_abc(3,3,1)*qp_i(3,3))
                         ef1_j(2) = ef1_j(2) + fac*(tij_abc(1,1,2)*qp_i(1,1)+&
                                                    tij_abc(2,1,2)*qp_i(2,1)+&
                                                    tij_abc(3,1,2)*qp_i(3,1)+&
                                                    tij_abc(1,2,2)*qp_i(1,2)+&
                                                    tij_abc(2,2,2)*qp_i(2,2)+&
                                                    tij_abc(3,2,2)*qp_i(3,2)+&
                                                    tij_abc(1,3,2)*qp_i(1,3)+&
                                                    tij_abc(2,3,2)*qp_i(2,3)+&
                                                    tij_abc(3,3,2)*qp_i(3,3))
                         ef1_j(3) = ef1_j(3) + fac*(tij_abc(1,1,3)*qp_i(1,1)+&
                                                    tij_abc(2,1,3)*qp_i(2,1)+&
                                                    tij_abc(3,1,3)*qp_i(3,1)+&
                                                    tij_abc(1,2,3)*qp_i(1,2)+&
                                                    tij_abc(2,2,3)*qp_i(2,2)+&
                                                    tij_abc(3,2,3)*qp_i(3,2)+&
                                                    tij_abc(1,3,3)*qp_i(1,3)+&
                                                    tij_abc(2,3,3)*qp_i(2,3)+&
                                                    tij_abc(3,3,3)*qp_i(3,3))
                      END IF
                      ! Electric field gradient
                      IF (do_efield2) THEN
                         tmp11 =   fac *(tij_abcd(1,1,1,1)*qp_j(1,1)+&
                                         tij_abcd(2,1,1,1)*qp_j(2,1)+&
                                         tij_abcd(3,1,1,1)*qp_j(3,1)+&
                                         tij_abcd(1,2,1,1)*qp_j(1,2)+&
                                         tij_abcd(2,2,1,1)*qp_j(2,2)+&
                                         tij_abcd(3,2,1,1)*qp_j(3,2)+&
                                         tij_abcd(1,3,1,1)*qp_j(1,3)+&
                                         tij_abcd(2,3,1,1)*qp_j(2,3)+&
                                         tij_abcd(3,3,1,1)*qp_j(3,3))
                         tmp12 =   fac *(tij_abcd(1,1,1,2)*qp_j(1,1)+&
                                         tij_abcd(2,1,1,2)*qp_j(2,1)+&
                                         tij_abcd(3,1,1,2)*qp_j(3,1)+&
                                         tij_abcd(1,2,1,2)*qp_j(1,2)+&
                                         tij_abcd(2,2,1,2)*qp_j(2,2)+&
                                         tij_abcd(3,2,1,2)*qp_j(3,2)+&
                                         tij_abcd(1,3,1,2)*qp_j(1,3)+&
                                         tij_abcd(2,3,1,2)*qp_j(2,3)+&
                                         tij_abcd(3,3,1,2)*qp_j(3,3))
                         tmp13 =   fac *(tij_abcd(1,1,1,3)*qp_j(1,1)+&
                                         tij_abcd(2,1,1,3)*qp_j(2,1)+&
                                         tij_abcd(3,1,1,3)*qp_j(3,1)+&
                                         tij_abcd(1,2,1,3)*qp_j(1,2)+&
                                         tij_abcd(2,2,1,3)*qp_j(2,2)+&
                                         tij_abcd(3,2,1,3)*qp_j(3,2)+&
                                         tij_abcd(1,3,1,3)*qp_j(1,3)+&
                                         tij_abcd(2,3,1,3)*qp_j(2,3)+&
                                         tij_abcd(3,3,1,3)*qp_j(3,3))
                         tmp22 =   fac *(tij_abcd(1,1,2,2)*qp_j(1,1)+&
                                         tij_abcd(2,1,2,2)*qp_j(2,1)+&
                                         tij_abcd(3,1,2,2)*qp_j(3,1)+&
                                         tij_abcd(1,2,2,2)*qp_j(1,2)+&
                                         tij_abcd(2,2,2,2)*qp_j(2,2)+&
                                         tij_abcd(3,2,2,2)*qp_j(3,2)+&
                                         tij_abcd(1,3,2,2)*qp_j(1,3)+&
                                         tij_abcd(2,3,2,2)*qp_j(2,3)+&
                                         tij_abcd(3,3,2,2)*qp_j(3,3))
                         tmp23 =   fac *(tij_abcd(1,1,2,3)*qp_j(1,1)+&
                                         tij_abcd(2,1,2,3)*qp_j(2,1)+&
                                         tij_abcd(3,1,2,3)*qp_j(3,1)+&
                                         tij_abcd(1,2,2,3)*qp_j(1,2)+&
                                         tij_abcd(2,2,2,3)*qp_j(2,2)+&
                                         tij_abcd(3,2,2,3)*qp_j(3,2)+&
                                         tij_abcd(1,3,2,3)*qp_j(1,3)+&
                                         tij_abcd(2,3,2,3)*qp_j(2,3)+&
                                         tij_abcd(3,3,2,3)*qp_j(3,3))
                         tmp33 =   fac *(tij_abcd(1,1,3,3)*qp_j(1,1)+&
                                         tij_abcd(2,1,3,3)*qp_j(2,1)+&
                                         tij_abcd(3,1,3,3)*qp_j(3,1)+&
                                         tij_abcd(1,2,3,3)*qp_j(1,2)+&
                                         tij_abcd(2,2,3,3)*qp_j(2,2)+&
                                         tij_abcd(3,2,3,3)*qp_j(3,2)+&
                                         tij_abcd(1,3,3,3)*qp_j(1,3)+&
                                         tij_abcd(2,3,3,3)*qp_j(2,3)+&
                                         tij_abcd(3,3,3,3)*qp_j(3,3))

                         ef2_i(1,1) = ef2_i(1,1) - tmp11
                         ef2_i(1,2) = ef2_i(1,2) - tmp12
                         ef2_i(1,3) = ef2_i(1,3) - tmp13
                         ef2_i(2,1) = ef2_i(2,1) - tmp12
                         ef2_i(2,2) = ef2_i(2,2) - tmp22
                         ef2_i(2,3) = ef2_i(2,3) - tmp23
                         ef2_i(3,1) = ef2_i(3,1) - tmp13
                         ef2_i(3,2) = ef2_i(3,2) - tmp23
                         ef2_i(3,3) = ef2_i(3,3) - tmp33

                         tmp11 =   fac *(tij_abcd(1,1,1,1)*qp_i(1,1)+&
                                         tij_abcd(2,1,1,1)*qp_i(2,1)+&
                                         tij_abcd(3,1,1,1)*qp_i(3,1)+&
                                         tij_abcd(1,2,1,1)*qp_i(1,2)+&
                                         tij_abcd(2,2,1,1)*qp_i(2,2)+&
                                         tij_abcd(3,2,1,1)*qp_i(3,2)+&
                                         tij_abcd(1,3,1,1)*qp_i(1,3)+&
                                         tij_abcd(2,3,1,1)*qp_i(2,3)+&
                                         tij_abcd(3,3,1,1)*qp_i(3,3))
                         tmp12 =   fac *(tij_abcd(1,1,1,2)*qp_i(1,1)+&
                                         tij_abcd(2,1,1,2)*qp_i(2,1)+&
                                         tij_abcd(3,1,1,2)*qp_i(3,1)+&
                                         tij_abcd(1,2,1,2)*qp_i(1,2)+&
                                         tij_abcd(2,2,1,2)*qp_i(2,2)+&
                                         tij_abcd(3,2,1,2)*qp_i(3,2)+&
                                         tij_abcd(1,3,1,2)*qp_i(1,3)+&
                                         tij_abcd(2,3,1,2)*qp_i(2,3)+&
                                         tij_abcd(3,3,1,2)*qp_i(3,3))
                         tmp13 =   fac *(tij_abcd(1,1,1,3)*qp_i(1,1)+&
                                         tij_abcd(2,1,1,3)*qp_i(2,1)+&
                                         tij_abcd(3,1,1,3)*qp_i(3,1)+&
                                         tij_abcd(1,2,1,3)*qp_i(1,2)+&
                                         tij_abcd(2,2,1,3)*qp_i(2,2)+&
                                         tij_abcd(3,2,1,3)*qp_i(3,2)+&
                                         tij_abcd(1,3,1,3)*qp_i(1,3)+&
                                         tij_abcd(2,3,1,3)*qp_i(2,3)+&
                                         tij_abcd(3,3,1,3)*qp_i(3,3))
                         tmp22 =   fac *(tij_abcd(1,1,2,2)*qp_i(1,1)+&
                                         tij_abcd(2,1,2,2)*qp_i(2,1)+&
                                         tij_abcd(3,1,2,2)*qp_i(3,1)+&
                                         tij_abcd(1,2,2,2)*qp_i(1,2)+&
                                         tij_abcd(2,2,2,2)*qp_i(2,2)+&
                                         tij_abcd(3,2,2,2)*qp_i(3,2)+&
                                         tij_abcd(1,3,2,2)*qp_i(1,3)+&
                                         tij_abcd(2,3,2,2)*qp_i(2,3)+&
                                         tij_abcd(3,3,2,2)*qp_i(3,3))
                         tmp23 =   fac *(tij_abcd(1,1,2,3)*qp_i(1,1)+&
                                         tij_abcd(2,1,2,3)*qp_i(2,1)+&
                                         tij_abcd(3,1,2,3)*qp_i(3,1)+&
                                         tij_abcd(1,2,2,3)*qp_i(1,2)+&
                                         tij_abcd(2,2,2,3)*qp_i(2,2)+&
                                         tij_abcd(3,2,2,3)*qp_i(3,2)+&
                                         tij_abcd(1,3,2,3)*qp_i(1,3)+&
                                         tij_abcd(2,3,2,3)*qp_i(2,3)+&
                                         tij_abcd(3,3,2,3)*qp_i(3,3))
                         tmp33 =   fac *(tij_abcd(1,1,3,3)*qp_i(1,1)+&
                                         tij_abcd(2,1,3,3)*qp_i(2,1)+&
                                         tij_abcd(3,1,3,3)*qp_i(3,1)+&
                                         tij_abcd(1,2,3,3)*qp_i(1,2)+&
                                         tij_abcd(2,2,3,3)*qp_i(2,2)+&
                                         tij_abcd(3,2,3,3)*qp_i(3,2)+&
                                         tij_abcd(1,3,3,3)*qp_i(1,3)+&
                                         tij_abcd(2,3,3,3)*qp_i(2,3)+&
                                         tij_abcd(3,3,3,3)*qp_i(3,3))

                         ef2_j(1,1) = ef2_j(1,1) - tmp11
                         ef2_j(1,2) = ef2_j(1,2) - tmp12
                         ef2_j(1,3) = ef2_j(1,3) - tmp13
                         ef2_j(2,1) = ef2_j(2,1) - tmp12
                         ef2_j(2,2) = ef2_j(2,2) - tmp22
                         ef2_j(2,3) = ef2_j(2,3) - tmp23
                         ef2_j(3,1) = ef2_j(3,1) - tmp13
                         ef2_j(3,2) = ef2_j(3,2) - tmp23
                         ef2_j(3,3) = ef2_j(3,3) - tmp33
                      END IF
                   END IF
                END IF
                IF (task(3,2)) THEN
                   ! Quadrupole - Dipole
                   fac = 1.0_dp/3.0_dp
                   ! Dipole i (locally B) - Quadrupole j (locally A)
                   tmp_ij = dp_i(1)*(tij_abc(1,1,1)*qp_j(1,1)+&
                                     tij_abc(2,1,1)*qp_j(2,1)+&
                                     tij_abc(3,1,1)*qp_j(3,1)+&
                                     tij_abc(1,2,1)*qp_j(1,2)+&
                                     tij_abc(2,2,1)*qp_j(2,2)+&
                                     tij_abc(3,2,1)*qp_j(3,2)+&
                                     tij_abc(1,3,1)*qp_j(1,3)+&
                                     tij_abc(2,3,1)*qp_j(2,3)+&
                                     tij_abc(3,3,1)*qp_j(3,3))+&
                            dp_i(2)*(tij_abc(1,1,2)*qp_j(1,1)+&
                                     tij_abc(2,1,2)*qp_j(2,1)+&
                                     tij_abc(3,1,2)*qp_j(3,1)+&
                                     tij_abc(1,2,2)*qp_j(1,2)+&
                                     tij_abc(2,2,2)*qp_j(2,2)+&
                                     tij_abc(3,2,2)*qp_j(3,2)+&
                                     tij_abc(1,3,2)*qp_j(1,3)+&
                                     tij_abc(2,3,2)*qp_j(2,3)+&
                                     tij_abc(3,3,2)*qp_j(3,3))+&
                            dp_i(3)*(tij_abc(1,1,3)*qp_j(1,1)+&
                                     tij_abc(2,1,3)*qp_j(2,1)+&
                                     tij_abc(3,1,3)*qp_j(3,1)+&
                                     tij_abc(1,2,3)*qp_j(1,2)+&
                                     tij_abc(2,2,3)*qp_j(2,2)+&
                                     tij_abc(3,2,3)*qp_j(3,2)+&
                                     tij_abc(1,3,3)*qp_j(1,3)+&
                                     tij_abc(2,3,3)*qp_j(2,3)+&
                                     tij_abc(3,3,3)*qp_j(3,3))

                   ! Dipole j (locally A) - Quadrupole i (locally B)
                   tmp_ji = dp_j(1)*(tij_abc(1,1,1)*qp_i(1,1)+&
                                     tij_abc(2,1,1)*qp_i(2,1)+&
                                     tij_abc(3,1,1)*qp_i(3,1)+&
                                     tij_abc(1,2,1)*qp_i(1,2)+&
                                     tij_abc(2,2,1)*qp_i(2,2)+&
                                     tij_abc(3,2,1)*qp_i(3,2)+&
                                     tij_abc(1,3,1)*qp_i(1,3)+&
                                     tij_abc(2,3,1)*qp_i(2,3)+&
                                     tij_abc(3,3,1)*qp_i(3,3))+&
                            dp_j(2)*(tij_abc(1,1,2)*qp_i(1,1)+&
                                     tij_abc(2,1,2)*qp_i(2,1)+&
                                     tij_abc(3,1,2)*qp_i(3,1)+&
                                     tij_abc(1,2,2)*qp_i(1,2)+&
                                     tij_abc(2,2,2)*qp_i(2,2)+&
                                     tij_abc(3,2,2)*qp_i(3,2)+&
                                     tij_abc(1,3,2)*qp_i(1,3)+&
                                     tij_abc(2,3,2)*qp_i(2,3)+&
                                     tij_abc(3,3,2)*qp_i(3,3))+&
                            dp_j(3)*(tij_abc(1,1,3)*qp_i(1,1)+&
                                     tij_abc(2,1,3)*qp_i(2,1)+&
                                     tij_abc(3,1,3)*qp_i(3,1)+&
                                     tij_abc(1,2,3)*qp_i(1,2)+&
                                     tij_abc(2,2,3)*qp_i(2,2)+&
                                     tij_abc(3,2,3)*qp_i(3,2)+&
                                     tij_abc(1,3,3)*qp_i(1,3)+&
                                     tij_abc(2,3,3)*qp_i(2,3)+&
                                     tij_abc(3,3,3)*qp_i(3,3))

                   tmp= fac * (tmp_ij - tmp_ji)
                   eloc = eloc + tmp
                   IF (do_forces.OR.do_stress) THEN
                      DO k = 1, 3
                         ! Dipole i (locally B) - Quadrupole j (locally A)
                         tmp_ij = dp_i(1)*(tij_abcd(1,1,1,k)*qp_j(1,1)+&
                                           tij_abcd(2,1,1,k)*qp_j(2,1)+&
                                           tij_abcd(3,1,1,k)*qp_j(3,1)+&
                                           tij_abcd(1,2,1,k)*qp_j(1,2)+&
                                           tij_abcd(2,2,1,k)*qp_j(2,2)+&
                                           tij_abcd(3,2,1,k)*qp_j(3,2)+&
                                           tij_abcd(1,3,1,k)*qp_j(1,3)+&
                                           tij_abcd(2,3,1,k)*qp_j(2,3)+&
                                           tij_abcd(3,3,1,k)*qp_j(3,3))+&
                                  dp_i(2)*(tij_abcd(1,1,2,k)*qp_j(1,1)+&
                                           tij_abcd(2,1,2,k)*qp_j(2,1)+&
                                           tij_abcd(3,1,2,k)*qp_j(3,1)+&
                                           tij_abcd(1,2,2,k)*qp_j(1,2)+&
                                           tij_abcd(2,2,2,k)*qp_j(2,2)+&
                                           tij_abcd(3,2,2,k)*qp_j(3,2)+&
                                           tij_abcd(1,3,2,k)*qp_j(1,3)+&
                                           tij_abcd(2,3,2,k)*qp_j(2,3)+&
                                           tij_abcd(3,3,2,k)*qp_j(3,3))+&
                                  dp_i(3)*(tij_abcd(1,1,3,k)*qp_j(1,1)+&
                                           tij_abcd(2,1,3,k)*qp_j(2,1)+&
                                           tij_abcd(3,1,3,k)*qp_j(3,1)+&
                                           tij_abcd(1,2,3,k)*qp_j(1,2)+&
                                           tij_abcd(2,2,3,k)*qp_j(2,2)+&
                                           tij_abcd(3,2,3,k)*qp_j(3,2)+&
                                           tij_abcd(1,3,3,k)*qp_j(1,3)+&
                                           tij_abcd(2,3,3,k)*qp_j(2,3)+&
                                           tij_abcd(3,3,3,k)*qp_j(3,3))

                         ! Dipole j (locally A) - Quadrupole i (locally B)
                         tmp_ji = dp_j(1)*(tij_abcd(1,1,1,k)*qp_i(1,1)+&
                                           tij_abcd(2,1,1,k)*qp_i(2,1)+&
                                           tij_abcd(3,1,1,k)*qp_i(3,1)+&
                                           tij_abcd(1,2,1,k)*qp_i(1,2)+&
                                           tij_abcd(2,2,1,k)*qp_i(2,2)+&
                                           tij_abcd(3,2,1,k)*qp_i(3,2)+&
                                           tij_abcd(1,3,1,k)*qp_i(1,3)+&
                                           tij_abcd(2,3,1,k)*qp_i(2,3)+&
                                           tij_abcd(3,3,1,k)*qp_i(3,3))+&
                                  dp_j(2)*(tij_abcd(1,1,2,k)*qp_i(1,1)+&
                                           tij_abcd(2,1,2,k)*qp_i(2,1)+&
                                           tij_abcd(3,1,2,k)*qp_i(3,1)+&
                                           tij_abcd(1,2,2,k)*qp_i(1,2)+&
                                           tij_abcd(2,2,2,k)*qp_i(2,2)+&
                                           tij_abcd(3,2,2,k)*qp_i(3,2)+&
                                           tij_abcd(1,3,2,k)*qp_i(1,3)+&
                                           tij_abcd(2,3,2,k)*qp_i(2,3)+&
                                           tij_abcd(3,3,2,k)*qp_i(3,3))+&
                                  dp_j(3)*(tij_abcd(1,1,3,k)*qp_i(1,1)+&
                                           tij_abcd(2,1,3,k)*qp_i(2,1)+&
                                           tij_abcd(3,1,3,k)*qp_i(3,1)+&
                                           tij_abcd(1,2,3,k)*qp_i(1,2)+&
                                           tij_abcd(2,2,3,k)*qp_i(2,2)+&
                                           tij_abcd(3,2,3,k)*qp_i(3,2)+&
                                           tij_abcd(1,3,3,k)*qp_i(1,3)+&
                                           tij_abcd(2,3,3,k)*qp_i(2,3)+&
                                           tij_abcd(3,3,3,k)*qp_i(3,3))

                         fr(k) = fr(k) - fac * (tmp_ij - tmp_ji)
                      END DO
                   END IF
                END IF
                IF (task(3,1)) THEN
                   ! Quadrupole - Charge
                   fac = 1.0_dp/3.0_dp

                   ! Quadrupole j (locally A) - Charge j (locally B)
                   tmp_ij = ch_i * (tij_ab(1,1)*qp_j(1,1)+&
                                    tij_ab(2,1)*qp_j(2,1)+&
                                    tij_ab(3,1)*qp_j(3,1)+&
                                    tij_ab(1,2)*qp_j(1,2)+&
                                    tij_ab(2,2)*qp_j(2,2)+&
                                    tij_ab(3,2)*qp_j(3,2)+&
                                    tij_ab(1,3)*qp_j(1,3)+&
                                    tij_ab(2,3)*qp_j(2,3)+&
                                    tij_ab(3,3)*qp_j(3,3))

                   ! Quadrupole i (locally B) - Charge j (locally A)
                   tmp_ji = ch_j * (tij_ab(1,1)*qp_i(1,1)+&
                                    tij_ab(2,1)*qp_i(2,1)+&
                                    tij_ab(3,1)*qp_i(3,1)+&
                                    tij_ab(1,2)*qp_i(1,2)+&
                                    tij_ab(2,2)*qp_i(2,2)+&
                                    tij_ab(3,2)*qp_i(3,2)+&
                                    tij_ab(1,3)*qp_i(1,3)+&
                                    tij_ab(2,3)*qp_i(2,3)+&
                                    tij_ab(3,3)*qp_i(3,3))

                   eloc = eloc + fac*(tmp_ij+tmp_ji)
                   IF (do_forces.OR.do_stress) THEN
                      DO k = 1, 3
                         ! Quadrupole j (locally A) - Charge i (locally B)
                         tmp_ij = ch_i * (tij_abc(1,1,k)*qp_j(1,1)+&
                                          tij_abc(2,1,k)*qp_j(2,1)+&
                                          tij_abc(3,1,k)*qp_j(3,1)+&
                                          tij_abc(1,2,k)*qp_j(1,2)+&
                                          tij_abc(2,2,k)*qp_j(2,2)+&
                                          tij_abc(3,2,k)*qp_j(3,2)+&
                                          tij_abc(1,3,k)*qp_j(1,3)+&
                                          tij_abc(2,3,k)*qp_j(2,3)+&
                                          tij_abc(3,3,k)*qp_j(3,3))

                         ! Quadrupole i (locally B) - Charge j (locally A)
                         tmp_ji = ch_j * (tij_abc(1,1,k)*qp_i(1,1)+&
                                          tij_abc(2,1,k)*qp_i(2,1)+&
                                          tij_abc(3,1,k)*qp_i(3,1)+&
                                          tij_abc(1,2,k)*qp_i(1,2)+&
                                          tij_abc(2,2,k)*qp_i(2,2)+&
                                          tij_abc(3,2,k)*qp_i(3,2)+&
                                          tij_abc(1,3,k)*qp_i(1,3)+&
                                          tij_abc(2,3,k)*qp_i(2,3)+&
                                          tij_abc(3,3,k)*qp_i(3,3))

                         fr(k) = fr(k) - fac *(tmp_ij + tmp_ji)
                      END DO
                   END IF
                END IF
#if 0
                energy = energy + eloc
#endif
#if 0
                IF (do_forces) THEN
                   forces(1,atom_a) = forces(1,atom_a) - fr(1)
                   forces(2,atom_a) = forces(2,atom_a) - fr(2)
                   forces(3,atom_a) = forces(3,atom_a) - fr(3)
                   forces(1,atom_b) = forces(1,atom_b) + fr(1)
                   forces(2,atom_b) = forces(2,atom_b) + fr(2)
                   forces(3,atom_b) = forces(3,atom_b) + fr(3)
                END IF
#endif
                ! Electric fields
                IF (do_efield) THEN
                   ! Potential
                   IF (do_efield0) THEN
                      efield0(  atom_a) = efield0(  atom_a) + ef0_j

                      efield0(  atom_b) = efield0(  atom_b) + ef0_i
                   END IF
                   ! Electric field
                   IF (do_efield1) THEN
                      efield1(1,atom_a) = efield1(1,atom_a) + ef1_j(1)
                      efield1(2,atom_a) = efield1(2,atom_a) + ef1_j(2)
                      efield1(3,atom_a) = efield1(3,atom_a) + ef1_j(3)

                      efield1(1,atom_b) = efield1(1,atom_b) + ef1_i(1)
                      efield1(2,atom_b) = efield1(2,atom_b) + ef1_i(2)
                      efield1(3,atom_b) = efield1(3,atom_b) + ef1_i(3)
                   END IF
                   ! Electric field gradient
                   IF (do_efield2) THEN
                      efield2(1,atom_a) = efield2(1,atom_a) + ef2_j(1,1)
                      efield2(2,atom_a) = efield2(2,atom_a) + ef2_j(1,2)
                      efield2(3,atom_a) = efield2(3,atom_a) + ef2_j(1,3)
                      efield2(4,atom_a) = efield2(4,atom_a) + ef2_j(2,1)
                      efield2(5,atom_a) = efield2(5,atom_a) + ef2_j(2,2)
                      efield2(6,atom_a) = efield2(6,atom_a) + ef2_j(2,3)
                      efield2(7,atom_a) = efield2(7,atom_a) + ef2_j(3,1)
                      efield2(8,atom_a) = efield2(8,atom_a) + ef2_j(3,2)
                      efield2(9,atom_a) = efield2(9,atom_a) + ef2_j(3,3)

                      efield2(1,atom_b) = efield2(1,atom_b) + ef2_i(1,1)
                      efield2(2,atom_b) = efield2(2,atom_b) + ef2_i(1,2)
                      efield2(3,atom_b) = efield2(3,atom_b) + ef2_i(1,3)
                      efield2(4,atom_b) = efield2(4,atom_b) + ef2_i(2,1)
                      efield2(5,atom_b) = efield2(5,atom_b) + ef2_i(2,2)
                      efield2(6,atom_b) = efield2(6,atom_b) + ef2_i(2,3)
                      efield2(7,atom_b) = efield2(7,atom_b) + ef2_i(3,1)
                      efield2(8,atom_b) = efield2(8,atom_b) + ef2_i(3,2)
                      efield2(9,atom_b) = efield2(9,atom_b) + ef2_i(3,3)
                   END IF
                END IF
                IF (do_stress) THEN
                   ptens11 = ptens11 + rab(1) * fr(1)
                   ptens21 = ptens21 + rab(2) * fr(1)
                   ptens31 = ptens31 + rab(3) * fr(1)
                   ptens12 = ptens12 + rab(1) * fr(2)
                   ptens22 = ptens22 + rab(2) * fr(2)
                   ptens32 = ptens32 + rab(3) * fr(2)
                   ptens13 = ptens13 + rab(1) * fr(3)
                   ptens23 = ptens23 + rab(2) * fr(3)
                   ptens33 = ptens33 + rab(3) * fr(3)
                END IF
