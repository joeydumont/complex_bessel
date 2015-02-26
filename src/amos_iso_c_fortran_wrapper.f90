! -------------------------------------------------------------------
! - Author: 		Joey Dumont <joey.dumont@gmail.com>         -
!                       Denis Gagnon <gagnon88@gmail.com>           -
! - Date created:	2013-11-18                                  -
! - Date modded:	2015-02-26                                  -
! - Description:        ISO C Binding wrapper for the FORTRAN       -
!                       subroutines contained in D. E. Amos' Bessel -
!                       functions library.                          -
! -------------------------------------------------------------------

! Bessel function of the first kind.
subroutine zbesj_wrap(zr, zi, order, kode, N, cyr, cyi, nz, ierr) bind(C)

                ! State ISO_C_BINDING and no implicit variable types.
                use iso_c_binding
                implicit none

                ! State intent of FORTRAN variables.
                real(c_double), value, intent(in)     :: zr, zi, order
                integer(c_int), value, intent(in)     :: kode, N
                real(c_double), intent(out)           :: cyr, cyi
                integer(c_int), intent(out)           :: nz, ierr

                ! Interface to original FORTRAN subroutine.
                interface
                        SUBROUTINE ZBESJ(ZR, ZI, ORDER, KODE, N, CYR, CYI, NZ, IERR)
                                INTEGER KODE, N, NZ, IERR
                                DOUBLE PRECISION ZR, ZI, ORDER, CYR, CYI
                        END SUBROUTINE ZBESJ
                end interface

                call ZBESJ(zr, zi, order, kode, N, cyr, cyi, nz, ierr)
end subroutine zbesj_wrap

! Bessel function of the second kind.
subroutine zbesy_wrap(zr, zi, order, kode, N, cyr, cyi, nz, cwrkr, cwrki, ierr) bind(C)

                ! State ISO_C_BINDING and no implicit variable types.
                use iso_c_binding
                implicit none

                ! State intent of FORTRAN variables.
                real(c_double), value, intent(in)      :: zr, zi, order
                integer(c_int), value, intent(in)      :: kode, N
                real(c_double), intent(out)            :: cyr, cyi
                real(c_double), intent(in)             :: cwrkr, cwrki
                integer(c_int), intent(out)            :: nz, ierr

                ! Interface to original FORTRAN subroutine.
                interface
                        SUBROUTINE ZBESY(ZR, ZI, ORDER, KODE, N, CYR, CYI, NZ, CWRKR, CWRKI, IERR)
                                INTEGER KODE, N, NZ, IERR
                                DOUBLE PRECISION ZR, ZI, ORDER, CYR, CYI, CWRKR, CWRKI
                        END SUBROUTINE ZBESY
                end interface

                call ZBESY(zr, zi, order, kode, N, cyr, cyi, nz, cwrkr, cwrki, ierr)
end subroutine zbesy_wrap

! Modified Bessel function of the first kind.
subroutine zbesi_wrap(zr, zi, order, kode, N, cyr, cyi, nz, ierr) bind(C)

                ! State ISO_C_BINDING and no implicit variable types.
                use iso_c_binding
                implicit none

                ! State intent of FORTRAN variables.
                real(c_double), value, intent(in)      :: zr, zi, order
                integer(c_int), value, intent(in)       :: kode, N
                real(c_double), intent(out)            :: cyr, cyi
                integer(c_int), intent(out)             :: nz, ierr

                ! Interface to original FORTRAN subroutine.
                interface
                        SUBROUTINE ZBESI(ZR, ZI, ORDER, KODE, N, CYR, CYI, NZ, IERR)
                                INTEGER KODE, N, NZ, IERR
                                DOUBLE PRECISION ZR, ZI, ORDER, CYR, CYI
                        END SUBROUTINE ZBESI
                end interface

                call ZBESI(zr, zi, order, kode, N, cyr, cyi, nz, ierr)
end subroutine zbesi_wrap

! Modified Bessel function of the second kind.
subroutine zbesk_wrap(zr, zi, order, kode, N, cyr, cyi, nz, ierr) bind(C)

                ! State ISO_C_BINDING and no implicit variable types.
                use iso_c_binding
                implicit none

                ! State intent of FORTRAN variables.
                real(c_double), value, intent(in)     :: zr, zi, order
                integer(c_int), value, intent(in)     :: kode, N
                real(c_double), intent(out)           :: cyr, cyi
                integer(c_int), intent(out)           :: nz, ierr

                ! Interface to original FORTRAN subroutine.
                interface
                        SUBROUTINE ZBESK(ZR, ZI, ORDER, KODE, N, CYR, CYI, NZ, IERR)
                                INTEGER KODE, N, NZ, IERR
                                DOUBLE PRECISION ZR, ZI, ORDER, CYR, CYI
                        END SUBROUTINE ZBESK
                end interface

                call ZBESK(zr, zi, order, kode, N, cyr, cyi, nz, ierr)
end subroutine zbesk_wrap

! Hankel functions of both kinds.
subroutine zbesh_wrap(zr, zi, order, kode, M, N, cyr, cyi, nz, ierr) bind(C)

                ! State ISO_C_BINDING and no implicit variable types.
                use iso_c_binding
                implicit none

                ! State intent of FORTRAN variables.
                real(c_double), value, intent(in)      :: zr, zi, order
                integer(c_int), value, intent(in)      :: kode, M, N
                real(c_double), intent(out)            :: cyr, cyi
                integer(c_int), intent(out)            :: nz, ierr

                ! Interface to original FORTRAN subroutine.
                interface
                        SUBROUTINE ZBESH(ZR, ZI, ORDER, KODE, M, N, CYR, CYI, NZ, IERR)
                                INTEGER KODE, M, N, NZ, IERR
                                DOUBLE PRECISION ZR, ZI, ORDER, CYR, CYI
                        END SUBROUTINE ZBESH
                end interface

                call ZBESH(zr, zi, order, kode, M, N, cyr, cyi, nz, ierr)
end subroutine zbesh_wrap

! Airy function of the first kind.
subroutine zairy_wrap(zr, zi, id, kode, air, aii, nz, ierr) bind(C)

                ! State ISO_C_BINDING and no implicit variable types.
                use iso_c_binding
                implicit none

                ! State intent of FORTRAN variables.
                real(c_double), value, intent(in)               :: zr, zi
                integer(c_int), value, intent(in)               :: id, kode
                real(c_double), intent(out)                     :: air, aii
                integer(c_int), intent(out)                     :: nz, ierr

                ! Interface to original FORTRAN subroutine.
                interface
                        SUBROUTINE ZAIRY(ZR, ZI, ID, KODE, AIR, AII, NZ, IERR)
                                INTEGER ID, KODE, NZ, IERR
                                DOUBLE PRECISION ZR, ZI, AIR, AII
                        END SUBROUTINE ZAIRY
                end interface

                call ZAIRY(zr, zi, id, kode, air, aii, nz, ierr)
end subroutine zairy_wrap

! Airy function of the second kind.
subroutine zbiry_wrap(zr, zi, id, kode, bir, bii, ierr) bind(C)

                ! State ISO_C_BINDING and no implicit variable types.
                use iso_c_binding
                implicit none

                ! State intent of FORTRAN variables.
                real(c_double), value, intent(in)               :: zr, zi
                integer(c_int), value, intent(in)               :: id, kode
                real(c_double), intent(out)                     :: bir, bii
                integer(c_int), intent(out)                     :: ierr

                ! Interface to original FORTRAN subroutine.
                interface
                        SUBROUTINE BAIRY(ZR, ZI, ID, KODE, BIR, BII, IERR)
                                INTEGER ID, KODE, IERR
                                DOUBLE PRECISION ZR, ZI, BIR, BII
                        END SUBROUTINE BAIRY
                end interface

                call ZBIRY(zr, zi, id, kode, bir, bii, ierr)
end subroutine zbiry_wrap
