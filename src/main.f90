! sasa
! Copyright (C) 2016 James W. Barnett <jbarnet4@tulane.edu>
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51
! Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
! The full license is located in a text file titled "LICENSE" in the root
! directory of the source.

module subs

    implicit none
    real(8), parameter :: pi = 2.0d0*acos(0.0d0)

contains


    subroutine init_seed()

        implicit none
        integer :: n, u
        integer, allocatable :: seed(:)

        call random_seed(size=n)

        allocate(seed(n))
        open(newunit=u, file="/dev/urandom", access="stream", form="unformatted", &
            action="read", status="old")
        read(u) seed(:)
        close(u)

        call random_seed(put=seed)

    end subroutine init_seed


    function gen_sphere_point(center, r)

        implicit none
        real(8) :: xi(2),  zeta(2),  zeta2, gen_sphere_point(3)
        real(8), intent(in) :: center(3), r

        zeta2 = 100.0d0

        do while (zeta2 .gt. 1.0d0)
            call random_number(xi)
            zeta = 1.0d0 - 2.0d0*xi
            zeta2 = sum(zeta**2) 
        end do

        gen_sphere_point = [ 2.0d0 * zeta(1) * sqrt(1.0d0 - zeta2), &
                             2.0d0 * zeta(2) * sqrt(1.0d0 - zeta2), &
                             1.0d0 - 2.0d0*zeta2 ]
        gen_sphere_point = gen_sphere_point * r + center

    end function gen_sphere_point


end module subs


program sasa

    use subs
    use gmxfort_trajectory
    use gmxfort_utils
    use json_module
    use omp_lib
    use, intrinsic :: iso_fortran_env

    implicit none
    real(8) :: r, area_sasa_avg, area_sasa_err, rvdw_tmp
    real(8), allocatable :: area_sasa(:), area_sasa_avg_block(:), rvdw(:)
    integer :: nsites, nrand, u, nblocks, block_size, n, nthread
    type(Trajectory) :: trj
    type(json_file) :: config
    character (len=32) :: start_time, finish_time
    character (len=2048) :: arg
    character (len=:), allocatable :: xtcfile, ndxfile, ndxgroup, outfile
    logical :: found
    real(8), parameter :: fourPi = 4.0d0 * pi
    integer(8) :: start_i, finish_i

    start_i = time8()

    if (command_argument_count() .ne. 1) then
        write(error_unit, "(a)") "ERROR: First argument should be config file."
        call abort()
    end if

    call get_command_argument(1, arg)

    call init_seed()

    call config%initialize()
    call config%load_file(filename=trim(arg))

    call config%get("input.xtcfile", xtcfile, found)
    if (.not. found) xtcfile = "traj.xtc"

    call config%get("input.ndxfile", ndxfile, found)
    if (.not. found) ndxfile = "index.ndx"

    call config%get("input.ndxgroup", ndxgroup, found)
    if (.not. found) ndxgroup = "Site"

    call config%get("config.nrand", nrand, found)
    if (.not. found) nrand = 1000

    call config%get("config.r", r, found)
    if (.not. found) r = 0.14d0

    call config%get("config.nblocks", nblocks, found)
    if (.not. found) nblocks = 5

    call config%get("output.file", outfile, found)
    if (.not. found) outfile = "sasa.txt"

    call trj%read(xtcfile, ndxfile, ndxgroup)

    n = trj%nframes
    nsites = trj%natoms()

    call config%get("config.rvdw", rvdw, found)

    if (.not. found) then
        write(error_unit, "(a)") "WARNING: No VDW radii in input file. Using 0.2 nm for each radius."
        allocate(rvdw(nsites))
        rvdw = 0.2d0
    end if

    if (size(rvdw) .ne. nsites) then
        if (size(rvdw) .eq. 1) then
            rvdw_tmp = rvdw(1)
            deallocate(rvdw)
            allocate(rvdw(nsites))
            rvdw = rvdw_tmp
            write(error_unit, "(a,f9.6,a)") "NOTE: Setting all VDW radii to ", rvdw_tmp, " nm."
        else
            write(error_unit, "(a,a,a)") "ERROR: Number of VDW radii does not equal number of atoms in ", ndxgroup, "."
            call abort()
        end if
    end if

    allocate(area_sasa(n))
    area_sasa = 0.0d0

    !$omp parallel

    !$omp master
    nthread = omp_get_num_threads()
    write(error_unit, "(a,i0,a)") "Using ", nthread, " OpenMP threads"
    !$omp end master

    mainloop : block

        integer :: i, j, k, l, naccept_sasa
        real(8) :: rand_point_sasa(3), rsasa(nsites), rsasa2(nsites)

        rsasa = r+rvdw
        rsasa2 = rsasa**2

        !$omp do 
        do l = 1, n

            if (modulo(l,100) .eq. 0) write(error_unit, "(i0)") l

            do i = 1, nsites

                naccept_sasa = 0

                do j = 1, nrand

                    rand_point_sasa = gen_sphere_point(dble(trj%x(l, i)), rsasa(i))

                    do k = 1, nsites

                        if (i .ne. k) then
                            ! Point is too close to one of the sites: reject
                            if (distance2(dble(trj%x(l, k)), rand_point_sasa, dble(trj%box(l))) .lt. rsasa2(i)) goto 100
                        end if

                    end do

                    naccept_sasa = naccept_sasa + 1

100                 continue

                end do

                area_sasa(l) = area_sasa(l) + rsasa2(i) * dble(naccept_sasa)

            end do

        end do
        !$omp end do

    end block mainloop

    !$omp end parallel

    area_sasa = area_sasa * fourPi / dble(nrand)
    area_sasa_avg = sum(area_sasa) / dble(trj%nframes)

    block_size = trj%nframes / nblocks
    allocate(area_sasa_avg_block(nblocks))

    !$omp parallel

    blockavgs : block

        integer :: i, start_frame, end_frame

        !$omp do
        do i = 1, nblocks
            start_frame = (i-1)*block_size + 1
            end_frame = i*block_size + 1
            if (i .eq. nblocks) end_frame = trj%nframes
            area_sasa_avg_block(i) = sum(area_sasa(start_frame:end_frame)) / dble(end_frame - start_frame + 1)
        end do
        !$omp end do

    end block blockavgs

    !$omp end parallel

    area_sasa_err = dsqrt( sum((area_sasa_avg_block - area_sasa_avg)**2) / dble(nblocks-1) )
    deallocate(area_sasa_avg_block)

    finish_i = time8()
    call ctime(start_i, start_time)
    call ctime(finish_i, finish_time)

    open(newunit=u, file=trim(outfile))
    write(u, "(a,a,a)") "# Program started on ", trim(start_time), "."
    write(u, "(a,a,a)") "# Program completed on ", trim(finish_time), "."
    write(u, "(a)") "# Configuration file read in:"
    call config%print_file(u)
    write(u, "(a)") "# Solvent accessible surface area (nm²)"
    write(u, "(f12.6,a,f12.6)") area_sasa_avg, " ± ", area_sasa_err
    close(u)

    call config%destroy()

end program sasa
