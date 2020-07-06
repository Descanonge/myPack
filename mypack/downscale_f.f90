! Various functions

module downscale
  implicit none
  public

contains

  subroutine downscale_average(x_i, y_i, nx_i, ny_i, &
       x_o, y_o, nx_o, ny_o, in, out)
    ! Down sample by average.
    ! nx_i > nx_o
    ! x_* and y_* increasing
    implicit none
    integer, intent(in) :: nx_i, ny_i, nx_o, ny_o
    double precision, intent(in) :: x_i(nx_i), y_i(ny_i)
    double precision, intent(in) :: x_o(nx_o), y_o(ny_o)
    double precision, intent(in) :: in(ny_i, nx_i)
    double precision, intent(out) :: out(ny_o, nx_o)

    integer :: map_x(nx_i), map_y(ny_i)
    integer :: n(ny_o, nx_o)
    integer :: ix, iy
    integer :: jx, jy

    call down_coord(x_i, nx_i, x_o, nx_o, map_x)
    call down_coord(y_i, ny_i, y_o, ny_o, map_y)

    out = 0.
    n = 0
    do ix=1, nx_i
       jx = map_x(ix)
       do iy=1, ny_i
          jy = map_y(iy)
          out(jy, jx) = out(jy, jx) + in(iy, ix)
          n(jy, jx) = n(jy, jx) + 1
       end do
    end do
    out = out / n

  end subroutine downscale_average

  subroutine downscale_average_stack(x_i, y_i, nx_i, ny_i, nt, &
       x_o, y_o, nx_o, ny_o, in, out)
    implicit none
    integer, intent(in) :: nx_i, ny_i, nt, nx_o, ny_o
    double precision, intent(in) :: x_i(nx_i), y_i(ny_i)
    double precision, intent(in) :: x_o(nx_o), y_o(ny_o)
    double precision, intent(in) :: in(nt, ny_i, nx_i)
    double precision, intent(out) :: out(nt, ny_o, nx_o)

    integer :: map_x(nx_i), map_y(ny_i)
    integer :: n(nt, ny_o, nx_o)
    integer :: ix, iy, it
    integer :: jx, jy

    call down_coord(x_i, nx_i, x_o, nx_o, map_x)
    call down_coord(y_i, ny_i, y_o, ny_o, map_y)

    out = 0.
    n = 0
    do it=1, nt
        do ix=1, nx_i
            jx = map_x(ix)
            do iy=1, ny_i
                jy = map_y(iy)
                out(it, jy, jx) = out(it, jy, jx) + in(it, iy, ix)
                n(it, jy, jx) = n(it, jy, jx) + 1
            end do
        end do
    end do
    out = out / n

  end subroutine downscale_average_stack

  subroutine down_coord(x_i, n_i, x_o, n_o, map)
    implicit none
    integer, intent(in) :: n_i, n_o
    double precision, intent(in) :: x_i(n_i), x_o(n_o)
    integer, intent(out) :: map(n_i)

    integer :: i, j
    double precision :: step_i, bin_end_i

    step_i = 0
    do i=1, n_i-1
       step_i = ((i-1)*step_i + x_i(i+1)-x_i(i)) / i
    enddo

    do i=1, n_i
       j = 1
       do while (j < n_o)
          if (x_i(i) < x_o(j+1)) then
             exit
          endif
          j = j + 1
       enddo
       ! If last bin, no way to know its size,
       ! unless we assume bin length is regular
       if (i == n_i) then
          bin_end_i = x_i(i) + step_i
       else
          bin_end_i = x_i(i+1)
       endif
       if (j < n_o .and. bin_end_i-x_o(j+1) > x_o(j+1)-x_i(i) ) then
          j = j + 1
       endif
       map(i) = j
    enddo

    ! ! check
    ! do i=1, n_i
    !    write(*, *) i, x_i(i), map(i), x_o(map(i))
    ! enddo

  end subroutine down_coord
 
end module downscale
