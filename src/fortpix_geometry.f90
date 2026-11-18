module fortpix_geometry
    !! Crop, quarter turns, flips, free rotation and resampling.
    use iso_fortran_env, only: int8
    use fortpix_base, only: u8, to_byte_real
    implicit none
    private
    public :: pix_crop, pix_rotate90, pix_flip, pix_rotate, pix_straighten_scale, pix_resize

    real, parameter :: PI = 3.14159265358979

contains

    subroutine pix_crop(src, x0, y0, w, h, dst)
        !! Pixels x0..x0+w-1, y0..y0+h-1 (1-based), clipped to the image.
        integer(int8), intent(in) :: src(:, :, :)
        integer, intent(in) :: x0, y0, w, h
        integer(int8), allocatable, intent(out) :: dst(:, :, :)
        integer :: a, b, c, d
        a = max(1, x0)
        b = max(1, y0)
        c = min(size(src, 2), x0 + w - 1)
        d = min(size(src, 3), y0 + h - 1)
        allocate(dst(size(src, 1), max(0, c - a + 1), max(0, d - b + 1)))
        if (c >= a .and. d >= b) dst = src(:, a:c, b:d)
    end subroutine pix_crop

    subroutine pix_rotate90(src, turns, dst)
        !! Clockwise quarter turns (any integer; negative turns counter-clockwise).
        integer(int8), intent(in) :: src(:, :, :)
        integer, intent(in) :: turns
        integer(int8), allocatable, intent(out) :: dst(:, :, :)
        integer :: w, h, x, y, q
        w = size(src, 2)
        h = size(src, 3)
        q = modulo(turns, 4)
        select case (q)
        case (0)
            dst = src
        case (1)
            allocate(dst(size(src, 1), h, w))
            do y = 1, h
                do x = 1, w
                    dst(:, h - y + 1, x) = src(:, x, y)
                end do
            end do
        case (2)
            allocate(dst(size(src, 1), w, h))
            do y = 1, h
                do x = 1, w
                    dst(:, w - x + 1, h - y + 1) = src(:, x, y)
                end do
            end do
        case default
            allocate(dst(size(src, 1), h, w))
            do y = 1, h
                do x = 1, w
                    dst(:, y, w - x + 1) = src(:, x, y)
                end do
            end do
        end select
    end subroutine pix_rotate90

    subroutine pix_flip(src, horizontal, dst)
        integer(int8), intent(in) :: src(:, :, :)
        logical, intent(in) :: horizontal
        integer(int8), allocatable, intent(out) :: dst(:, :, :)
        integer :: w, h
        w = size(src, 2)
        h = size(src, 3)
        if (horizontal) then
            dst = src(:, w:1:-1, :)
        else
            dst = src(:, :, h:1:-1)
        end if
    end subroutine pix_flip

    pure real function pix_straighten_scale(w, h, degrees) result(s)
        !! Largest scale of a centred w x h rectangle that stays inside the
        !! image after rotating it by `degrees`.
        integer, intent(in) :: w, h
        real, intent(in) :: degrees
        real :: c, sn
        c = abs(cos(degrees * PI / 180.0))
        sn = abs(sin(degrees * PI / 180.0))
        s = min(real(w) / (w * c + h * sn), real(h) / (w * sn + h * c))
    end function pix_straighten_scale

    subroutine pix_rotate(src, degrees, crop, fill, dst)
        !! Rotation by `degrees` clockwise about the centre with bicubic
        !! sampling. With `crop` the result keeps the image's aspect ratio
        !! and contains no corners from outside; otherwise it is the bounding
        !! box filled with `fill` (bytes) where no source pixel lands.
        integer(int8), intent(in) :: src(:, :, :)
        real, intent(in) :: degrees
        logical, intent(in) :: crop
        integer, intent(in) :: fill(3)
        integer(int8), allocatable, intent(out) :: dst(:, :, :)
        integer :: w, h, nw, nh, x, y, k
        real :: c, sn, cx, cy, ox, oy, u, v, s, val(3)
        logical :: inside

        w = size(src, 2)
        h = size(src, 3)
        c = cos(degrees * PI / 180.0)
        sn = sin(degrees * PI / 180.0)
        if (crop) then
            s = pix_straighten_scale(w, h, degrees)
            nw = max(1, int(w * s))
            nh = max(1, int(h * s))
        else
            nw = ceiling(w * abs(c) + h * abs(sn) - 1e-3)
            nh = ceiling(w * abs(sn) + h * abs(c) - 1e-3)
        end if
        allocate(dst(3, nw, nh))
        cx = 0.5 * real(w)
        cy = 0.5 * real(h)
        ox = 0.5 * real(nw)
        oy = 0.5 * real(nh)
        !$omp parallel do private(x, u, v, val, inside, k)
        do y = 1, nh
            do x = 1, nw
                ! Pixel centres, inverse rotation back into the source.
                u = c * (x - 0.5 - ox) + sn * (y - 0.5 - oy) + cx
                v = -sn * (x - 0.5 - ox) + c * (y - 0.5 - oy) + cy
                inside = u >= 0.0 .and. v >= 0.0 .and. u <= real(w) .and. v <= real(h)
                if (inside) then
                    call bicubic(src, u, v, val)
                    dst(:, x, y) = to_byte_real(val)
                else
                    do k = 1, 3
                        dst(k, x, y) = to_byte_real(real(fill(k)))
                    end do
                end if
            end do
        end do
        !$omp end parallel do
    end subroutine pix_rotate

    pure real function cubic(t)
        !! Catmull-Rom kernel (a = -0.5).
        real, intent(in) :: t
        real :: x
        x = abs(t)
        if (x < 1.0) then
            cubic = (1.5 * x - 2.5) * x * x + 1.0
        else if (x < 2.0) then
            cubic = ((-0.5 * x + 2.5) * x - 4.0) * x + 2.0
        else
            cubic = 0.0
        end if
    end function cubic

    pure subroutine bicubic(src, u, v, val)
        !! Sample at continuous position (u, v), pixel centres at k - 0.5,
        !! edges clamped.
        integer(int8), intent(in) :: src(:, :, :)
        real, intent(in) :: u, v
        real, intent(out) :: val(3)
        integer :: i0, j0, i, j, ii, jj
        real :: fx, fy, wx(4), wy(4)
        fx = u - 0.5
        fy = v - 0.5
        i0 = floor(fx)
        j0 = floor(fy)
        do i = 1, 4
            wx(i) = cubic(fx - real(i0 + i - 2))
            wy(i) = cubic(fy - real(j0 + i - 2))
        end do
        val = 0.0
        do j = 1, 4
            jj = max(1, min(size(src, 3), j0 + j - 1))
            do i = 1, 4
                ii = max(1, min(size(src, 2), i0 + i - 1))
                val = val + wx(i) * wy(j) * real(u8(src(1:3, ii, jj)))
            end do
        end do
    end subroutine bicubic

    pure real function lanczos3(t)
        real, intent(in) :: t
        real :: x
        x = abs(t)
        if (x < 1e-6) then
            lanczos3 = 1.0
        else if (x < 3.0) then
            lanczos3 = 3.0 * sin(PI * x) * sin(PI * x / 3.0) / (PI * PI * x * x)
        else
            lanczos3 = 0.0
        end if
    end function lanczos3

    subroutine weights(n_in, n_out, first, count, w)
        !! Lanczos-3 taps for each output sample; the kernel widens when
        !! shrinking so that it also filters out frequencies the result
        !! cannot hold.
        integer, intent(in) :: n_in, n_out
        integer, allocatable, intent(out) :: first(:), count(:)
        real, allocatable, intent(out) :: w(:, :)
        real :: scale, support, centre, total
        integer :: k, j, lo, hi, taps
        scale = real(n_in) / real(n_out)
        support = 3.0 * max(1.0, scale)
        taps = 2 * ceiling(support) + 2
        allocate(first(n_out), count(n_out), w(taps, n_out))
        w = 0.0
        do k = 1, n_out
            centre = (real(k) - 0.5) * scale
            lo = max(1, floor(centre - support + 0.5))
            hi = min(n_in, ceiling(centre + support + 0.5))
            hi = min(hi, lo + taps - 1)
            first(k) = lo
            count(k) = hi - lo + 1
            total = 0.0
            do j = lo, hi
                w(j - lo + 1, k) = lanczos3((real(j) - 0.5 - centre) / max(1.0, scale))
                total = total + w(j - lo + 1, k)
            end do
            if (total /= 0.0) w(1:count(k), k) = w(1:count(k), k) / total
        end do
    end subroutine weights

    subroutine pix_resize(src, nw, nh, dst)
        !! Separable Lanczos-3 resampling to nw x nh.
        integer(int8), intent(in) :: src(:, :, :)
        integer, intent(in) :: nw, nh
        integer(int8), allocatable, intent(out) :: dst(:, :, :)
        integer, allocatable :: fx(:), cx(:), fy(:), cy(:)
        real, allocatable :: wx(:, :), wy(:, :), tmp(:, :, :)
        real :: acc(3)
        integer :: w, h, x, y, j
        w = size(src, 2)
        h = size(src, 3)
        allocate(dst(3, max(1, nw), max(1, nh)))
        if (w == 0 .or. h == 0) return
        call weights(w, nw, fx, cx, wx)
        call weights(h, nh, fy, cy, wy)
        allocate(tmp(3, nw, h))
        !$omp parallel do private(x, j, acc)
        do y = 1, h
            do x = 1, nw
                acc = 0.0
                do j = 1, cx(x)
                    acc = acc + wx(j, x) * real(u8(src(1:3, fx(x) + j - 1, y)))
                end do
                tmp(:, x, y) = acc
            end do
        end do
        !$omp end parallel do
        !$omp parallel do private(x, j, acc)
        do y = 1, nh
            do x = 1, nw
                acc = 0.0
                do j = 1, cy(y)
                    acc = acc + wy(j, y) * tmp(:, x, fy(y) + j - 1)
                end do
                dst(:, x, y) = to_byte_real(acc)
            end do
        end do
        !$omp end parallel do
    end subroutine pix_resize
end module fortpix_geometry
