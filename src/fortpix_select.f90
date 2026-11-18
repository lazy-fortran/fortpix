module fortpix_select
    !! Selection masks: coverage from 0 to 1 per pixel, antialiased at
    !! geometric edges. Coordinates are continuous, pixel (i, j) covering
    !! [i-1, i] x [j-1, j].
    use iso_fortran_env, only: int8
    use fortpix_base, only: u8, to_byte_real
    implicit none
    private
    public :: pix_mask_rect, pix_mask_ellipse, pix_mask_polygon, pix_magic_wand, &
        pix_mask_combine, pix_mask_bounds, pix_blend_masked, &
        MASK_REPLACE, MASK_ADD, MASK_SUBTRACT, MASK_INTERSECT

    integer, parameter :: MASK_REPLACE = 0, MASK_ADD = 1, MASK_SUBTRACT = 2, MASK_INTERSECT = 3
    integer, parameter :: SUB = 4

contains

    subroutine pix_mask_rect(w, h, x0, y0, x1, y1, mask)
        integer, intent(in) :: w, h
        real, intent(in) :: x0, y0, x1, y1
        real, allocatable, intent(out) :: mask(:, :)
        real :: a, b, c, d
        integer :: x, y
        a = min(x0, x1)
        b = min(y0, y1)
        c = max(x0, x1)
        d = max(y0, y1)
        allocate(mask(w, h))
        do y = 1, h
            do x = 1, w
                mask(x, y) = overlap(real(x - 1), real(x), a, c) * &
                    overlap(real(y - 1), real(y), b, d)
            end do
        end do
    end subroutine pix_mask_rect

    pure real function overlap(p0, p1, q0, q1)
        real, intent(in) :: p0, p1, q0, q1
        overlap = max(0.0, min(p1, q1) - max(p0, q0))
    end function overlap

    subroutine pix_mask_ellipse(w, h, x0, y0, x1, y1, mask)
        !! Ellipse inscribed in the box, SUB x SUB samples per edge pixel.
        integer, intent(in) :: w, h
        real, intent(in) :: x0, y0, x1, y1
        real, allocatable, intent(out) :: mask(:, :)
        real :: cx, cy, rx, ry, u, v, d
        integer :: x, y, i, j, hit
        cx = 0.5 * (x0 + x1)
        cy = 0.5 * (y0 + y1)
        rx = max(1e-3, 0.5 * abs(x1 - x0))
        ry = max(1e-3, 0.5 * abs(y1 - y0))
        allocate(mask(w, h))
        do y = 1, h
            do x = 1, w
                d = ((x - 0.5 - cx) / rx)**2 + ((y - 0.5 - cy) / ry)**2
                ! Far from the edge: whole or empty without sampling.
                if (sqrt(d) < 1.0 - 1.5 / min(rx, ry)) then
                    mask(x, y) = 1.0
                else if (sqrt(d) > 1.0 + 1.5 / min(rx, ry)) then
                    mask(x, y) = 0.0
                else
                    hit = 0
                    do j = 1, SUB
                        do i = 1, SUB
                            u = (x - 1 + (i - 0.5) / SUB - cx) / rx
                            v = (y - 1 + (j - 0.5) / SUB - cy) / ry
                            if (u * u + v * v <= 1.0) hit = hit + 1
                        end do
                    end do
                    mask(x, y) = real(hit) / SUB**2
                end if
            end do
        end do
    end subroutine pix_mask_ellipse

    subroutine pix_mask_polygon(w, h, px, py, mask)
        !! Even-odd polygon (a lasso). Each pixel row is sampled on SUB
        !! scanlines with exact horizontal coverage of the spans.
        integer, intent(in) :: w, h
        real, intent(in) :: px(:), py(:)
        real, allocatable, intent(out) :: mask(:, :)
        real, allocatable :: xs(:)
        real :: sy, x0, x1, ya, yb, t
        integer :: n, y, s, k, m, i, j, a, b
        n = size(px)
        allocate(mask(w, h), xs(n))
        mask = 0.0
        if (n < 3) return
        do y = 1, h
            do s = 1, SUB
                sy = y - 1 + (s - 0.5) / SUB
                m = 0
                do k = 1, n
                    j = modulo(k, n) + 1
                    ya = py(k)
                    yb = py(j)
                    if ((ya <= sy .and. yb > sy) .or. (yb <= sy .and. ya > sy)) then
                        m = m + 1
                        xs(m) = px(k) + (sy - ya) / (yb - ya) * (px(j) - px(k))
                    end if
                end do
                ! Insertion sort: crossings per scanline are few.
                do i = 2, m
                    t = xs(i)
                    j = i - 1
                    do while (j >= 1)
                        if (xs(j) <= t) exit
                        xs(j + 1) = xs(j)
                        j = j - 1
                    end do
                    xs(j + 1) = t
                end do
                do i = 1, m - 1, 2
                    x0 = max(0.0, xs(i))
                    x1 = min(real(w), xs(i + 1))
                    if (x1 <= x0) cycle
                    a = int(x0) + 1
                    b = min(w, int(x1) + 1)
                    do k = a, b
                        mask(k, y) = mask(k, y) + overlap(real(k - 1), real(k), x0, x1) / SUB
                    end do
                end do
            end do
        end do
        mask = min(1.0, mask)
    end subroutine pix_mask_polygon

    subroutine pix_magic_wand(src, x, y, tolerance, contiguous, mask)
        !! Pixels whose largest channel difference from the seed pixel (x, y)
        !! is at most `tolerance` (0..255); with `contiguous` only those
        !! connected to the seed through 4-neighbours.
        integer(int8), intent(in) :: src(:, :, :)
        integer, intent(in) :: x, y, tolerance
        logical, intent(in) :: contiguous
        real, allocatable, intent(out) :: mask(:, :)
        logical, allocatable :: near(:, :)
        integer, allocatable :: stack(:, :)
        integer :: w, h, seed(3), i, j, top, l, r, k
        w = size(src, 2)
        h = size(src, 3)
        allocate(mask(w, h), near(w, h))
        mask = 0.0
        if (x < 1 .or. y < 1 .or. x > w .or. y > h) return
        seed = u8(src(1:3, x, y))
        do j = 1, h
            do i = 1, w
                near(i, j) = maxval(abs(u8(src(1:3, i, j)) - seed)) <= tolerance
            end do
        end do
        if (.not. contiguous) then
            where (near) mask = 1.0
            return
        end if
        ! Scanline flood fill: each stack entry is the start of a row span.
        allocate(stack(2, 64))
        top = 1
        stack(:, 1) = [x, y]
        do while (top > 0)
            i = stack(1, top)
            j = stack(2, top)
            top = top - 1
            if (.not. near(i, j) .or. mask(i, j) > 0.0) cycle
            l = i
            do while (l > 1)
                if (.not. near(l - 1, j) .or. mask(l - 1, j) > 0.0) exit
                l = l - 1
            end do
            r = i
            do while (r < w)
                if (.not. near(r + 1, j) .or. mask(r + 1, j) > 0.0) exit
                r = r + 1
            end do
            mask(l:r, j) = 1.0
            do k = l, r
                if (j > 1) call push(k, j - 1)
                if (j < h) call push(k, j + 1)
            end do
        end do
    contains
        subroutine push(a, b)
            integer, intent(in) :: a, b
            integer, allocatable :: grown(:, :)
            if (.not. near(a, b) .or. mask(a, b) > 0.0) return
            if (top == size(stack, 2)) then
                allocate(grown(2, 2 * top))
                grown(:, 1:top) = stack
                call move_alloc(grown, stack)
            end if
            top = top + 1
            stack(:, top) = [a, b]
        end subroutine push
    end subroutine pix_magic_wand

    subroutine pix_mask_combine(mask, shape, mode)
        !! Updates `mask` with `shape` (same size) by MASK_* mode.
        real, allocatable, intent(inout) :: mask(:, :)
        real, intent(in) :: shape(:, :)
        integer, intent(in) :: mode
        if (.not. allocated(mask) .or. mode == MASK_REPLACE) then
            mask = shape
            return
        end if
        select case (mode)
        case (MASK_ADD)
            mask = max(mask, shape)
        case (MASK_SUBTRACT)
            mask = min(mask, 1.0 - shape)
        case (MASK_INTERSECT)
            mask = min(mask, shape)
        end select
    end subroutine pix_mask_combine

    subroutine pix_mask_bounds(mask, x0, y0, x1, y1)
        !! Smallest pixel box (1-based, inclusive) with nonzero coverage;
        !! x1 < x0 for an empty mask.
        real, intent(in) :: mask(:, :)
        integer, intent(out) :: x0, y0, x1, y1
        integer :: i, j
        x0 = size(mask, 1) + 1
        y0 = size(mask, 2) + 1
        x1 = 0
        y1 = 0
        do j = 1, size(mask, 2)
            do i = 1, size(mask, 1)
                if (mask(i, j) > 0.0) then
                    x0 = min(x0, i)
                    x1 = max(x1, i)
                    y0 = min(y0, j)
                    y1 = max(y1, j)
                end if
            end do
        end do
    end subroutine pix_mask_bounds

    subroutine pix_blend_masked(base, edited, mask, dst)
        !! edited where the mask covers, base elsewhere, mixed at partial
        !! coverage: confines an adjustment to a selection.
        integer(int8), intent(in) :: base(:, :, :), edited(:, :, :)
        real, intent(in) :: mask(:, :)
        integer(int8), allocatable, intent(out) :: dst(:, :, :)
        integer :: x, y
        allocate(dst, mold=base)
        do y = 1, size(base, 3)
            do x = 1, size(base, 2)
                dst(1:3, x, y) = to_byte_real(real(u8(base(1:3, x, y))) * (1.0 - mask(x, y)) + &
                    real(u8(edited(1:3, x, y))) * mask(x, y))
            end do
        end do
    end subroutine pix_blend_masked
end module fortpix_select
