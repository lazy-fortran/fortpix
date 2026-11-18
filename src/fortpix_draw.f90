module fortpix_draw
    !! Rasterizes a variable-width ink stroke into an RGB image.
    use iso_fortran_env, only: int8
    use fortpix_base, only: u8, to_byte_real
    implicit none
    private
    public :: pix_stroke

contains

    subroutine pix_stroke(img, x, y, r, rgb)
        !! Stroke through centreline samples (x, y) with radius r, all in pixel
        !! units with pixel (i, j) covering [i-1, i] x [j-1, j]. The coverage
        !! is the union of round segments whose radius varies linearly, with
        !! a one-pixel analytic edge, and the colour is applied once: a pen
        !! resting in one place does not darken its deposit.
        integer(int8), intent(inout) :: img(:, :, :)
        real, intent(in) :: x(:), y(:), r(:)
        integer, intent(in) :: rgb(3)
        real, allocatable :: cover(:, :)
        integer :: n, k, i, j, i0, i1, j0, j1, b0, b1, c0, c1
        real :: ax, ay, bx, by, ra, rb, dx, dy, len2, t, px, py, d
        n = size(x)
        if (n == 0) return
        b0 = size(img, 2) + 1
        c0 = size(img, 3) + 1
        b1 = 0
        c1 = 0
        allocate(cover(size(img, 2), size(img, 3)))
        cover = 0.0
        do k = 1, max(1, n - 1)
            ax = x(k)
            ay = y(k)
            ra = max(0.5, r(k))
            bx = x(min(n, k + 1))
            by = y(min(n, k + 1))
            rb = max(0.5, r(min(n, k + 1)))
            i0 = max(1, floor(min(ax - ra, bx - rb)) - 1)
            i1 = min(size(img, 2), ceiling(max(ax + ra, bx + rb)) + 1)
            j0 = max(1, floor(min(ay - ra, by - rb)) - 1)
            j1 = min(size(img, 3), ceiling(max(ay + ra, by + rb)) + 1)
            dx = bx - ax
            dy = by - ay
            len2 = dx * dx + dy * dy
            do j = j0, j1
                py = j - 0.5
                do i = i0, i1
                    px = i - 0.5
                    if (len2 > 0.0) then
                        t = max(0.0, min(1.0, ((px - ax) * dx + (py - ay) * dy) / len2))
                    else
                        t = 0.0
                    end if
                    d = sqrt((px - ax - t * dx)**2 + (py - ay - t * dy)**2) - (ra + t * (rb - ra))
                    cover(i, j) = max(cover(i, j), max(0.0, min(1.0, 0.5 - d)))
                end do
            end do
            b0 = min(b0, i0)
            b1 = max(b1, i1)
            c0 = min(c0, j0)
            c1 = max(c1, j1)
        end do
        do j = c0, c1
            do i = b0, b1
                if (cover(i, j) <= 0.0) cycle
                img(1:3, i, j) = to_byte_real(real(u8(img(1:3, i, j))) * (1.0 - cover(i, j)) + &
                    real(rgb) * cover(i, j))
            end do
        end do
    end subroutine pix_stroke
end module fortpix_draw
