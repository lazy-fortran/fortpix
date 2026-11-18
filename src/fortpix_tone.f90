module fortpix_tone
    !! Histograms, levels and colour adjustments of RGB images.
    use iso_fortran_env, only: int8, int64
    use fortpix_base, only: u8, to_byte_real, srgb_to_linear, linear_to_srgb, luma
    implicit none
    private
    public :: pix_histogram, pix_levels_t, pix_levels, pix_auto_levels, pix_adjust_t, &
        pix_adjust, pix_neutral_gains, pix_is_identity

    type :: pix_levels_t
        !! Input black and white points (0..255) and gamma per channel R, G, B,
        !! applied after the master values of the same meaning.
        real :: black(0:3) = 0.0, white(0:3) = 255.0, gamma(0:3) = 1.0
    end type pix_levels_t

    type :: pix_adjust_t
        !! brightness, contrast, saturation, temperature and tint on -1..1;
        !! exposure in stops; gains multiply linear R, G, B (white balance).
        real :: exposure = 0.0, brightness = 0.0, contrast = 0.0
        real :: saturation = 0.0, temperature = 0.0, tint = 0.0
        real :: gains(3) = 1.0
    end type pix_adjust_t

    interface pix_is_identity
        module procedure levels_identity, adjust_identity
    end interface pix_is_identity

contains

    subroutine pix_histogram(src, hist)
        !! Counts per byte value for R, G, B (columns 1-3) and Rec. 709 luma (4).
        integer(int8), intent(in) :: src(:, :, :)
        integer(int64), intent(out) :: hist(0:255, 4)
        integer :: x, y, r, g, b
        hist = 0
        do y = 1, size(src, 3)
            do x = 1, size(src, 2)
                r = u8(src(1, x, y))
                g = u8(src(2, x, y))
                b = u8(src(3, x, y))
                hist(r, 1) = hist(r, 1) + 1
                hist(g, 2) = hist(g, 2) + 1
                hist(b, 3) = hist(b, 3) + 1
                hist(nint(luma(real(r), real(g), real(b))), 4) = &
                    hist(nint(luma(real(r), real(g), real(b))), 4) + 1
            end do
        end do
    end subroutine pix_histogram

    pure logical function levels_identity(p)
        type(pix_levels_t), intent(in) :: p
        levels_identity = all(p%black == 0.0) .and. all(p%white == 255.0) .and. all(p%gamma == 1.0)
    end function levels_identity

    pure logical function adjust_identity(p)
        type(pix_adjust_t), intent(in) :: p
        adjust_identity = p%exposure == 0.0 .and. p%brightness == 0.0 .and. p%contrast == 0.0 &
            .and. p%saturation == 0.0 .and. p%temperature == 0.0 .and. p%tint == 0.0 &
            .and. all(p%gains == 1.0)
    end function adjust_identity

    pure real function level(v, black, white, gamma)
        real, intent(in) :: v, black, white, gamma
        real :: t
        t = (v - black) / max(1.0, white - black)
        t = max(0.0, min(1.0, t))
        level = 255.0 * t**(1.0 / max(0.01, gamma))
    end function level

    subroutine pix_levels(src, p, dst)
        integer(int8), intent(in) :: src(:, :, :)
        type(pix_levels_t), intent(in) :: p
        integer(int8), allocatable, intent(out) :: dst(:, :, :)
        integer(int8) :: lut(0:255, 3)
        integer :: k, v, x, y
        do k = 1, 3
            do v = 0, 255
                lut(v, k) = to_byte_real(level(level(real(v), p%black(0), p%white(0), &
                    p%gamma(0)), p%black(k), p%white(k), p%gamma(k)))
            end do
        end do
        allocate(dst, mold=src)
        !$omp parallel do private(x, k)
        do y = 1, size(src, 3)
            do x = 1, size(src, 2)
                do k = 1, 3
                    dst(k, x, y) = lut(u8(src(k, x, y)), k)
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine pix_levels

    function pix_auto_levels(hist, clip) result(p)
        !! Per-channel black and white points that clip the fraction `clip` of
        !! pixels at each end (0.001 is a common choice).
        integer(int64), intent(in) :: hist(0:255, 4)
        real, intent(in) :: clip
        type(pix_levels_t) :: p
        integer(int64) :: total, acc, limit
        integer :: k, v
        do k = 1, 3
            total = sum(hist(:, k))
            if (total == 0) cycle
            limit = int(clip * real(total), int64)
            acc = 0
            do v = 0, 255
                acc = acc + hist(v, k)
                if (acc > limit) exit
            end do
            p%black(k) = real(min(v, 254))
            acc = 0
            do v = 255, 0, -1
                acc = acc + hist(v, k)
                if (acc > limit) exit
            end do
            p%white(k) = real(max(v, int(p%black(k)) + 1))
        end do
    end function pix_auto_levels

    function pix_neutral_gains(r, g, b) result(gains)
        !! Linear gains that make the sampled colour (bytes) neutral grey at
        !! its own luminance.
        integer, intent(in) :: r, g, b
        real :: gains(3), lin(3), y
        lin = srgb_to_linear(max(1.0, real([r, g, b])) / 255.0)
        y = 0.2126 * lin(1) + 0.7152 * lin(2) + 0.0722 * lin(3)
        gains = y / lin
    end function pix_neutral_gains

    subroutine pix_adjust(src, p, dst)
        !! Exposure, white balance, temperature and tint act on linear light;
        !! contrast, brightness and saturation on encoded values, as in common
        !! photo editors.
        integer(int8), intent(in) :: src(:, :, :)
        type(pix_adjust_t), intent(in) :: p
        integer(int8), allocatable, intent(out) :: dst(:, :, :)
        real :: lin(0:255), mul(3), c(3), y, contrast
        integer :: x, yy, v
        do v = 0, 255
            lin(v) = srgb_to_linear(real(v) / 255.0)
        end do
        mul = p%gains * 2.0**p%exposure
        mul(1) = mul(1) * (1.0 + 0.3 * p%temperature) * (1.0 + 0.15 * p%tint)
        mul(2) = mul(2) * (1.0 - 0.3 * p%tint)
        mul(3) = mul(3) * (1.0 - 0.3 * p%temperature) * (1.0 + 0.15 * p%tint)
        contrast = tan((0.25 + 0.2499 * max(-1.0, min(1.0, p%contrast))) * 3.14159265)
        allocate(dst, mold=src)
        !$omp parallel do private(x, c, y, v)
        do yy = 1, size(src, 3)
            do x = 1, size(src, 2)
                do v = 1, 3
                    c(v) = lin(u8(src(v, x, yy))) * mul(v)
                end do
                c = linear_to_srgb(c)
                c = (c - 0.5) * contrast + 0.5 + 0.5 * p%brightness
                y = luma(c(1), c(2), c(3))
                c = y + (c - y) * (1.0 + p%saturation)
                dst(1:3, x, yy) = to_byte_real(255.0 * c)
            end do
        end do
        !$omp end parallel do
    end subroutine pix_adjust
end module fortpix_tone
