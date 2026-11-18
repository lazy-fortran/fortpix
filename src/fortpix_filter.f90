module fortpix_filter
    !! Gaussian blur and unsharp masking of RGB images.
    use iso_fortran_env, only: int8
    use fortpix_base, only: u8, to_byte_real
    implicit none
    private
    public :: pix_blur, pix_unsharp

contains

    subroutine blur_real(img, sigma)
        !! Separable Gaussian with clamped edges, in place on (3, w, h).
        real, intent(inout) :: img(:, :, :)
        real, intent(in) :: sigma
        real, allocatable :: k(:), line(:, :)
        integer :: r, w, h, x, y, j
        real :: acc(3)
        if (sigma <= 0.0) return
        r = max(1, ceiling(3.0 * sigma))
        allocate(k(-r:r))
        do j = -r, r
            k(j) = exp(-0.5 * (real(j) / sigma)**2)
        end do
        k = k / sum(k)
        w = size(img, 2)
        h = size(img, 3)
        !$omp parallel do private(line, x, j, acc)
        do y = 1, h
            line = img(:, :, y)
            do x = 1, w
                acc = 0.0
                do j = -r, r
                    acc = acc + k(j) * line(:, max(1, min(w, x + j)))
                end do
                img(:, x, y) = acc
            end do
        end do
        !$omp end parallel do
        !$omp parallel do private(line, y, j, acc)
        do x = 1, w
            line = img(:, x, :)
            do y = 1, h
                acc = 0.0
                do j = -r, r
                    acc = acc + k(j) * line(:, max(1, min(h, y + j)))
                end do
                img(:, x, y) = acc
            end do
        end do
        !$omp end parallel do
    end subroutine blur_real

    subroutine pix_blur(src, sigma, dst)
        integer(int8), intent(in) :: src(:, :, :)
        real, intent(in) :: sigma
        integer(int8), allocatable, intent(out) :: dst(:, :, :)
        real, allocatable :: img(:, :, :)
        img = real(u8(src(1:3, :, :)))
        call blur_real(img, sigma)
        allocate(dst(3, size(src, 2), size(src, 3)))
        dst = to_byte_real(img)
    end subroutine pix_blur

    subroutine pix_unsharp(src, radius, amount, threshold, dst)
        !! out = in + amount (in - blur(in)) where |in - blur| exceeds
        !! `threshold` (0..255); radius is the Gaussian sigma in pixels.
        integer(int8), intent(in) :: src(:, :, :)
        real, intent(in) :: radius, amount, threshold
        integer(int8), allocatable, intent(out) :: dst(:, :, :)
        real, allocatable :: orig(:, :, :), soft(:, :, :), detail(:, :, :)
        orig = real(u8(src(1:3, :, :)))
        soft = orig
        call blur_real(soft, radius)
        detail = orig - soft
        where (abs(detail) < threshold) detail = 0.0
        allocate(dst(3, size(src, 2), size(src, 3)))
        dst = to_byte_real(orig + amount * detail)
    end subroutine pix_unsharp
end module fortpix_filter
