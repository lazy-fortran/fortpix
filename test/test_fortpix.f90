program test_fortpix
    !! Properties of the operations on generated images.
    use iso_fortran_env, only: int8, int64
    use fortpix
    implicit none
    integer(int8), allocatable :: img(:, :, :), a(:, :, :), b(:, :, :)
    real, allocatable :: m(:, :), m2(:, :)
    integer(int64) :: hist(0:255, 4)
    type(pix_levels_t) :: lv
    type(pix_adjust_t) :: adj
    integer :: x, y, x0, y0, x1, y1, failures
    real :: g(3)
    failures = 0

    ! Gradient test image: R = x, G = y, B = constant.
    allocate(img(3, 40, 30))
    do y = 1, 30
        do x = 1, 40
            img(:, x, y) = to_byte([5 * x, 7 * y, 100])
        end do
    end do

    call pix_crop(img, 3, 4, 10, 5, a)
    call check(size(a, 2) == 10 .and. size(a, 3) == 5, 'crop size')
    call check(all(a(:, 1, 1) == img(:, 3, 4)) .and. all(a(:, 10, 5) == img(:, 12, 8)), 'crop pixels')
    call pix_crop(img, 35, 28, 10, 10, a)
    call check(size(a, 2) == 6 .and. size(a, 3) == 3, 'crop clipped at the edge')

    call pix_rotate90(img, 1, a)
    call check(size(a, 2) == 30 .and. size(a, 3) == 40, 'quarter turn swaps sides')
    call check(all(a(:, 30, 1) == img(:, 1, 1)), 'clockwise: top-left goes top-right')
    call pix_rotate90(a, 3, b)
    call check(all(b == img), 'four quarter turns are the identity')
    call pix_flip(img, .true., a)
    call check(all(a(:, 1, :) == img(:, 40, :)), 'horizontal flip')
    call pix_flip(a, .true., b)
    call check(all(b == img), 'flip twice')

    call pix_rotate(img, 0.0, .true., [0, 0, 0], a)
    call check(all(shape(a) == shape(img)), 'zero rotation keeps the size')
    call check(maxval(abs(u8(a) - u8(img))) <= 1, 'zero rotation is the identity')
    call pix_rotate(img, 10.0, .true., [255, 0, 255], a)
    call check(abs(real(size(a, 2)) / size(a, 3) - 40.0 / 30.0) < 0.06, 'straightening keeps the aspect')
    call check(.not. any(u8(a(1, :, :)) == 255 .and. u8(a(2, :, :)) == 0 .and. u8(a(3, :, :)) == 255), &
        'straightened image has no fill corners')
    call check(abs(pix_straighten_scale(100, 100, 45.0) - 1.0 / sqrt(2.0)) < 1e-5, 'square at 45 degrees')
    call pix_rotate(img, 90.0, .false., [0, 0, 0], a)
    call check(size(a, 2) == 30 .and. size(a, 3) == 40, 'free rotation by 90 degrees')

    a = img
    a(:, :, :) = to_byte(77)
    call pix_resize(a, 13, 9, b)
    call check(all(u8(b) == 77), 'resampling preserves a constant')
    do y = 1, 30
        do x = 1, 40
            a(:, x, y) = to_byte(merge(0, 200, mod(x + y, 2) == 0))
        end do
    end do
    call pix_resize(a, 10, 7, b)
    call check(maxval(abs(u8(b) - 100)) <= 6, 'downsampling averages a checkerboard')
    call pix_resize(img, 80, 60, b)
    call check(abs(u8(b(1, 41, 30)) - u8(img(1, 21, 15))) <= 4, 'upsampling follows a ramp')

    call pix_histogram(img, hist)
    call check(all(sum(hist, dim=1) == 1200_int64), 'histogram counts every pixel')
    call check(hist(100, 3) == 1200_int64, 'constant channel in one bin')
    call pix_levels(img, lv, a)
    call check(all(a == img), 'default levels are the identity')
    lv = pix_auto_levels(hist, 0.0)
    call check(lv%black(1) == 5.0 .and. lv%white(1) == 200.0, 'auto levels find the range')
    call pix_levels(img, lv, a)
    call check(u8(a(1, 1, 1)) == 0 .and. u8(a(1, 40, 1)) == 255, 'levels stretch to full range')
    lv = pix_levels_t()
    lv%gamma(0) = 2.0
    call pix_levels(img, lv, a)
    call check(u8(a(3, 1, 1)) > 100, 'gamma above one brightens mid tones')

    call pix_adjust(img, adj, a)
    call check(maxval(abs(u8(a) - u8(img))) <= 1, 'neutral adjustment is the identity')
    adj%exposure = 1.0
    call pix_adjust(img, adj, a)
    call check(u8(a(3, 5, 5)) > 130, 'one stop brightens')
    adj = pix_adjust_t()
    adj%saturation = -1.0
    call pix_adjust(img, adj, a)
    call check(maxval(abs(u8(a(1, :, :)) - u8(a(2, :, :)))) <= 1, 'full desaturation is grey')
    g = pix_neutral_gains(200, 100, 50)
    adj = pix_adjust_t()
    adj%gains = g
    a = img
    a(:, 1, 1) = to_byte([200, 100, 50])
    call pix_adjust(a, adj, b)
    call check(maxval(u8(b(:, 1, 1))) - minval(u8(b(:, 1, 1))) <= 2, 'white balance neutralizes the sample')

    a = img
    a(:, :, :) = to_byte(90)
    call pix_unsharp(a, 2.0, 1.0, 0.0, b)
    call check(all(u8(b) == 90), 'sharpening a constant changes nothing')
    a(:, 21:, :) = to_byte(160)
    call pix_unsharp(a, 2.0, 1.0, 0.0, b)
    call check(u8(b(1, 20, 10)) < 90 .and. u8(b(1, 21, 10)) > 160, 'sharpening overshoots at an edge')
    call pix_blur(a, 3.0, b)
    call check(u8(b(1, 20, 10)) > 90 .and. u8(b(1, 21, 10)) < 160, 'blurring softens an edge')

    call pix_mask_rect(40, 30, 2.5, 3.0, 10.5, 8.0, m)
    call check(abs(sum(m) - 40.0) < 1e-3, 'rectangle area')
    call check(abs(m(3, 4) - 0.5) < 1e-6, 'rectangle edge coverage')
    call pix_mask_ellipse(40, 30, 0.0, 0.0, 40.0, 30.0, m)
    call check(abs(sum(m) - 3.14159265 * 20 * 15) < 6.0, 'ellipse area')
    call pix_mask_polygon(40, 30, [5.0, 35.0, 5.0], [5.0, 5.0, 25.0], m)
    call check(abs(sum(m) - 300.0) < 1.5, 'lasso triangle area')
    call pix_mask_bounds(m, x0, y0, x1, y1)
    call check(x0 == 6 .and. y0 == 6 .and. x1 == 35 .and. y1 == 25, 'mask bounds')
    call pix_mask_rect(40, 30, 0.0, 0.0, 20.0, 30.0, m2)
    call pix_mask_combine(m, m2, MASK_SUBTRACT)
    call check(abs(sum(m) - 75.0) < 1.5, 'subtracting a selection')

    a = img
    a(:, :, :) = to_byte(10)
    a(:, 5:15, 5:15) = to_byte(200)
    a(:, 30:35, 20:25) = to_byte(205)
    call pix_magic_wand(a, 8, 8, 10, .true., m)
    call check(nint(sum(m)) == 121, 'contiguous magic wand')
    call pix_magic_wand(a, 8, 8, 10, .false., m)
    call check(nint(sum(m)) == 157, 'global magic wand')
    call pix_magic_wand(a, 1, 1, 0, .true., m)
    call check(nint(sum(m)) == 1200 - 157, 'wand around holes')

    b = a
    b(:, :, :) = to_byte(250)
    call pix_mask_rect(40, 30, 0.0, 0.0, 10.0, 30.0, m)
    call pix_blend_masked(a, b, m, img)
    call check(all(u8(img(:, 1:10, :)) == 250) .and. all(u8(img(:, 11:, 1:4)) == 10), 'masked adjustment')

    a(:, :, :) = to_byte(255)
    call pix_stroke(a, [5.0, 35.0], [15.0, 15.0], [2.0, 2.0], [0, 0, 0])
    call check(u8(a(1, 20, 15)) == 0 .and. u8(a(1, 20, 10)) == 255, 'stroke centre and outside')
    call check(abs(real(count(u8(a(1, 20, :)) < 128)) - 4.0) <= 1.0, 'stroke width')
    b = a
    call pix_stroke(a, [20.0, 20.0, 20.0, 20.0], [15.0, 15.0, 15.0, 15.0], [3.0, 3.0, 3.0, 3.0], [0, 0, 0])
    a(:, :, :) = to_byte(255)
    call pix_stroke(a, [20.0, 20.0, 20.0, 20.0], [15.0, 15.0, 15.0, 15.0], [3.0, 3.0, 3.0, 3.0], [200, 0, 0])
    call check(u8(a(1, 22, 15)) == 200 .and. u8(a(2, 23, 17)) > 0, 'resting pen keeps one deposit')

    if (failures > 0) error stop 'fortpix tests failed'
    print '(a)', 'fortpix: all checks passed'
contains
    subroutine check(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            print '(2a)', 'FAIL: ', what
            failures = failures + 1
        end if
    end subroutine check
end program test_fortpix
