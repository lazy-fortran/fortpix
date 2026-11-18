module fortpix_base
    !! Byte conversion shared by the operations. Images are integer(int8)
    !! arrays (3, width, height) holding sRGB bytes; masks are real (width,
    !! height) coverage from 0 to 1.
    use iso_fortran_env, only: int8
    implicit none
    private
    public :: u8, to_byte, to_byte_real, srgb_to_linear, linear_to_srgb, luma

contains

    elemental integer function u8(b)
        integer(int8), intent(in) :: b
        u8 = iand(int(b), 255)
    end function u8

    elemental integer(int8) function to_byte(v)
        integer, intent(in) :: v
        integer :: c
        c = max(0, min(255, v))
        if (c > 127) c = c - 256
        to_byte = int(c, int8)
    end function to_byte

    elemental integer(int8) function to_byte_real(v)
        !! Rounds a value on the 0..255 scale.
        real, intent(in) :: v
        to_byte_real = to_byte(nint(max(0.0, min(255.0, v))))
    end function to_byte_real

    elemental real function srgb_to_linear(c)
        !! c and the result on 0..1.
        real, intent(in) :: c
        if (c <= 0.04045) then
            srgb_to_linear = c / 12.92
        else
            srgb_to_linear = ((c + 0.055) / 1.055)**2.4
        end if
    end function srgb_to_linear

    elemental real function linear_to_srgb(c)
        real, intent(in) :: c
        real :: x
        x = max(0.0, min(1.0, c))
        if (x <= 0.0031308) then
            linear_to_srgb = 12.92 * x
        else
            linear_to_srgb = 1.055 * x**(1.0 / 2.4) - 0.055
        end if
    end function linear_to_srgb

    elemental real function luma(r, g, b)
        !! Rec. 709 weights on encoded values, as image editors show luminosity.
        real, intent(in) :: r, g, b
        luma = 0.2126 * r + 0.7152 * g + 0.0722 * b
    end function luma
end module fortpix_base
