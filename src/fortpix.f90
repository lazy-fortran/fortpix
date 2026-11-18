module fortpix
    !! Raster image operations on (3, width, height) sRGB byte arrays and
    !! (width, height) coverage masks.
    use fortpix_base, only: u8, to_byte, to_byte_real, srgb_to_linear, linear_to_srgb, luma
    use fortpix_geometry
    use fortpix_tone
    use fortpix_filter
    use fortpix_select
    use fortpix_draw
    implicit none
    public
end module fortpix
