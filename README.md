# FortPix

Raster image operations in Fortran for photo editing.
[Source](https://github.com/lazy-fortran/fortpix).

Images are `integer(int8)` arrays `(3, width, height)` of sRGB bytes. Masks are
`real` arrays `(width, height)` of coverage from 0 to 1, with pixel `(i, j)`
covering `[i-1, i] x [j-1, j]`.

| Area | Procedures |
|---|---|
| Geometry | `pix_crop`, `pix_rotate90`, `pix_flip`, `pix_rotate` (bicubic, optional straightening crop), `pix_resize` (Lanczos-3) |
| Tone | `pix_histogram`, `pix_levels`, `pix_auto_levels`, `pix_adjust` (exposure, white balance, temperature, tint, brightness, contrast, saturation), `pix_neutral_gains` |
| Filters | `pix_blur`, `pix_unsharp` |
| Selections | `pix_mask_rect`, `pix_mask_ellipse`, `pix_mask_polygon`, `pix_magic_wand`, `pix_mask_combine`, `pix_mask_bounds`, `pix_blend_masked` |
| Ink | `pix_stroke`: union of round segments with varying radius, coloured once |

The library has no dependencies and builds with gfortran and LFortran, natively
and for WebAssembly. Loops are parallel with OpenMP when enabled. File formats
belong to [FortImg](https://github.com/lazy-fortran/fortimg).

```sh
fo test        # or: fpm test
```

MIT licence.
