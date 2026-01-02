# Overall prirorities

# libs
- [X] c++20
- [X] nanobind 

# Algorithms?

- [X] Du algorithm for large complex refractive indices?
    - works with all examples from Du paper



# Core
- [ ] complete SIMD operations for far-field calculations
  - [X] Qsca\Qext\Qback
  - [ ] S1\S2
- [X] complete SIMD operations for nearfield calculations 

- [ ] verify SIMD speedup for WASM/Emscripten builds (both nearfield and far-field)

- [ ] implemnt multicore processing for spectra\nearfield\far-field
  calculations
  - No need to paralleize sing mie calculation over several hosts, so
    multicore should be sufficient -> OpenMP?
  - OpenMP seems to be supported by WASM\Emscripten?

 - [ ] verify multicore speedup for WASM/Emscripten builds (both nearfield and far-field)


# Testing

- [ ] use all examples from original Ovidio paper
- [ ] use table 2 from Du paper (S1/S2)
- [ ] add adapter for miepython API, reuse miepython  [test suite](https://github.com/scottprahl/miepython/tree/main/tests)
- [ ] compare performance with miepython



# Ideas for Improvements


## Performance Improvements

- ?      


## Akima Spline

Interpolation Strategy: Currently, you use cubic-spline-ts. Cubic splines can sometimes "overshoot" (Runge's phenomenon) for materials with sharp absorption edges (like some organics in the UV). Adding an option for Akima splines or simple linear interpolation would give users more control over material accuracy.