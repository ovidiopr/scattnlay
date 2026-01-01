# Overall prirorities

# libs
- [ ] c++20
- [ ] nanobind 

# Algorithms?

- [X] Du algorithm for large complex refractive indices?
    - works with all examples from Du paper



# Core
- [X] complete SIMD operations for far-field calculations
- [ ] complete SIMD operations for nearfield calculations 
- [ ] implemnt multicore processing for spectra\nearfield\far-field
  calculations
  - No need to paralleize sing mie calculation over several hosts, so
    multicore should be sufficient -> OpenMP?
  - OpenMP seems to be supported by WASM\Emscripten?



# Testing

- [ ] use all examples from original Ovidio paper
- [ ] use table 2 from Du paper (S1/S2)
- [ ] add adapter for miepython API, reuse miepython  [test suite](https://github.com/scottprahl/miepython/tree/main/tests)
- [ ] compare performance with miepython



# Ideas for Improvements

## Performance Improvements
Remove np.atleast_1d inside loops: In scattnlay/main.py, you
repeatedly call np.vstack. For large spectra (e.g., 10,000 points),
this leads to O(N2) memory copying. It is better to pre-allocate the
results array:
code Python

        
    # Instead of vstack
    Qext = np.empty(x.shape[0])
    for i, xi in enumerate(x):
        # ... compute ...
        Qext[i] = res

      

## Dependency Management

    pyproject.toml: The dynamic = ["dependencies"] tag is present, but usually, it's better to explicitly list numpy and pybind11 under dependencies to ensure a clean install from PyPI without needing a requirements.txt.

## Akima Spline

Interpolation Strategy: Currently, you use cubic-spline-ts. Cubic splines can sometimes "overshoot" (Runge's phenomenon) for materials with sharp absorption edges (like some organics in the UV). Adding an option for Akima splines or simple linear interpolation would give users more control over material accuracy.