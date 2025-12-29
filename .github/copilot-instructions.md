# Scattnlay Project Instructions

## Tech Stack
- **Core** C++ with Python bindings via pybind11.
  - **Build System:** CMake.
  - **WASM Compilation:** Emscripten with SIMD support.
  - **SIMD Library:** Google Highway.
  - **Multiprecision:** Boost Multiprecision library.
- **Testing:** Vitest for frontend, pytest for Python bindings, and Google Test for C++ core.
- **Frontend:** Vue 3 (Composition API), Quasar v2, Vite, TypeScript, using core WASM module.

## Coding Standards
- **Minimalism:** Touch as few files as possible when fixing or adding features.
- **TypeScript:** Use strict typing. Avoid `any`, especially for physical parameters.

## Testing Rules
Run
./clean-all.sh && ./build-all.sh && ./test-all.sh
before implementing changes to ensure no regressions and after changes to verify correctness. Fix all issues.

During development build and run tests frequently, use specific build/test commands for faster feedback loops, e.g. for C++ core:
```bash
cmake --build build --target test_bulk_sphere -j 8
ctest --output-on-failure -R test_bulk_sphere
```
for Python only development install C++ extension in editable mode and run specific file to verify changes:
```bash
python3 -m pip install -e .  # Once
python3 tests/test_py.py
```

## Other Guidelines

Never commit. Never use `rm` with any flags, use `./clean-all.sh` instead or `git clean -fd` to delete untracked files and directories.
