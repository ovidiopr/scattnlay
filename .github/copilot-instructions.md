# Scattnlay Project Instructions

## Tech Stack
- **Frontend:** Vue 3 (Composition API), Quasar v2, Vite, TypeScript.
- **Engine:** C++ compiled to WASM via Emscripten, with SIMD support using Google Highway, Boost for Multiprecision and pybind11 for Python extension.
- **State Management:** Transitioning from Vuex to Pinia (prioritize Pinia for new code).
- **Testing:** Vitest (Node environment for WASM/Physics logic).


## Coding Standards
- **Minimalism:** Touch as few files as possible when fixing or adding features.
- **TypeScript:** Use strict typing. Avoid `any`, especially for physical parameters.

## Testing Rules
Run
./clean-all.sh && ./build-all.sh && ./test-all.sh
before implementing changes to ensure no regressions and after changes to verify correctness. Fix all issues.

## Other Guidelines

Never commit. Never use `rm` with any flags, use `./clean-all.sh` instead or `git clean -fd` for untracked files.
