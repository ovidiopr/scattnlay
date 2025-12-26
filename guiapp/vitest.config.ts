import { defineConfig } from 'vitest/config';

export default defineConfig({
  test: {
    environment: 'node', // WASM engine logic is math-heavy, no DOM needed
    include: ['test/**/*.test.ts'],
  },
});
