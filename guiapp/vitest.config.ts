import { defineConfig } from 'vitest/config';
import { resolve } from 'path';

export default defineConfig({
  resolve: {
    alias: {
      src: resolve(__dirname, './src'),
      components: resolve(__dirname, './src/components'),
    },
  },
  test: {
    environment: 'node', // WASM engine logic is math-heavy, no DOM needed
    include: ['test/**/*.test.ts'],
  },
});
