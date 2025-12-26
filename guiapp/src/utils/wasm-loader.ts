import createNmiejsModule from '../../public/wasm/nmiejs.js';
import { nmieModule } from '../nmiejs';

let instance: nmieModule | null = null;

export async function getWasmModule(): Promise<nmieModule> {
  if (instance) return instance;

  // LocateFile tells Emscripten where to find the .wasm file at runtime
  instance = await createNmiejsModule({
    locateFile: (path: string) => {
      if (path.endsWith('.wasm')) {
        return `${process.env.publicPath || '/'}wasm/${path}`;
      }
      return path;
    }
  });

  return instance!;
}
