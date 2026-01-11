import createModule from '../../build_wasm/src/nmiejs.js';

async function run() {
  const module = await createModule();
  const WL = 800.0;
  const res = 256; // Grid resolution

  const runBench = (EngineClass, label) => {
    if (!EngineClass) {
        console.log(`${label}: Not available`);
        return null;
    }
    const nmie = new EngineClass();
    nmie.SetWavelength(WL);
    // Si-Ag-Si model (cumulative radii)
    nmie.AddTargetLayerReIm(0.139, 3.69, 0.006);  // Si
    nmie.AddTargetLayerReIm(0.322, 0.14, 5.29);   // Ag
    nmie.AddTargetLayerReIm(0.503, 3.69, 0.006);  // Si

    const start = performance.now();
    // 256x256, Plane Ek, relative size 2.0
    nmie.RunFieldCalculationCartesian(res, res, 2.0, 0, 0, 0, 0, true);
    const end = performance.now();

    return { time: (end - start), data: nmie.GetFieldEabs() };
  };

  const scalar = runBench(module.nmie_scalar, "Scalar");
  const simd = runBench(module.nmie, "SIMD");

  if (scalar && simd) {
      console.log(`Scalar: ${scalar.time.toFixed(2)}ms`);
      console.log(`SIMD:   ${simd.time.toFixed(2)}ms`);
      console.log(`Speedup: ${(scalar.time / simd.time).toFixed(2)}x`);

      // Parity check
      let diff = 0;
      for(let i=0; i<scalar.data.length; i++) diff += Math.abs(scalar.data[i] - simd.data[i]);
      console.log(`Total numerical diff: ${diff.toExponential(2)}`);
      console.log(`Average per element: ${(diff / scalar.data.length).toExponential(2)}`);
  }
}
run();
