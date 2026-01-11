import createNmiejsModule from '../../build_wasm/src/nmiejs.js';

async function test() {
  const instance = await createNmiejsModule();
  const nmie = new instance.nmie();
  nmie.SetWavelength(500);
  nmie.AddTargetLayerReIm(100, 1.5, 0.01);
  nmie.RunMieCalculation();
  console.log("Qsca:", nmie.GetQsca());
  console.log("WASM Module Loaded and Executed Successfully.");
}
test().catch(console.error);
