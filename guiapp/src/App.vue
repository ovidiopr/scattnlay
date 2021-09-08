<template>
  <router-view />
</template>
<script lang="ts">
import { defineComponent } from 'vue';


// const fs = require("fs");
// const compiled = new WebAssembly.Module("./nmiejs.wasm"));
// const imports = {};
// Object.defineProperty(module, "exports", {
//   get: () => new WebAssembly.Instance(compiled, imports).exports
// });


// import VueWasm from 'vue-wasm';
// import nmiejs from './nmiejs.wasm';
// const nmieModule = require("./nmiejs.wasm");
//
// const nmiejs = nmieModule().then(({ instance }:any) => {
//   return instance.exports
// });
//
// const nmie = new nmiejs.nmie();


import nmiejs from './nmiejs.js';
// Test nmiejs if working
void (async () => {
  // eslint-disable-next-line @typescript-eslint/no-unsafe-assignment
  const module = await nmiejs();
  // eslint-disable-next-line @typescript-eslint/no-unsafe-assignment
  const nmie = new module.nmie();
  nmie.ClearTarget();
  let R = 100.0;
  let reN = 4.0;
  let imN = 0.01;
  nmie.AddTargetLayerReIm(R, reN, imN)
  nmie.SetModeNmaxAndType(-1, -1);
  let WL = 800;
  nmie.SetWavelength(WL);
  nmie.RunMieCalculation();
  console.log(nmie.GetQsca());
  // outer_arc_points, radius_points, from_Rho, to_Rho,
  // from_Theta, to_Theta, from_Phi, to_Phi, isIgnoreAvailableNmax
  nmie.RunFieldCalculationPolar(2, 2, 0.1, 1.5, 0, 3.1415, 0, 3.1415, 0);
  console.log('Field Eabs:', nmie.GetFieldEabs());


})();


// nmie.ClearTarget();
// let R = 100.0;
// let reN = 4.0;
// let imN = 0.01;
// nmie.AddTargetLayerReIm(R, reN, imN)
// nmie.SetModeNmaxAndType(-1, -1);
// let WL = 800;
// nmie.SetWavelength(WL);
// nmie.RunMieCalculation();
// console.log(nmie.GetQsca());

// const extractModule = async (module) => {
//   const { instance } = await module();
//   return instance.exports;
// };
// const nmie = extractModule(nmieModule);
// const nmie = nmieModule().then(({ instance }) => {
//   return instance.exports
// });

export default defineComponent({
  name: 'App'
})
</script>
