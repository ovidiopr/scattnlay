<template>
  <div class="row items-baseline">
    <div class="col-xs-12 col-sm-auto text-weight-bold text-center q-px-md q-py-sm">
        <q-btn :loading="isRunning"
               :disable="isRunning||!isNmieLoaded"
               color="primary"
               no-caps
               :label="isNmieLoaded ? 'Run simulation' : 'Loading...'"
               @click="runSpectrumSimulation"
      />
    </div>
    <div class="col-xs-grow col-sm q-px-xs">
      <div class="row justify-xs-center justify-sm-start items-baseline">

        <div class="col-auto">
        </div>

      </div>
    </div>
  </div>
</template>

<script lang="ts">
import {
  defineComponent,
  ref,
  computed,
  } from 'vue'
import { useStore } from 'src/store'
import { range, rangeInt } from 'components/utils'

export default defineComponent({
  name: 'RunSimulationSpectrum',

  setup() {
    const $store = useStore()

    const isRunning = computed({
      get: ()=> $store.state.simulationSetup.isNmieRunning,
      set: val => {
        val ? $store.commit('simulationSetup/markNmieAsStarted') : $store.commit('simulationSetup/markNmieAsFinished')
      }
    })

    const isNmieLoaded = computed(()=>{ return $store.state.simulationSetup.isNmieLoaded })

      //-----------------------------------------------------------------------//
      //-------------------  Main  --------------------------------------------//
      //-----------------------------------------------------------------------//
    function runSpectrumSimulation() {
        isRunning.value = true
        $store.commit('simulationSetup/copySetupFromGuiToCurrent')

        const host = $store.state.simulationSetup.current.hostIndex

        const fromWL = $store.state.simulationSetup.current.fromWL
        const toWL = $store.state.simulationSetup.current.toWL
        const pointsWL = $store.state.simulationSetup.current.pointsWL
        const stepWL = (toWL-fromWL)/(pointsWL-1)
        const WLs = range(fromWL, toWL, stepWL);

        const total_mode_n = $store.state.simulationSetup.current.numberOfModesToPlot
        const mode_n = rangeInt(total_mode_n, 1);
        const mode_types = range(0, 1);

        let Qsca:number[] = [], Qabs:number[] = [], Qext:number[] = []
        let Qsca_n:number[][][] = [[], []]
        let Qabs_n:number[][][] = [[], []]
        let Qext_n:number[][][] = [[], []]

        mode_types.forEach(function (mode_type) {
          mode_n.forEach(function () {
            Qsca_n[mode_type].push([])
            Qabs_n[mode_type].push([])
            Qext_n[mode_type].push([])
          })
        })
        try {
          if (!$store.state.simulationSetup.nmie) throw 'ERROR! Scattnlay module was not loaded'
          const nmie = $store.state.simulationSetup.nmie

          for (const WL of WLs) {
            nmie.ClearTarget()
            for (const layer of $store.state.simulationSetup.current.layers) {
              let reN = layer.n
              let imN = layer.k
              if (layer.nSpline) reN = layer.nSpline.at(WL)
              if (layer.kSpline) imN = layer.kSpline.at(WL)
              nmie.AddTargetLayerReIm(  layer.layerWidth*host,
                  reN/host, imN/host);
            }
            nmie.SetModeNmaxAndType(-1, -1);
            nmie.SetWavelength(WL);
            nmie.RunMieCalculation();
            Qsca.push(nmie.GetQsca());
            Qabs.push(nmie.GetQabs());
            Qext.push(nmie.GetQsca()+nmie.GetQabs());
            mode_types.forEach(function (mode_type) {
              mode_n.forEach(function (n) {
                nmie.SetModeNmaxAndType(n, mode_type);
                nmie.RunMieCalculation();
                Qsca_n[mode_type][n - 1].push(nmie.GetQsca());
                Qabs_n[mode_type][n - 1].push(nmie.GetQabs());
                Qext_n[mode_type][n - 1].push(nmie.GetQext());
              });
            });

          }
        } catch (e) {
          console.log(e)
        }
        console.log(Qsca)
        isRunning.value = false
      }
    return { isRunning, isNmieLoaded,
      runSpectrumSimulation
    }
  },
})
</script>
