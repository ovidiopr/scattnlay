<template>
  <div class="row items-baseline">
    <div class="col-xs-12 col-sm-auto text-weight-bold text-center q-px-md q-py-sm">
        <q-btn :loading="isRunning"
               :disable="isRunning||!isNmieLoaded"
               color="primary"
               no-caps
               :label="isNmieLoaded ? 'Run simulation' : 'Loading...'"
               @click="runSpectrumSimulation">
          <template #loading>
            <q-spinner-gears />
          </template>
        </q-btn>
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
  computed, watch,
} from 'vue'
import { useStore } from 'src/store'
import { range, rangeInt } from 'components/utils'
import { cloneDeep } from 'lodash'

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

    function getWLs(){
      const fromWL = $store.state.simulationSetup.current.fromWL
      const toWL = $store.state.simulationSetup.current.toWL
      const pointsWL = $store.state.simulationSetup.current.pointsWL
      const stepWL = (toWL-fromWL)/(pointsWL-1)
      const WLs = range(fromWL, toWL, stepWL);
      return WLs
    }

    function initQ(mode_n:number[], mode_types:number[]) {
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
      return {Qsca, Qabs, Qext, Qsca_n, Qabs_n, Qext_n}
    }

    //-------------------------------------------------------------------------//
    //---------------------  Main  --------------------------------------------//
    //-------------------------------------------------------------------------//
    function runSpectrumSimulation() {
      if (isRunning.value) {
        console.log('Some Nmie is already running!')
        return
      }
      isRunning.value = true
      setTimeout(()=> {

        $store.commit('simulationSetup/copySetupFromGuiToCurrent')

        const host = $store.state.simulationSetup.current.hostIndex
        const WLs = getWLs()
        const mode_n = rangeInt($store.state.simulationSetup.current.numberOfModesToPlot, 1);
        const mode_types = range(0, 1);
        let {Qsca, Qabs, Qext, Qsca_n, Qabs_n, Qext_n} = initQ(mode_n, mode_types)

        try {
          if (!$store.state.simulationSetup.nmie) throw 'ERROR! Scattnlay module was not loaded'
          const nmie = $store.state.simulationSetup.nmie
          const layers = cloneDeep($store.state.simulationSetup.current.layers)
          const nmieStartedTime = performance.now()
          for (const WL of WLs) {
            nmie.SetWavelength(WL)

            nmie.ClearTarget()
            for (const layer of layers) {
              if (layer.nSpline) layer.n = layer.nSpline.at(WL)
              if (layer.kSpline) layer.k = layer.kSpline.at(WL)
              nmie.AddTargetLayerReIm(layer.layerWidth * host, layer.n / host, layer.k / host)
            }

            nmie.SetModeNmaxAndType(-1, -1)
            nmie.RunMieCalculation()
            Qsca.push(nmie.GetQsca())
            Qabs.push(nmie.GetQabs())
            Qext.push(nmie.GetQext())

            for (const mode_type of mode_types) {
              for (const n of mode_n) {
                nmie.SetModeNmaxAndType(n, mode_type)
                nmie.RunMieCalculation()
                Qsca_n[mode_type][n - 1].push(nmie.GetQsca())
                Qabs_n[mode_type][n - 1].push(nmie.GetQabs())
                Qext_n[mode_type][n - 1].push(nmie.GetQext())
              }
            }
          }
          const nmieTotalRunTime = (performance.now()-nmieStartedTime)/1000
          // console.log('Total simulation time:', nmieTotalRunTime, 's')
          $store.commit('simulationSetup/setNmieTotalRunTime', nmieTotalRunTime)
          $store.commit('plotRuntime/setQ', {WLs, Qsca, Qabs, Qext, Qsca_n, Qabs_n, Qext_n})
          $store.commit('plotRuntime/updateNumberOfPlotsFromPreviousSimulations')
          $store.commit('plotRuntime/setCommonLabel', $store.state.simulationSetup.current.plotLabel)
          $store.commit('simulationSetup/setPlotLabel', '')
          $store.commit('plotRuntime/updateSpectraPlot')
        } catch (e) {
          console.log(e)
        }
        isRunning.value = false
      },50)
    }

    watch(isNmieLoaded, ()=>{
      if (isNmieLoaded.value) runSpectrumSimulation()
    })

    return { isRunning, isNmieLoaded,
      runSpectrumSimulation
    }
  },
})
</script>
