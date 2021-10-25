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

    return { isRunning, isNmieLoaded,
      //-----------------------------------------------------------------------//
      //-------------------  Main  --------------------------------------------//
      //-----------------------------------------------------------------------//
      runSpectrumSimulation() {
        isRunning.value = true
        $store.commit('simulationSetup/copySetupFromGuiToCurrent')

        const host = $store.state.simulationSetup.current.hostIndex

        const fromWL = $store.state.simulationSetup.current.fromWL
        const toWL = $store.state.simulationSetup.current.toWL
        const pointsWL = $store.state.simulationSetup.current.pointsWL
        const stepWL = (toWL-fromWL)/(pointsWL-1)
        let WLs = range(fromWL, toWL, stepWL);

        const total_mode_n = $store.state.simulationSetup.current.numberOfModesToPlot

        isRunning.value = false
      }
    }
  },
})
</script>
