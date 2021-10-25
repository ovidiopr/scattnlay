<template>
  <div class="row items-baseline">
    <div class="col-xs-12 col-sm-auto text-weight-bold text-center q-px-md q-py-sm">
        <q-btn :loading="isRunning"
               color="primary"
               no-caps
               label="Run simulation"
               @click="runSimulation"
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
import { flexRowTitleStyle } from 'components/utils'

export default defineComponent({
  name: 'RunSimulationSpectrum',

  setup() {
    const $store = useStore()

    const isRunning = ref(false)

    const numberOfModesToPlot = computed({
      get: () => $store.state.simulationSetup.gui.numberOfModesToPlot,
      set: val => {
        const intVal = parseInt(val.toString())
        if (!isNaN(intVal)) $store.commit('simulationSetup/setNumberOfModesToPlot', intVal)
      }
    })

    function runSimulation() {
      isRunning.value = true
// simulate a delay
      setTimeout(() => {
        // we're done, we reset loading state
        isRunning.value = false
      }, 3000)
    }

    return { flexRowTitleStyle,
      numberOfModesToPlot, runSimulation,
      isRunning
      }
  },
})
</script>
