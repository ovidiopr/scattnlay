<template>
  <div class="row items-baseline">
    <div class="col-xs-12 col-sm-auto text-weight-bold text-center q-px-md q-py-sm">
      <div :style="flexRowTitleStyle">
        Modes to plot
      </div>
    </div>
    <div class="col-xs-grow col-sm q-px-xs">
      <div class="row justify-xs-center justify-sm-start items-baseline">

        <div class="col-auto">
              <q-input
                  v-model.number="numberOfModesToPlot"
                  outlined
                  type="number"
                  min=1
                  max=100
                  step=1
                  dense
                  style="width: 6em"
              />
        </div>

      </div>
    </div>
  </div>
</template>

<script lang="ts">
import {
  defineComponent,
  computed,
  } from 'vue'
import { useStore } from 'src/store'
import { flexRowTitleStyle } from 'components/utils'

export default defineComponent({
  name: 'GetPlotSettings',

  setup() {
    const $store = useStore()

    const numberOfModesToPlot = computed({
      get: () => $store.state.simulationSetup.gui.numberOfModesToPlot,
      set: val => {
        let intVal = parseInt(val.toString())
        if (isNaN(intVal)) return
        if (intVal<1) intVal = 1
        if (intVal>10) intVal = 10
        $store.commit('simulationSetup/setNumberOfModesToPlot', intVal)
        $store.commit('plotRuntime/resizeIsPlotMode',intVal)
      }
    })

    return { flexRowTitleStyle,
      numberOfModesToPlot,
      }
  },
})
</script>
