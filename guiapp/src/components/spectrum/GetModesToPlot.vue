<template>
  <div class="row items-baseline">
    <div class="col-xs-12 col-sm-auto text-weight-bold text-center q-px-md q-py-sm">
      <div :style="flexRowTitleStyle">
        Modes to plot
      </div>
    </div>
    <div class="col-xs-grow col-sm q-px-xs">
      <div class="row justify-xs-center justify-sm-start items-center">

        <div class="col-auto">
              <q-input
                  v-model.number="numberOfModesToPlot"
                  outlined
                  type="number"
                  min=1
                  :max='maxModes'
                  step=1
                  dense
                  :style="'width: '+basicSelectorWidthStyle"
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
  watch,
  } from 'vue'

import { useStore } from 'src/store'
import { flexRowTitleStyle,
  maxNumberOfModesToPlot as maxModes,
  basicSelectorWidthStyle,
    basicWidthStyle
} from 'components/config'

export default defineComponent({
  name: 'GetModesToPlot',

  setup() {
    const $store = useStore()
    const numberOfModesToPlot = computed({
      get: () => $store.state.simulationSetup.gui.numberOfModesToPlot,
      set: val => {
        let intVal = parseInt(val.toString())
        if (isNaN(intVal)) return
        if (intVal<1) intVal = 1
        if (intVal>maxModes) intVal = maxModes
        $store.commit('simulationSetup/setNumberOfModesToPlot', intVal)
        $store.commit('plotRuntime/resizeSelectorIsPlotMode',intVal)
      }
    })

    watch(numberOfModesToPlot, ()=>{
      const intVal = parseInt(numberOfModesToPlot.value.toString())
      if (isNaN(intVal)) return
      if (intVal<1) numberOfModesToPlot.value = 1
      if (intVal>maxModes) numberOfModesToPlot.value = maxModes
    })

    return { flexRowTitleStyle, basicSelectorWidthStyle, basicWidthStyle,
      numberOfModesToPlot, maxModes,

    }
  },
})
</script>
