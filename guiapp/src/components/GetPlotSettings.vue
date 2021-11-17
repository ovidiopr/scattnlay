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
  <div class="row items-center">
    <div class="col-xs-12 col-sm-auto text-center q-px-md q-py-sm">
      <div :style="flexRowTitleStyle">
        <span class="text-weight-bold">Plot label</span> (optional)
      </div>
    </div>
    <div class="col-xs-grow col-sm q-px-xs">
      <div class="row justify-xs-center justify-sm-start items-center">
        <div class="col-auto">
          <q-input  v-model="plotLabel"
                    :shadow-text="plotLabel=='' ? shadowText+' (Tab to complete)' : ''"
                    outlined
                    dense
                    class="q-py-sm"
                    :style="'width: '+basicWidthStyle"
                    @keydown="processInputFill"
                    @focus="processInputFill"
          />
        </div>
        <div class="col-auto q-pa-sm">
          <q-checkbox v-model="isAppendPlots" size="sm">
            append new spectra to the plot
          </q-checkbox>
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
import { event } from 'quasar'
const { stopAndPrevent } = event

import { useStore } from 'src/store'
import { flexRowTitleStyle,
  maxNumberOfModesToPlot as maxModes,
  basicSelectorWidthStyle,
    basicWidthStyle
} from 'components/config'

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
        if (intVal>maxModes) intVal = maxModes
        $store.commit('simulationSetup/setNumberOfModesToPlot', intVal)
        $store.commit('plotRuntime/resizeSelectorIsPlotMode',intVal)
      }
    })

    const plotLabel = computed({
      get: ()=> $store.state.simulationSetup.gui.plotLabel,
      set: val => $store.commit('simulationSetup/setPlotLabel', val)
    })

    const isAppendPlots = computed({
      get: ()=> $store.state.plotRuntime.isAppendPlots,
      set: val => $store.commit('plotRuntime/setIsAppendPlots', val)
    })

    const shadowText = computed(()=>{
      const numberOfLayers = $store.state.simulationSetup.gui.layers.length
      let particleType = 'bulk'
      if (numberOfLayers == 2) return 'core-shell'
      if (numberOfLayers > 2) return 'multilayer'
      const R = $store.state.simulationSetup.gui.layers[0].layerWidth.toString()
      const units = $store.state.guiRuntime.units
      return particleType+' R='+R+units
    })

    watch(numberOfModesToPlot, ()=>{
      const intVal = parseInt(numberOfModesToPlot.value.toString())
      if (isNaN(intVal)) return
      if (intVal<1) numberOfModesToPlot.value = 1
      if (intVal>maxModes) numberOfModesToPlot.value = maxModes
    })

    return { flexRowTitleStyle, basicSelectorWidthStyle, basicWidthStyle,
      numberOfModesToPlot, maxModes,
      plotLabel, isAppendPlots,
      shadowText,

      processInputFill (e:KeyboardEvent) {
        if (e === void 0) {
          return
        }
        if (e.code === 'Tab') {
          if (shadowText.value.length > 0 && plotLabel.value == '') {
            stopAndPrevent(e)
            plotLabel.value = shadowText.value
          }
        }
      },

    }
  },
})
</script>
