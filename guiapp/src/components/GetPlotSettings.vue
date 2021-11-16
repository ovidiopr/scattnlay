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
          <!--          TODO move to config.ts -> min max
           -->
              <q-input
                  v-model.number="numberOfModesToPlot"
                  outlined
                  type="number"
                  min=1
                  max=10
                  step=1
                  dense
                  style="width: 6em"
              />
        </div>
        <div class="col-auto">
<!--          TODO move to config.ts ->                     style="width: 17.7em"
 -->
          <q-input  v-model="plotLabel"
                    label="plot label"
                    :shadow-text="plotLabel=='' ? shadowText+' (Tab to complete)' : ''"
                    outlined
                    dense
                    class="q-pa-sm"
                    style="width: 17.7em"
                    @keydown="processInputFill"
                    @focus="processInputFill"
          />
        </div>
        <div class="col-auto">
          <q-checkbox v-model="isAppendPlots">
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
  ref
  } from 'vue'
import { event } from 'quasar'
const { stopAndPrevent } = event

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
      if (intVal>10) numberOfModesToPlot.value = 10
    })

    const inputFillCancelled = ref(false)

    return { flexRowTitleStyle,
      numberOfModesToPlot,
      plotLabel, isAppendPlots,
      shadowText,
      processInputFill (e:any) {
        if (e === void 0) {
          return
        }

        if (e.keyCode === 27) {
          if (inputFillCancelled.value !== true) {
            inputFillCancelled.value = true
          }
        }
        else if (e.keyCode === 9) {
          if (inputFillCancelled.value !== true && shadowText.value.length > 0) {
            stopAndPrevent(e)
            if (plotLabel.value == '') plotLabel.value = shadowText.value
          }
        }
        else if (inputFillCancelled.value === true) {
          inputFillCancelled.value = false
        }
      },
      }
  },
})
</script>
