<template>
  <q-page class="column q-px-md">
    <div class="q-ma-md"/>
    <GetHostIndex/>
    <div class="q-ma-xs"/>
    <GetUnits/>
    <div class="q-ma-xs"/>
    <GetSourceParameters/>
    <div class="q-ma-xs"/>
    <GetParticleParameters/>
    <div class="q-ma-xs"/>
    <GetPlotSettings/>
    <div class="q-ma-xs"/>
    <RunSimulationSpectrum/>
    <div class="q-ma-xs"/>
    <PlotSelector/>
    <div class="q-ma-xs"/>
    <ReactiveChart :chart="$store.state.plotRuntime.spectraPlot"/>
    <div class="col-auto q-pa-md">
      Input result: {{$store.state.plotRuntime.isQscaToggle}}
    </div>
  </q-page>
</template>

<script lang='ts'>
import {
  defineComponent,
} from 'vue'

import GetUnits from 'components/GetUnits.vue'
import GetHostIndex from 'components/GetHostIndex.vue'
import GetSourceParameters from 'components/GetSourceParameters.vue'
import GetParticleParameters from 'components/GetParticleParameters.vue'
import GetPlotSettings from 'components/GetPlotSettings.vue'
import RunSimulationSpectrum from 'components/RunSimulationSpectrum.vue'
import ReactiveChart from 'components/ReactiveChart.vue'
import PlotSelector from 'components/PlotSelector.vue'
// import { useStore } from 'src/store'


export default defineComponent({
  name: 'PageIndex',
  components: {
    RunSimulationSpectrum, GetUnits, GetHostIndex,
    GetSourceParameters, GetParticleParameters,
    GetPlotSettings, ReactiveChart, PlotSelector
  },
  // setup() {
  //
  // }
})
</script>
