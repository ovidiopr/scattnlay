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
    <div class="col-auto q-pa-md">
      Input result: {{$store.state.simulationSetup.gui.plotLabel}}
    </div>
    <div class="q-ma-sm"/>
    <PlotSelector/>
    <div class="q-ma-xs"/>
    <PlotSpectra/>
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
import PlotSelector from 'components/PlotSelector.vue'
import PlotSpectra from 'components/PlotSpectra.vue'
// import { useStore } from 'src/store'


export default defineComponent({
  name: 'PageIndex',
  components: {
    RunSimulationSpectrum, GetUnits, GetHostIndex,
    GetSourceParameters, GetParticleParameters,
    GetPlotSettings, PlotSelector,
    PlotSpectra
  },
  // setup() {
  //   const $store = useStore()
  // }
})
</script>
