<template>
  <q-page class="column q-px-md">
    <div class="q-ma-md"/>
    <GetHostIndex/>
    <div class="q-ma-xs"/>
    <GetUnits/>
    <div id="GetSourceParametersHref" class="q-ma-xs"/>
    <GetSourceParameters/>
    <div class="q-ma-xs"/>
    <GetParticleParameters/>
    <div class="q-ma-xs"/>
    <GetModesToPlot/>
    <div class="q-ma-xs"/>
    <GetPlotSettings/>
    <div class="q-ma-xs"/>
    <RunSimulationSpectrum/>
<!--    <div class="col-auto q-pa-md">-->
<!--      Input result: {{$store.state.plotRuntime.isPlotQabs}}-->
<!--    </div>-->
    <div class="q-ma-sm"/>
    <PlotSelector/>
    <div class="q-ma-xs"/>
    <PlotSpectra/>

    <!-- place QPageScroller at end of page -->
    <div class="q-ma-lg"/>
    <q-page-scroller position="bottom-right" :scroll-offset="150" :offset="[18, 18]">
      <q-btn size="xs" fab icon="keyboard_arrow_up" color="primary" />
    </q-page-scroller>
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
import GetModesToPlot from 'components/GetModesToPlot.vue'
import RunSimulationSpectrum from 'components/RunSimulationSpectrum.vue'
import GetPlotSettings from 'components/GetPlotSettings.vue'
import PlotSelector from 'components/PlotSelector.vue'
import PlotSpectra from 'components/PlotSpectra.vue'
// import { useStore } from 'src/store'


export default defineComponent({
  name: 'SpectrumPage',
  components: {
    GetModesToPlot,
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
