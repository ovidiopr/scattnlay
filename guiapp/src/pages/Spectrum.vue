<template>
  <q-page class="column q-px-md">
    <div class="q-ma-md" />
    <GetHostIndex />
    <div class="q-ma-xs" />
    <GetUnits />
    <div id="GetSourceParametersHref" class="q-ma-xs" />
    <GetSourceParameters />
    <div class="q-ma-xs" />
    <GetParticleParameters />
    <div class="q-ma-xs" />
    <GetModesToPlot />
    <div class="q-ma-xs" />
    <GetPlotSettings />
    <div class="q-ma-xs" />
    <RunSimulationSpectrum />
    <!--    <div class="col-auto q-pa-md">-->
    <!--      Input result: {{$store.state.plotRuntime.isPlotQabs}}-->
    <!--    </div>-->
    <div class="q-ma-sm" />
    <PlotSelector />
    <div class="q-ma-xs" />
    <PlotSpectra />

    <!-- place QPageScroller at end of page -->
    <div class="q-ma-lg" />
    <q-page-scroller
      position="bottom-right"
      :scroll-offset="150"
      :offset="[18, 18]"
    >
      <q-btn size="xs" fab icon="keyboard_arrow_up" color="primary" />
    </q-page-scroller>
  </q-page>
</template>

<script lang="ts">
import { defineComponent } from 'vue';

import GetUnits from 'components/spectrum/GetUnits.vue';
import GetHostIndex from 'components/spectrum/GetHostIndex.vue';
import GetSourceParameters from 'components/spectrum/GetSourceParameters.vue';
import GetParticleParameters from 'components/spectrum/GetParticleParameters.vue';
import GetModesToPlot from 'components/spectrum/GetModesToPlot.vue';
import RunSimulationSpectrum from 'components/spectrum/RunSimulationSpectrum.vue';
import GetPlotSettings from 'components/spectrum/GetPlotSettings.vue';
import PlotSelector from 'components/spectrum/PlotSelector.vue';
import PlotSpectra from 'components/spectrum/PlotSpectra.vue';
// import { useStore } from 'src/store'

export default defineComponent({
  name: 'SpectrumPage',
  components: {
    GetModesToPlot,
    RunSimulationSpectrum,
    GetUnits,
    GetHostIndex,
    GetSourceParameters,
    GetParticleParameters,
    GetPlotSettings,
    PlotSelector,
    PlotSpectra,
  },
  // setup() {
  //   const $store = useStore()
  // }
});
</script>
