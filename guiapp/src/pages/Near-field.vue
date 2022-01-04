<template>
  <q-page class="column q-px-md">
    <div class="q-ma-md"/>
    <ShowNearFieldWarning/>
    <RunSimulationSpectrum v-show="false"/>
    <GetWlFromPlot/>
    <div class="q-ma-xs"/>
    <GetNearFieldSettings/>
    <div class="q-ma-xs"/>
    <RunSimulationNearField/>
    <div class="q-ma-xs"/>
    <PlotNearField/>
    <div class="q-ma-xs"/>
    <GetNearFieldColorScale/>
<!--    <div class="col-auto">-->
<!--      Input result: {{$store.state.simulationSetup.gui.nearFieldSetup.plotSideResolution}}-->
<!--    </div>-->

    <div class="q-ma-lg"/>
    <!-- place QPageScroller at end of page -->
    <q-page-scroller position="bottom-right" :scroll-offset="150" :offset="[18, 18]">
      <q-btn size="xs" fab icon="keyboard_arrow_up" color="primary" />
    </q-page-scroller>
  </q-page>
</template>

<script lang='ts'>
import {
  defineComponent
} from 'vue'
import GetWlFromPlot from 'components/nearfield/GetWlFromPlot.vue'
import GetNearFieldSettings from 'components/nearfield/GetNearFieldSettings.vue'
import RunSimulationSpectrum from 'components/spectrum/RunSimulationSpectrum.vue'
import RunSimulationNearField from 'components/nearfield/RunSimulationNearField.vue'
import PlotNearField from 'components/nearfield/PlotNearField.vue'
import ShowNearFieldWarning from 'components/nearfield/ShowNearFieldWarning.vue'
import GetNearFieldColorScale from 'components/nearfield/GetNearFieldColorScale.vue'

// import { useStore } from 'src/store'


export default defineComponent({
  name: 'NearField',
  components: {
    GetNearFieldColorScale,
    ShowNearFieldWarning,
    PlotNearField,
    RunSimulationNearField,
    GetWlFromPlot, RunSimulationSpectrum,
    GetNearFieldSettings,
  },
  setup() {

    // const $store = useStore()

    return {
    }
  }
})
</script>
