<template>
  <q-page class="column q-px-md">
    <div class="q-ma-md"/>
    <RunSimulationSpectrum v-show="false"/>
    <GetWlFromPlot/>
    <div class="q-ma-xs"/>
    <GetNearFieldSettings/>
    <div class="q-ma-xs"/>
<!--    <RunSimulationNearField/>-->
    <div class="col-auto">
      Input result: {{$store.state.simulationSetup.gui.nearFieldSetup.crossSection}}
    </div>

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
import GetWlFromPlot from 'components/GetWlFromPlot.vue'
import GetNearFieldSettings from 'components/GetNearFieldSettings.vue'
import RunSimulationSpectrum from 'components/RunSimulationSpectrum.vue'
// import RunSimulationNearField from 'components/RunSimulationNearField.vue'

// import { useStore } from 'src/store'


export default defineComponent({
  name: 'NearField',
  components: {
    // RunSimulationNearField,
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
