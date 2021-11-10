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
    <GetNumberOfModes/>
    <div class="q-ma-xs"/>
    <RunSimulationSpectrum/>
    <div class="q-ma-xs"/>
    <ReactiveChart :chart="chart"/>
    <div class="col-auto q-pa-md">
      Input result: {{$store.state.plotRuntime.Qsca}}
    </div>
  </q-page>
</template>

<script lang='ts'>
import {
  defineComponent,
    ref
} from 'vue'
import GetUnits from 'components/GetUnits.vue'
import GetHostIndex from 'components/GetHostIndex.vue'
import GetSourceParameters from 'components/GetSourceParameters.vue'
import GetParticleParameters from 'components/GetParticleParameters.vue'
import GetNumberOfModes from 'components/GetNumberOfModes.vue'
import RunSimulationSpectrum from 'components/RunSimulationSpectrum.vue'
import ReactiveChart from 'components/ReactiveChart.vue'
// import { useStore } from 'src/store'


export default defineComponent({
  name: 'PageIndex',
  components: {
    RunSimulationSpectrum, GetUnits, GetHostIndex, GetSourceParameters,
    GetParticleParameters, GetNumberOfModes,
    ReactiveChart
  },
  setup() {
    const chart = ref({
      data: [{
        y: [],
        line: {
          color: '#5e9e7e',
          width: 4,
          shape: 'line'
        }
      }],
      layout: {
        margin: {
          l: 0,
          r: 40,
          b: 50,
          t: 0
        },
        // paper_bgcolor: '#7f7f7f',
        // plot_bgcolor: '#c7c7c7',
        // title: 'reactive charts',
        xaxis: {
          // will be set on mount
          title: ''
        },
        yaxis: {
          title: 'Normalized cross-sections'
        },
        showlegend: true,
        legend: {
          orientation: 'h',
          x: -.1,
          y: 1.05
        },
      }
    })
    setTimeout(()=>{
      chart.value.layout.yaxis.title='some new title'
    }, 3000)
    // const $store = useStore()
    return {chart}
  }
})
</script>
