<template>
  <div>
    <ReactiveChart :chart="$store.state.plotRuntime.spectraPlot"/>
  </div>
</template>

<script>
import ReactiveChart from 'components/ReactiveChart.vue'
import { useStore } from 'src/store'
import {
  defineComponent,
  // ref,
  // reactive,
  computed,
  watch
} from 'vue'

export default defineComponent({
  name: 'PlotSpectra',
  components: {
    ReactiveChart,
  },
  setup () {
    const $store = useStore()
    const isPlotQsca = computed( ()=>$store.state.plotRuntime.isPlotQsca)
    watch(isPlotQsca, ()=>$store.commit('plotRuntime/updateSpectraPlot'))
    const isPlotQabs = computed( ()=>$store.state.plotRuntime.isPlotQabs)
    watch(isPlotQabs, ()=>$store.commit('plotRuntime/updateSpectraPlot'))
    const isPlotQext = computed( ()=>$store.state.plotRuntime.isPlotQext)
    watch(isPlotQext, ()=>$store.commit('plotRuntime/updateSpectraPlot'))

    const isPlotQscaTotal = computed( ()=>$store.state.plotRuntime.isPlotQscaTotal)
    watch(isPlotQscaTotal, ()=>$store.commit('plotRuntime/updateSpectraPlot'))
    const isPlotQabsTotal = computed( ()=>$store.state.plotRuntime.isPlotQabsTotal)
    watch(isPlotQabsTotal, ()=>$store.commit('plotRuntime/updateSpectraPlot'))
    const isPlotQextTotal = computed( ()=>$store.state.plotRuntime.isPlotQextTotal)
    watch(isPlotQextTotal, ()=>$store.commit('plotRuntime/updateSpectraPlot'))

    const isPlotModeE = computed( ()=>$store.state.plotRuntime.isPlotModeE)
    watch(isPlotModeE, ()=>$store.commit('plotRuntime/updateSpectraPlot'), { deep: true })
    const isPlotModeH = computed( ()=>$store.state.plotRuntime.isPlotModeH)
    watch(isPlotModeH, ()=>$store.commit('plotRuntime/updateSpectraPlot'), { deep: true })

    return {}
  }
})
</script>
