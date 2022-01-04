<template>
  <div>
    <ReactiveChart :chart="$store.state.plotRuntime.spectrumPlots"/>
  </div>
</template>

<script>
import ReactiveChart from 'guiapp/src/components/ReactiveChart.vue'
import { useStore } from 'guiapp/src/store'
import {
  defineComponent,
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
    watch(isPlotQsca, ()=>$store.commit('plotRuntime/updateSpectrumPlots'))
    const isPlotQabs = computed( ()=>$store.state.plotRuntime.isPlotQabs)
    watch(isPlotQabs, ()=>$store.commit('plotRuntime/updateSpectrumPlots'))
    const isPlotQext = computed( ()=>$store.state.plotRuntime.isPlotQext)
    watch(isPlotQext, ()=>$store.commit('plotRuntime/updateSpectrumPlots'))

    const isPlotQscaTotal = computed( ()=>$store.state.plotRuntime.isPlotQscaTotal)
    watch(isPlotQscaTotal, ()=>$store.commit('plotRuntime/updateSpectrumPlots'))
    const isPlotQabsTotal = computed( ()=>$store.state.plotRuntime.isPlotQabsTotal)
    watch(isPlotQabsTotal, ()=>$store.commit('plotRuntime/updateSpectrumPlots'))
    const isPlotQextTotal = computed( ()=>$store.state.plotRuntime.isPlotQextTotal)
    watch(isPlotQextTotal, ()=>$store.commit('plotRuntime/updateSpectrumPlots'))

    const isPlotModeE = computed( ()=>$store.state.plotRuntime.isPlotModeE)
    watch(isPlotModeE, ()=>$store.commit('plotRuntime/updateSpectrumPlots'), { deep: true })
    const isPlotModeH = computed( ()=>$store.state.plotRuntime.isPlotModeH)
    watch(isPlotModeH, ()=>$store.commit('plotRuntime/updateSpectrumPlots'), { deep: true })

    const sourceUnits = computed( ()=>$store.state.guiRuntime.sourceUnits)
    function setPlotTitle() {
      let title=''
      if (sourceUnits.value.endsWith('Hz')) {
        title = 'Frequency [' + sourceUnits.value + ']'
      } else if (sourceUnits.value.endsWith('eV')) {
        title = 'Energy [' + sourceUnits.value + ']'
      } else if (sourceUnits.value.endsWith('s')) {
        title = 'Period [' + sourceUnits.value + ']'
      } else {
        title = 'Wavelength [' + sourceUnits.value + ']'
      }
      $store.commit('plotRuntime/updateXAxisTitle', title)
    }

    function updateSpectraPlotsUnits(){
      setPlotTitle()
      $store.commit('plotRuntime/setWLsInUnits', sourceUnits.value)
      $store.commit('plotRuntime/updateSpectrumPlots')
    }

    updateSpectraPlotsUnits()
    watch(sourceUnits, ()=> updateSpectraPlotsUnits())

    return {}
  }
})
</script>
