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
    const isPlotQabs = computed( ()=>$store.state.plotRuntime.isPlotQabs)
    watch(isPlotQabs, ()=>{
      console.log('isPlotQabs = ', isPlotQabs.value)
      $store.commit('plotRuntime/updateSpectraPlot')
    })
    return {}
  }
})
</script>
