<template>
  <div :id="chartUUID"></div>
</template>

<script lang="ts">
import { newPlot, react, Data, Layout, Config } from 'plotly.js-dist-min'
import {
  defineComponent,
  PropType,
  watch,
  onMounted,
  onUnmounted,
} from 'vue'
import {cloneDeep} from 'lodash'
import { v4 as uuidv4 } from 'uuid'

interface PlotlyChart {
  data: Data[],
  layout: Layout,
  config: Config|undefined
}
export default defineComponent({
  name: 'ReactiveChart',
  props: {
    chart: {
      type: Object as PropType<PlotlyChart>,
      required: true,
    },
  },
  setup(props) {
    const chartUUID = uuidv4()
    let chartLocal = cloneDeep(props.chart)

    // Update (or add if absent) width and height of the layout to fit current window
    // and replot it.
    function plotlyReact () {
      const width = window.innerWidth
      const height = window.innerHeight*0.8
      chartLocal.layout.width = width * 0.92
      chartLocal.layout.height = height * 0.95
      if (height < 600) chartLocal.layout.height = height

      // react(...) is a promise, but we do not care to await it, so mark it with `void` keyword
      if (chartLocal.config == undefined) {
        void react(chartUUID, chartLocal.data, chartLocal.layout)
      } else {
        void react(chartUUID, chartLocal.data, chartLocal.layout, chartLocal.config)
      }
    }

    onMounted( async () => {
      if (chartLocal.config == undefined) {
        await newPlot(chartUUID, chartLocal.data, chartLocal.layout)
      } else {
        await newPlot(chartUUID, chartLocal.data, chartLocal.layout, chartLocal.config)
      }
    })

    window.addEventListener('resize', plotlyReact)
    onUnmounted(()=>{
      window.removeEventListener('resize', plotlyReact)
    })

    watch(props, ()=>{
      chartLocal = cloneDeep(props.chart)
      plotlyReact()
    })

    return {
      chartUUID
    }
  }
})
</script>
