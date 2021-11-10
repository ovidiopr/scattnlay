<template>
  <div id="plotly_chart"></div>
</template>

<script lang="ts">
// import * as Plotly from 'plotly.js'
import { newPlot, Data, Layout, Config } from 'plotly.js'

import {
  defineComponent,
  PropType,
  onMounted,
  onUnmounted,
} from 'vue'

interface PlotlyChart {
  uuid: string,
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

    function  handleResize() {
      const width = window.innerWidth
      const height = window.innerHeight*0.8
      // props.chart.layout.width = width * 0.92
      // props.chart.layout.height = height * 0.95
      // // if (width < 600) props.chart.layout.width = width
      // if (height < 600) props.chart.layout.height = height
      console.log(width, height)
      console.log('layout', props.chart.layout.width, props.chart.layout.height)
    }

    window.addEventListener('resize', handleResize)
    handleResize()
    onUnmounted(()=>{
      window.removeEventListener('resize', handleResize)
    })
    onMounted( async () => {
      if (props.chart.config == undefined) {
        await newPlot('plotly_chart', props.chart.data, props.chart.layout)
      } else {
        await newPlot('plotly_chart', props.chart.data, props.chart.layout, props.chart.config)
      }
    })
    }
})
</script>
