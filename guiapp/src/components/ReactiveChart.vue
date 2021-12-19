<template>
  <div :id="chartUUID"></div>
</template>

<script lang="ts">
import { newPlot, react } from 'plotly.js-dist-min'
import { plotlyChart } from 'src/store/plot-runtime/state'
import {
  defineComponent,
  PropType,
  onActivated,
  watch,
  onMounted,
  onUnmounted,
} from 'vue'
import {cloneDeep} from 'lodash'
import { v4 as uuidv4 } from 'uuid'

export default defineComponent({
  name: 'ReactiveChart',
  props: {
    chart: {
      type: Object as PropType<plotlyChart>,
      required: true,
    },
  },
  emits: [
    'settingID',
  ],

  setup(props, {emit}) {
    const chartUUID = uuidv4()
    // const chartUUID = 'plotly chart'
    emit('settingID', chartUUID)
    let chartLocal = cloneDeep(props.chart)

    // Update (or add if absent) width and height of the layout to fit current window
    // and replot it.
    function plotlyReact () {
      if (!document.getElementById(chartUUID)) return
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

    function plotlyNew() {
      if (chartLocal.config == undefined) {
        void newPlot(chartUUID, chartLocal.data, chartLocal.layout)
      } else {
        void newPlot(chartUUID, chartLocal.data, chartLocal.layout, chartLocal.config)
      }
    }
    onMounted( () => {
      plotlyNew()
    })

    window.addEventListener('resize', plotlyReact)
    onUnmounted(()=>{
      window.removeEventListener('resize', plotlyReact)
    })

    watch(props, ()=>{
      chartLocal = cloneDeep(props.chart)
      plotlyReact()
    })

    onActivated(()=>{
      chartLocal = cloneDeep(props.chart)
      plotlyReact()
    })

    return {
      chartUUID
    }
  }
})
</script>
