<template>
  <div>
    <ReactiveChart :chart="nearFieldPlot"/>
  </div>
</template>

<script lang="ts">
import ReactiveChart from 'components/ReactiveChart.vue'
import { useStore } from 'src/store'
import {
  defineComponent,
  // onActivated,
  // ref,
    reactive,
  computed,
  watch
} from 'vue'
// import { flexRowTitleStyle } from 'components/config'
import { plotlyChart } from 'src/store/plot-runtime/state'
import { PlotData, /*, DataTitle*/ } from 'plotly.js-dist-min'
import {nearFieldPlane} from 'src/store/simulation-setup/state';


export default defineComponent({
  name: 'PlotNearField',
  components: {
    ReactiveChart,
  },
  setup () {
    const $store = useStore()
    const nearFieldPlotInit:plotlyChart = {
      data: [],
      layout: {
        margin: {
          l: 0,
          r: 40,
          b: 50,
          t: 30
        },
        xaxis: {
          title: ''
        },
        yaxis: {
          scaleanchor: 'x',
          title: ''
        },
        showlegend: false,
      },
      config: {responsive: true,
        // showEditInChartStudio: true,
        displaylogo: false}
    }

    const nearFieldPlot = reactive(nearFieldPlotInit)
    const crossSection = computed(()=>$store.state.simulationSetup.current.nearFieldSetup.crossSection)
    const nearFieldEk = computed( ()=>$store.state.plotRuntime.nearFieldEk)
    const nearFieldHk = computed( ()=>$store.state.plotRuntime.nearFieldHk)
    const nearFieldEH = computed( ()=>$store.state.plotRuntime.nearFieldEH)
    watch([nearFieldEk, nearFieldHk, nearFieldEH], ()=>{
      nearFieldPlot.data.length = 0
      const heatMapSettings: Partial<PlotData> = {type: 'heatmap',  colorscale: 'Jet'}
      let heatMapData: Partial<PlotData> = {}
      if (crossSection.value == nearFieldPlane.Ek) heatMapData  = { z: nearFieldEk.value}
      if (crossSection.value == nearFieldPlane.Hk) heatMapData  = { z: nearFieldHk.value}
      if (crossSection.value == nearFieldPlane.EH) heatMapData  = { z: nearFieldEH.value}
      nearFieldPlot.data.push({...heatMapData, ...heatMapSettings})

    })
    // const xaxisTitle = computed(()=>{
    //   let title:string|Partial<DataTitle> = ''
    //   if ($store.state.plotRuntime.spectrumPlots.layout.xaxis?.title) {
    //     title = $store.state.plotRuntime.spectrumPlots.layout.xaxis.title
    //   }
    //   return title
    // })
    // if (nearFieldPlot.layout.xaxis) nearFieldPlot.layout.xaxis.title = xaxisTitle.value
    // watch( xaxisTitle, ()=>{
    //   if (nearFieldPlot.layout.xaxis) nearFieldPlot.layout.xaxis.title = xaxisTitle.value
    // })
    //
    // const plotRange=ref('material data')
    // const isPlotReN = ref(true)
    // const isPlotImN = ref(true)
    // const isPlotInterpolation = ref(true)
    //
    // const fromWL = computed(()=>$store.state.simulationSetup.gui.fromWL)
    // const toWL = computed(()=>$store.state.simulationSetup.gui.toWL)
    // const pointsWL = computed(()=>$store.state.simulationSetup.gui.pointsWL)
    //
    // // updatenearFieldPlot used to be a Vuex mutation with a single argument. As
    // // long as it is in the component feel free as many arguments as needed.
    // function updatenearFieldPlot( val:{activatedMaterials:material[], sourceUnits:string,
    //     fromWL:number, toWL:number, pointsWL:number, plotRange: string,
    //     isPlotReN: boolean, isPlotImN: boolean, isPlotInterpolation: boolean })
    // {
    //   nearFieldPlot.data.length = 0
    //   for (const material of val.activatedMaterials) {
    //     if (!material.isPlot) continue
    //
    //     if (!material.nSpline) continue
    //     if (val.isPlotReN) {
    //       const traceDataReN: Partial<Data> = {
    //         x: material.nSpline.xs.map(x => toUnits(x, val.sourceUnits)),
    //         y: material.nSpline.ys,
    //         mode: 'markers',
    //         marker:{size:5},
    //         type: 'scatter',
    //         name: 'Re(n) ' + material.name + ' data'
    //       }
    //       nearFieldPlot.data.push(traceDataReN)
    //     }
    //
    //     if (!material.kSpline) continue
    //     if (val.isPlotImN) {
    //       const traceDataImN: Partial<Data> = {
    //         x: material.kSpline.xs.map(x => toUnits(x, val.sourceUnits)),
    //         y: material.kSpline.ys,
    //         mode: 'markers',
    //         marker:{size:5},
    //         type: 'scatter',
    //         name: 'Im(n) ' + material.name + ' data'
    //       }
    //       nearFieldPlot.data.push(traceDataImN)
    //     }
    //
    //     if (val.isPlotInterpolation) {
    //       let fromWL = val.fromWL
    //       let toWL = val.toWL
    //       let pointsWL = val.pointsWL-1
    //       if (nearFieldPlot.layout.xaxis) nearFieldPlot.layout.xaxis.range = [fromWL, toWL]
    //       if (val.plotRange == 'material data') {
    //         fromWL = material.nSpline.xs[0]
    //         toWL = material.nSpline.xs[material.nSpline.xs.length-1]
    //         pointsWL = 1000
    //         if (nearFieldPlot.layout.xaxis) nearFieldPlot.layout.xaxis.range = undefined
    //       }
    //       const stepWL = (toWL-fromWL)/pointsWL
    //
    //       const WLs:number[] =[]
    //       const nSpline:number[] = []
    //       const kSpline:number[] = []
    //       for (let i=0; i<pointsWL; ++i) {
    //         WLs.push(fromWL+i*stepWL)
    //         nSpline.push(material.nSpline.at(fromWL+i*stepWL))
    //         kSpline.push(material.kSpline.at(fromWL+i*stepWL))
    //       }
    //       WLs.push(toWL)
    //       nSpline.push(material.nSpline.at(toWL))
    //       kSpline.push(material.kSpline.at(toWL))
    //
    //       if (val.isPlotReN) {
    //         if (!material.nSpline) continue
    //         const traceDataReNi: Partial<Data> = {
    //           x: WLs.map(x => toUnits(x, val.sourceUnits)),
    //           y: nSpline,
    //           type: 'scatter',
    //           name: 'Re(n) ' + material.name + ' interpolate'
    //         }
    //         nearFieldPlot.data.push(traceDataReNi)
    //       }
    //
    //       if (val.isPlotImN) {
    //         if (!material.nSpline) continue
    //         const traceDataImNi: Partial<Data> = {
    //           x: WLs.map(x => toUnits(x, val.sourceUnits)),
    //           y: kSpline,
    //           type: 'scatter',
    //           name: 'Im(n) ' + material.name + ' interpolate'
    //         }
    //         nearFieldPlot.data.push(traceDataImNi)
    //       }
    //
    //     }
    //   }
    // }
    //
    // function updateSpectraPlot() {
    //   updatenearFieldPlot(
    //       {
    //         activatedMaterials: $store.state.guiRuntime.activatedMaterials,
    //         sourceUnits: sourceUnits.value,
    //         fromWL: fromWL.value,
    //         toWL: toWL.value,
    //         pointsWL: pointsWL.value,
    //         plotRange: plotRange.value,
    //         isPlotReN: isPlotReN.value,
    //         isPlotImN: isPlotImN.value,
    //         isPlotInterpolation: isPlotInterpolation.value
    //       })
    // }
    //
    // updateSpectraPlot()
    //
    // const materialsToPlot=computed(()=>$store.state.guiRuntime.activatedMaterials
    //     .filter((val)=>val.isPlot)
    //     .map(val=>val.name))
    //
    // watch ([materialsToPlot, plotRange, isPlotReN, isPlotImN, isPlotInterpolation], ()=>{
    //   updateSpectraPlot()
    // })
    //
    // onActivated(()=>updateSpectraPlot())

    return {
      nearFieldPlot,

    }
  }
})
</script>
