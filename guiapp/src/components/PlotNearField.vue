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
import { toUnits,
  getMaxFromHeatmap,
  getMinFromHeatmap,
  limitMap
} from 'components/utils'
import { plotlyChart } from 'src/store/plot-runtime/state'
import { PlotData, DataTitle } from 'plotly.js-dist-min'
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
        shapes: [],
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
    const relativePlotSize = computed(()=>$store.state.simulationSetup.current.nearFieldSetup.relativePlotSize)
    const plotSideResolution = computed(()=>$store.state.simulationSetup.current.nearFieldSetup.plotSideResolution)
    const layerWidths = computed(()=>$store.state.simulationSetup.current.layers.map(x=>x.layerWidth))
    const totalR = computed(()=>layerWidths.value.reduce((a, b) => a + b, 0))
    const x0 = computed(()=>totalR.value*relativePlotSize.value)
    const dx = computed(()=>x0.value*2.0/(plotSideResolution.value-1))
    const units = computed(()=>$store.state.guiRuntime.units)
    const xy = computed(()=> {
      let x:number[] = []
      let y:number[] = []
      for (let i = 0; i< plotSideResolution.value; ++i) {
        for (let j = 0; j< plotSideResolution.value; ++j) {
          x.push(toUnits(-x0.value + j * dx.value, units.value))
          y.push(toUnits(-x0.value + i * dx.value, units.value))
        }
      }
      return {x:x, y:y}
    })

    const limitFrom = computed( ()=>$store.state.plotRuntime.nearFieldLimitFrom )
    const limitTo = computed(()=>$store.state.plotRuntime.nearFieldLimitTo )

    const nearFieldEk = computed( ()=>$store.state.plotRuntime.nearFieldEk)
    const nearFieldHk = computed( ()=>$store.state.plotRuntime.nearFieldHk)
    const nearFieldEH = computed( ()=>$store.state.plotRuntime.nearFieldEH)
    const nearFieldProc = computed( ()=>{
      let nearFieldStore = nearFieldEk.value
      if (crossSection.value == nearFieldPlane.Hk) nearFieldStore = nearFieldHk.value
      if (crossSection.value == nearFieldPlane.EH) nearFieldStore = nearFieldEH.value
      if (!nearFieldStore) return nearFieldStore
      $store.commit('plotRuntime/setNearFieldDataTo',getMaxFromHeatmap(nearFieldStore))
      $store.commit('plotRuntime/setNearFieldDataFrom',getMinFromHeatmap(nearFieldStore))
      return limitMap(nearFieldStore, limitFrom.value, limitTo.value)
      // return nearFieldStore.map(x=>x>limitTo.value?limitTo.value:x)
    })
    watch([nearFieldProc, xy], ()=>{
      nearFieldPlot.data.length = 0
      const heatMapSettings: Partial<PlotData> = {type: 'heatmap',
        colorscale: 'Jet', colorbar:{title:'|𝐸|∕|𝐸𝜊|'},
        z: nearFieldProc.value
      }
      nearFieldPlot.data.push({...xy.value, ...heatMapSettings})

      if (nearFieldPlot.layout.shapes) {
        nearFieldPlot.layout.shapes.length = 0
        let r = 0
        for (let widthLayer of layerWidths.value) {
          r += widthLayer
          nearFieldPlot.layout.shapes.push({
            type: 'circle',
            xref: 'x',
            yref: 'y',
            x0: toUnits(-r, units.value),
            y0: toUnits(-r, units.value),
            x1: toUnits(r, units.value),
            y1: toUnits(r, units.value),
            line: {
              color: 'rgba(255, 255, 255, 0.7)',
              width: 3
            }
          })
        }
      }
    })

    const xaxisTitle = computed(()=>{
      let title:string|Partial<DataTitle> = 'x ['+units.value+']'
      return title
    })
    if (nearFieldPlot.layout.xaxis) nearFieldPlot.layout.xaxis.title = xaxisTitle.value
    watch( xaxisTitle, ()=>{
      if (nearFieldPlot.layout.xaxis) nearFieldPlot.layout.xaxis.title = xaxisTitle.value
    })

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
      nearFieldPlot, totalR

    }
  }
})
</script>
