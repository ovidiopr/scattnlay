<template>
  <div>
    <div class="row items-baseline">
      <div class="col-xs-12 col-sm-auto text-weight-bold text-center q-px-md q-py-sm">
        <div :style="flexRowTitleStyle">
          Plot options
        </div>
      </div>
      <div class="col-xs-grow col-sm">
        <div class="row justify-xs-center justify-sm-start items-center">
          <div class="col-auto" >select data: </div>
          <div class="col-auto" ><q-toggle v-model="isPlotReN">Re(n)</q-toggle></div>
          <div class="col-auto" ><q-toggle v-model="isPlotImN">Im(n)</q-toggle></div>
          <div class="col-auto" ><q-toggle v-model="isPlotInterpolation">Interpolation</q-toggle></div>

        </div>
      </div>
    </div>

    <div class="row items-baseline">
      <div class="col-xs-12 col-sm-auto text-weight-bold text-center q-px-md q-py-sm">
        <div :style="flexRowTitleStyle">

        </div>
      </div>
      <div class="col-xs-grow col-sm">
        <div class="row justify-xs-center justify-sm-start items-baseline">
          <div class="col-auto q-pr-md" >interpolation range from: </div>
          <div class="col-auto">
            <div class="q-gutter-x-md">
              <q-radio v-model="plotRange" dense size='sm' val="material data" label="material data" />
              <q-radio v-model="plotRange" dense size='sm' val="simulation settings" label="simulation settings" />
            </div>
          </div>
        </div>
      </div>
    </div>

    <ReactiveChart :chart="materialPlots"/>
  </div>
</template>

<script lang="ts">
import ReactiveChart from 'components/ReactiveChart.vue'
import { useStore } from 'src/store'
import {
  defineComponent,
  onActivated,
  ref,
    reactive,
  computed,
  watch
} from 'vue'
import { flexRowTitleStyle } from 'components/config'
import { plotlyChart } from 'src/store/plot-runtime/state'
import { material } from 'src/store/simulation-setup/state'
import { Data, DataTitle } from 'plotly.js-dist-min'
import { toUnits } from 'components/utils'


export default defineComponent({
  name: 'PlotMaterials',
  components: {
    ReactiveChart,
  },
  setup () {
    const $store = useStore()
    const materialPlotsInit:plotlyChart = {
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
          title: 'Refractive index'
        },
        showlegend: true,
        legend: {
          orientation: 'h',
          x: -.1,
          y: 1.05
        },
      },
      config: {responsive: true,
        // showEditInChartStudio: true,
        displaylogo: false}
    }

    const materialPlots = reactive(materialPlotsInit)

    const sourceUnits = computed( ()=>$store.state.guiRuntime.sourceUnits)

    const xaxisTitle = computed(()=>{
      let title:string|Partial<DataTitle> = ''
      if ($store.state.plotRuntime.spectrumPlots.layout.xaxis?.title) {
        title = $store.state.plotRuntime.spectrumPlots.layout.xaxis.title
      }
      return title
    })
    if (materialPlots.layout.xaxis) materialPlots.layout.xaxis.title = xaxisTitle.value
    watch( xaxisTitle, ()=>{
      if (materialPlots.layout.xaxis) materialPlots.layout.xaxis.title = xaxisTitle.value
    })

    const plotRange=ref('material data')
    const isPlotReN = ref(true)
    const isPlotImN = ref(true)
    const isPlotInterpolation = ref(true)

    const fromWL = computed(()=>$store.state.simulationSetup.gui.fromWL)
    const toWL = computed(()=>$store.state.simulationSetup.gui.toWL)
    const pointsWL = computed(()=>$store.state.simulationSetup.gui.pointsWL)

    // updateMaterialPlots used to be a Vuex mutation with a single argument. As
    // long as it is in the component feel free as many arguments as needed.
    function updateMaterialPlots( val:{activatedMaterials:material[], sourceUnits:string,
        fromWL:number, toWL:number, pointsWL:number, plotRange: string,
        isPlotReN: boolean, isPlotImN: boolean, isPlotInterpolation: boolean })
    {
      materialPlots.data.length = 0
      for (const material of val.activatedMaterials) {
        if (!material.isPlot) continue

        if (!material.nSpline) continue
        if (val.isPlotReN) {
          const traceDataReN: Partial<Data> = {
            x: material.nSpline.xs.map(x => toUnits(x, val.sourceUnits)),
            y: material.nSpline.ys,
            mode: 'markers',
            marker:{size:5},
            type: 'scatter',
            name: 'Re(n) ' + material.name + ' data'
          }
          materialPlots.data.push(traceDataReN)
        }

        if (!material.kSpline) continue
        if (val.isPlotImN) {
          const traceDataImN: Partial<Data> = {
            x: material.kSpline.xs.map(x => toUnits(x, val.sourceUnits)),
            y: material.kSpline.ys,
            mode: 'markers',
            marker:{size:5},
            type: 'scatter',
            name: 'Im(n) ' + material.name + ' data'
          }
          materialPlots.data.push(traceDataImN)
        }

        if (val.isPlotInterpolation) {
          let fromWL = val.fromWL
          let toWL = val.toWL
          let pointsWL = val.pointsWL-1
          if (materialPlots.layout.xaxis) materialPlots.layout.xaxis.range = [fromWL, toWL]
          if (val.plotRange == 'material data') {
            fromWL = material.nSpline.xs[0]
            toWL = material.nSpline.xs[material.nSpline.xs.length-1]
            pointsWL = 1000
            if (materialPlots.layout.xaxis) materialPlots.layout.xaxis.range = undefined
          }
          const stepWL = (toWL-fromWL)/pointsWL

          const WLs:number[] =[]
          const nSpline:number[] = []
          const kSpline:number[] = []
          for (let i=0; i<pointsWL; ++i) {
            WLs.push(fromWL+i*stepWL)
            nSpline.push(material.nSpline.at(fromWL+i*stepWL))
            kSpline.push(material.kSpline.at(fromWL+i*stepWL))
          }
          WLs.push(toWL)
          nSpline.push(material.nSpline.at(toWL))
          kSpline.push(material.kSpline.at(toWL))

          if (val.isPlotReN) {
            if (!material.nSpline) continue
            const traceDataReNi: Partial<Data> = {
              x: WLs.map(x => toUnits(x, val.sourceUnits)),
              y: nSpline,
              type: 'scatter',
              name: 'Re(n) ' + material.name + ' interpolate'
            }
            materialPlots.data.push(traceDataReNi)
          }

          if (val.isPlotImN) {
            if (!material.nSpline) continue
            const traceDataImNi: Partial<Data> = {
              x: WLs.map(x => toUnits(x, val.sourceUnits)),
              y: kSpline,
              type: 'scatter',
              name: 'Im(n) ' + material.name + ' interpolate'
            }
            materialPlots.data.push(traceDataImNi)
          }

        }
      }
    }

    function updateSpectraPlot() {
      updateMaterialPlots(
          {
            activatedMaterials: $store.state.guiRuntime.activatedMaterials,
            sourceUnits: sourceUnits.value,
            fromWL: fromWL.value,
            toWL: toWL.value,
            pointsWL: pointsWL.value,
            plotRange: plotRange.value,
            isPlotReN: isPlotReN.value,
            isPlotImN: isPlotImN.value,
            isPlotInterpolation: isPlotInterpolation.value
          })
    }

    updateSpectraPlot()

    const materialsToPlot=computed(()=>$store.state.guiRuntime.activatedMaterials
        .filter((val)=>val.isPlot)
        .map(val=>val.name))

    watch ([materialsToPlot, plotRange, isPlotReN, isPlotImN, isPlotInterpolation], ()=>{
      updateSpectraPlot()
    })

    onActivated(()=>updateSpectraPlot())

    return { flexRowTitleStyle,
      materialPlots,
    isPlotReN, isPlotImN, isPlotInterpolation,
    plotRange}
  }
})
</script>
