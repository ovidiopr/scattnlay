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

    <ReactiveChart :chart="$store.state.plotRuntime.materialPlots"/>
  </div>
</template>

<script>
import ReactiveChart from 'components/ReactiveChart.vue'
import { useStore } from 'src/store'
import {
  defineComponent,
  onActivated,
  ref,
  computed,
  watch
} from 'vue'
import { flexRowTitleStyle } from 'components/config'


export default defineComponent({
  name: 'PlotMaterials',
  components: {
    ReactiveChart,
  },
  setup () {
    const $store = useStore()

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


    const plotRange=ref('material data')
    const isPlotReN = ref(true)
    const isPlotImN = ref(true)
    const isPlotInterpolation = ref(true)

    const fromWL = computed(()=>$store.state.simulationSetup.gui.fromWL)
    const toWL = computed(()=>$store.state.simulationSetup.gui.toWL)
    const pointsWL = computed(()=>$store.state.simulationSetup.gui.pointsWL)
    function updateSpectraPlot() {
      setPlotTitle()

      $store.commit('plotRuntime/updateMaterialPlots',
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
    isPlotReN, isPlotImN, isPlotInterpolation,
    plotRange}
  }
})
</script>
