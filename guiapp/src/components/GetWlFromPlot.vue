<template>
  <div>
    <ReactiveChart
        :chart="$store.state.plotRuntime.spectrumPlots"
                   @settingID="mangeID($event)"/>
  </div>
</template>

<script lang="ts">
import ReactiveChart from 'components/ReactiveChart.vue'
import { useStore } from 'src/store'
import {
  defineComponent,
  computed,
  // watch,
    ref,
} from 'vue'
import { fromUnits, toUnits } from 'components/utils'


export default defineComponent({
  name: 'GetWlFromPlot',
  components: {
    ReactiveChart,
  },
  setup () {
    const $store = useStore()
    const chartID = ref('')

    const sourceUnits = computed( ()=>$store.state.guiRuntime.sourceUnits)

    const currentWavelengthInSourceUnits = computed({
      get: () => toUnits($store.state.simulationSetup.gui.nearFieldWL, sourceUnits.value),
      set: val => $store.commit('simulationSetup/setFarFieldWL', fromUnits(sourceUnits.value, val))
    })

    function mangeID(val:string) {
      chartID.value = val
      console.log(val)
    }
    return {currentWavelengthInSourceUnits,
      mangeID, chartID

    }
  }
})
</script>
