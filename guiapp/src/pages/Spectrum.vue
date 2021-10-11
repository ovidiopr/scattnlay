<template>
  <q-page class="row items-center justify-evenly">
    <!--    Your code-->
    <h4>Spectrum</h4>
          Wavelength
          <input-with-units
              v-model:input-result="simulationSetupGui.fromWL"
              v-model:is-showing-help="isShowingHelpForInputWithUnits"
              :initial-expression="simulationSetupGui.fromWL.toString()"
              title="from"
              units="nm"
          />
          <input-with-units
              v-model:input-result="someValue"
              v-model:is-showing-help="isShowingHelpForInputWithUnits"
              :initial-expression="someExpr"
              title=""
              units=""
              active
          />
    Input result: {{$store.state.simulationSetup.gui.fromWL}}
<!--    tooltip_text="help text"-->
  </q-page>
</template>

<script lang='ts'>
import {
  defineComponent, ref,
  computed, watch, reactive
} from 'vue'
import { useStore } from 'src/store'
import InputWithUnits from 'components/InputWithUnits.vue'
import { cloneDeep } from 'lodash'


export default defineComponent({
  name: 'PageIndex',
  components: {InputWithUnits },
  setup() {
    const $store = useStore()
    let someValue = ref(10)
    let someExpr = ref('10')
    // InputWithUnits component will disable showing help after first input
    const isShowingHelpForInputWithUnits = computed({
      get: () => $store.state.plotRuntime.isShowingHelpForInputWithUnits,
      set: val => $store.commit('plotRuntime/setIsShowingHelpForInputWithUnits', val)
    })

    const simulationSetupGui = reactive(cloneDeep($store.state.simulationSetup.gui))
    watch(simulationSetupGui, () => $store.commit('simulationSetup/setGuiState',simulationSetupGui))

    // console.log($store.state.simulationSetup.gui.toWL)
    return { someValue, someExpr, isShowingHelpForInputWithUnits, simulationSetupGui}
  }
})
</script>
