<template>
  <q-page class="row items-center justify-evenly">
    <!--    Your code-->

<!--    <div class="row justify-start items-center">-->
    <div class="column q-px-md">

      <div class="row">
        <div class="col-xs-12 col-sm-auto text-weight-bold text-left q-px-md q-py-sm">Host media</div>
        <div class="col-xs-auto col-sm">
          <div class="row justify-start items-center">
            <div class="col-auto" ><input-with-units
              v-model:input-result="simulationSetupGui.hostIndex"
              v-model:is-showing-help="isShowingHelpForInputWithUnits"
              :initial-expression="simulationSetupGui.hostIndex.toString()"
              title="Re(n)"
              units=""
          /></div>
          </div>
        </div>
      </div>


      <div class="row">
        <div class="col-xs-12 col-sm-auto text-weight-bold text-left q-px-md q-py-sm">Wavelength</div>
        <div class="col-xs-auto col-sm">
          <div class="row justify-start items-center">
            <div class="col-auto"><input-with-units
                v-model:input-result="simulationSetupGui.fromWL"
                v-model:is-showing-help="isShowingHelpForInputWithUnits"
                :initial-expression="simulationSetupGui.fromWL.toString()"
                title="from"
                units="nm"
            /></div>
            <div class="col-auto"><input-with-units
                v-model:input-result="simulationSetupGui.toWL"
                v-model:is-showing-help="isShowingHelpForInputWithUnits"
                :initial-expression="simulationSetupGui.toWL.toString()"
                title="to"
                units="nm"
                active
            /></div>
            <div class="col-auto"><input-with-units
                v-model:input-result="simulationSetupGui.pointsWL"
                v-model:is-showing-help="isShowingHelpForInputWithUnits"
                :initial-expression="simulationSetupGui.pointsWL.toString()"
                title="points"
                units=""
                active
            /></div>
          </div>
        </div>
      </div>

      <div class="row justify-start items-center">
        <div class="col-auto">
          Input result: {{$store.state.simulationSetup.gui.fromWL}}
        </div>
      </div>
    </div>
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
