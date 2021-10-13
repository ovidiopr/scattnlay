<template>
  <div class="row items-baseline">
    <div class="col-xs-12 col-sm-auto text-weight-bold text-center q-px-md q-py-sm">
      Wavelength
    </div>
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
</template>

<script lang="ts">
import {
  defineComponent,
  reactive,
  computed,
  watch,
  onBeforeUnmount,
  } from 'vue'
import { useStore } from 'src/store'
import { cloneDeep } from 'lodash'
import InputWithUnits from 'components/InputWithUnits.vue'

export default defineComponent({

  name: 'GetSourceParameters',
  components: {InputWithUnits,},

  setup() {
    const $store = useStore()

    let simulationSetupGui = reactive(cloneDeep($store.state.simulationSetup.gui))
    const unsubscribe = $store.subscribe((mutation, /*state*/) => {
      if (mutation.type === 'simulationSetup/setGuiState') {
        let key: keyof typeof simulationSetupGui
        for (key in simulationSetupGui) {
          simulationSetupGui[key] = $store.state.simulationSetup.gui[key]
        }
      }
    })
    onBeforeUnmount(()=>unsubscribe())
    watch(simulationSetupGui, () => $store.commit('simulationSetup/setGuiState',simulationSetupGui))

    const isShowingHelpForInputWithUnits = computed({
      get: () => $store.state.plotRuntime.isShowingHelpForInputWithUnits,
      set: val => $store.commit('plotRuntime/setIsShowingHelpForInputWithUnits', val)
    })
    return { simulationSetupGui, isShowingHelpForInputWithUnits }
  },
})
</script>
