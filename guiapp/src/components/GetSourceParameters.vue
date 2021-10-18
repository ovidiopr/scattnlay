<template>
  <div class="row items-baseline">
    <div class="col-xs-12 col-sm-auto text-weight-bold text-center q-px-md q-py-sm">
      <div :style="flexRowTitleStyle">
        Wavelength
      </div>
    </div>
    <div class="col-xs-grow col-sm">
      <div class="row justify-xs-center justify-sm-start items-center">

        <div class="col-auto"><input-with-units
            v-model:input-result="fromWL"
            v-model:is-showing-help="isShowingHelpForInputWithUnits"
            :initial-expression="fromWL.toString()"
            title="from"
            units="nm"
        /></div>
        <div class="col-auto"><input-with-units
            v-model:input-result="toWL"
            v-model:is-showing-help="isShowingHelpForInputWithUnits"
            :initial-expression="toWL.toString()"
            title="to"
            units="nm"
            active
        /></div>
        <div class="col-auto"><input-with-units
            v-model:input-result="pointsWL"
            v-model:is-showing-help="isShowingHelpForInputWithUnits"
            :initial-expression="pointsWL.toString()"
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
  computed,
  } from 'vue'
import { useStore } from 'src/store'
import InputWithUnits from 'components/InputWithUnits.vue'
import { flexRowTitleStyle } from 'components/utils'

export default defineComponent({

  name: 'GetSourceParameters',
  components: {InputWithUnits,},

  setup() {
    const $store = useStore()

    const isShowingHelpForInputWithUnits = computed({
      get: () => $store.state.plotRuntime.isShowingHelpForInputWithUnits,
      set: val => $store.commit('plotRuntime/setIsShowingHelpForInputWithUnits', val)
    })

    const fromWL = computed({
      get: () => $store.state.simulationSetup.gui.fromWL,
      set: val => $store.commit('simulationSetup/setFromWL', val)
    })

    const toWL = computed({
      get: () => $store.state.simulationSetup.gui.toWL,
      set: val => $store.commit('simulationSetup/setToWL', val)
    })

    const pointsWL = computed({
      get: () => $store.state.simulationSetup.gui.pointsWL,
      set: val => $store.commit('simulationSetup/setPointsWL', val)
    })

    return { fromWL, toWL, pointsWL, isShowingHelpForInputWithUnits, flexRowTitleStyle }
  },
})
</script>
