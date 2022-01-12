<template>
  <div class="row items-baseline">
    <div class="col-xs-12 col-sm-auto text-weight-bold text-center q-pr-md q-py-sm">
      <div :style="flexRowTitleStyle">
        Plot
      </div>
    </div>
    <div class="col-xs-grow col-sm">
      <div class="row justify-center items-baseline">
        <div class="col-auto text-center q-py-xs q-pr-md">
          <div :style="flexRowTitleStyle">
            <span  class="text-italic">2<sup>&thinsp;nd</sup></span> side resolution
          </div>
        </div>
        <div class="col-xs-grow col-sm">
          <div class="row justify-xs-center justify-sm-start items-center">
            <input-with-units
                v-model:input-result="plotYSideResolution"
                v-model:is-showing-help="isShowingHelpForInputWithUnits"
                :initial-expression="plotYSideResolution.toString()"
                title="points"
                units=""
            /></div>
        </div>
      </div>
      <div class="row justify-center items-baseline">
        <div class="col-auto text-center q-py-xs q-pr-md">
          <div :style="flexRowTitleStyle" >
        relative center
          </div>
        </div>
        <div class="col-xs-grow col-sm">
          <div class="row justify-xs-center justify-sm-start items-center">
            <div class="col-auto">
              <q-tooltip anchor="center start" self="center end" >
                E-field vector direction
              </q-tooltip><input-with-units
                  v-model:input-result="atRelativeX0"
                  v-model:is-showing-help="isShowingHelpForInputWithUnits"
                  :initial-expression="atRelativeX0.toString()"
                  title="𝑋𝑜&thinsp;/&hairsp;𝑅"
                  units=""
              />
            </div>
            <div class="col-auto">
              <q-tooltip anchor="center start" self="center end" >
                H-field vector direction
              </q-tooltip><input-with-units
                v-model:input-result="atRelativeY0"
                v-model:is-showing-help="isShowingHelpForInputWithUnits"
                :initial-expression="atRelativeY0.toString()"
                title="𝑌𝑜&thinsp;/&hairsp;𝑅"
                units=""
            />
            </div>
            <div class="col-auto">
              <q-tooltip anchor="center start" self="center end" >
                k-vector direction
              </q-tooltip><input-with-units
                v-model:input-result="atRelativeZ0"
                v-model:is-showing-help="isShowingHelpForInputWithUnits"
                :initial-expression="atRelativeZ0.toString()"
                title="𝑍𝑜&thinsp;/&hairsp;𝑅"
                units=""
            />
            </div>
          </div>
        </div>
      </div>
      <div class="q-ma-xs"/>

    </div>
  </div>
</template>

<script lang="ts">
import {
  defineComponent,
  computed,
  // watch
  } from 'vue'
import { useStore } from 'src/store'
import InputWithUnits from 'components/InputWithUnits.vue'
import { flexRowTitleStyle } from 'components/config'

export default defineComponent({

  name: 'GetNearFieldRefinedSettings',
  components: {InputWithUnits,},

  setup() {
    const $store = useStore()

    const isShowingHelpForInputWithUnits = computed({
      get: () => $store.state.guiRuntime.isShowingHelpForInputWithUnits,
      set: val => $store.commit('guiRuntime/setIsShowingHelpForInputWithUnits', val)
    })

    const atRelativeX0 = computed({
      get: () => $store.state.simulationSetup.gui.nearFieldSetup.atRelativeX0,
      set: val => $store.commit('simulationSetup/setNearFieldAtRelativeX0', val)
    })
    const atRelativeY0 = computed({
      get: () => $store.state.simulationSetup.gui.nearFieldSetup.atRelativeY0,
      set: val => $store.commit('simulationSetup/setNearFieldAtRelativeY0', val)
    })
    const atRelativeZ0 = computed({
      get: () => $store.state.simulationSetup.gui.nearFieldSetup.atRelativeZ0,
      set: val => $store.commit('simulationSetup/setNearFieldAtRelativeZ0', val)
    })

    const plotYSideResolution = computed({
      get: () => $store.state.simulationSetup.gui.nearFieldSetup.plotYSideResolution,
      set: val => $store.commit('simulationSetup/setNearFieldPlotYSideResolution', Math.floor(val))
    })

    return {isShowingHelpForInputWithUnits, flexRowTitleStyle,
      plotYSideResolution,
      atRelativeX0, atRelativeY0, atRelativeZ0}
  },
})
</script>
