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
          <div :style="flexRowTitleStyle" >
            max side resolution
          </div>
        </div>
        <div class="col-xs-grow col-sm">
          <div class="row justify-xs-center justify-sm-start items-center">
            <input-with-units
                v-model:input-result="plotSideResolution"
                v-model:is-showing-help="isShowingHelpForInputWithUnits"
                :initial-expression="plotSideResolution.toString()"
                title="points"
                units=""
            /></div>
        </div>
      </div>
      <div class="row justify-center items-baseline">
        <div class="col-auto text-center q-py-xs q-pr-md">
          <div :style="flexRowTitleStyle" >
        relative side length
          </div>
        </div>
        <div class="col-xs-grow col-sm">
          <div class="row justify-xs-center justify-sm-start items-center">
            <input-with-units
            v-model:input-result="relativePlotSize"
            v-model:is-showing-help="isShowingHelpForInputWithUnits"
            :initial-expression="relativePlotSize.toString()"
            title="𝐿&thinsp;/&hairsp;𝟐𝑅"
            units=""
            /></div>
        </div>
      </div>
      <div class="q-ma-xs"/>
      <div class="row justify-center items-center">
        <div class="col-auto text-center q-py-xs q-pr-md">
          <div :style="flexRowTitleStyle" >
            cross-section
          </div>
        </div>
        <div class="col-xs-grow col-sm">
          <div class="row justify-xs-center justify-sm-start items-center">
            <div class="q-gutter-md q-py-sm q-px-xs">
              <q-radio v-model="crossSection" dense size='sm' :val="nearFieldPlane.all" label="All" />
              <q-radio v-model="crossSection" dense size='sm' :val="nearFieldPlane.Ek" label="Ek" />
              <q-radio v-model="crossSection" dense size='sm' :val="nearFieldPlane.Hk" label="Hk" />
              <q-radio v-model="crossSection" dense size='sm' :val="nearFieldPlane.EH" label="EH" />
            </div>
          </div>
        </div>
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
import { flexRowTitleStyle } from 'components/config'
import { nearFieldPlane } from 'src/store/simulation-setup/state';

export default defineComponent({

  name: 'GetNearFieldSettings',
  components: {InputWithUnits,},

  setup() {
    const $store = useStore()

    const isShowingHelpForInputWithUnits = computed({
      get: () => $store.state.guiRuntime.isShowingHelpForInputWithUnits,
      set: val => $store.commit('guiRuntime/setIsShowingHelpForInputWithUnits', val)
    })

    const crossSection = computed({
      get: () => $store.state.simulationSetup.gui.nearFieldSetup.crossSection,
      set: val => $store.commit('simulationSetup/setNearFieldCrossSection', val)
    })

    const relativePlotSize = computed({
      get: () => $store.state.simulationSetup.gui.nearFieldSetup.relativePlotSize,
      set: val => $store.commit('simulationSetup/setNearFieldRelativePlotSize', val)
    })

    const maxComputeTime = computed({
      get: () => $store.state.simulationSetup.gui.nearFieldSetup.maxComputeTime,
      set: val => $store.commit('simulationSetup/setNearFieldMaxComputeTime', val)
    })

    const plotSideResolution = computed({
      get: () => $store.state.simulationSetup.gui.nearFieldSetup.plotSideResolution,
      // TODO: make InputWithUnits to handle integer input, so no need to use floor() in the next line.
      set: val => $store.commit('simulationSetup/setNearFieldPlotSideResolution', Math.floor(val))
    })

    return { crossSection, isShowingHelpForInputWithUnits, flexRowTitleStyle,
    relativePlotSize, maxComputeTime, plotSideResolution, nearFieldPlane}
  },
})
</script>
