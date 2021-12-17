<template>
  <div class="row items-baseline">

    <div class="col-xs-12 col-sm-auto text-weight-bold text-center q-px-md q-py-sm">
      <q-tooltip v-if="isInfoMode"
                 anchor="top middle"
                 self="center middle"
      >
        current  source settings<br>
      </q-tooltip>
      <div :style="flexRowTitleStyle"> {{rowTitle}} </div>
    </div>
    <div class="col-xs-grow col-sm">
      <div class="row justify-xs-center justify-sm-start items-baseline">

        <div class="col-auto"><input-with-units
            v-model:input-result="fromSource"
            v-model:is-showing-help="isShowingHelpForInputWithUnits"
            :initial-expression="fromSource.toString()"
            :units="sourceUnits"
            :is-info-mode="isInfoMode"
            :is-error="fromSource<$store.state.guiRuntime.safeFromWL"
            title="from"
        /></div>
        <div class="col-auto"><input-with-units
            v-model:input-result="toSource"
            v-model:is-showing-help="isShowingHelpForInputWithUnits"
            :initial-expression="toSource.toString()"
            :units="sourceUnits"
            :is-info-mode="isInfoMode"
            :is-error="toSource>$store.state.guiRuntime.safeToWL"
            title="to"
        /></div>
        <div v-if="!isInfoMode" class="col-auto">
          <input-with-units
            v-model:input-result="pointsSource"
            v-model:is-showing-help="isShowingHelpForInputWithUnits"
            :initial-expression="pointsSource.toString()"
            :title="isPointsToggle ? `points` : `step`"
            :units="isPointsToggle ? `` : sourceUnits"
          />
          <q-toggle v-model="isPointsToggle" dense class="q-ml-md"/>
        </div>

      </div>
    </div>
  </div>
</template>

<script lang="ts">
import {
  defineComponent,
  computed,
  ref
  } from 'vue'
import { useStore } from 'src/store'
import InputWithUnits from 'components/InputWithUnits.vue'
import { fromUnits, toUnits, isAlmostSame } from 'components/utils'
import { flexRowTitleStyle } from 'components/config'

export default defineComponent({

  name: 'GetSourceParameters',
  components: {InputWithUnits,},
  props: {
    isInfoMode: {
      type: Boolean,
      default: false
    },
  },

  setup() {
    const isPointsToggle = ref(true)

    const $store = useStore()

    const isShowingHelpForInputWithUnits = computed({
      get: () => $store.state.guiRuntime.isShowingHelpForInputWithUnits,
      set: val => $store.commit('guiRuntime/setIsShowingHelpForInputWithUnits', val)
    })

    const sourceUnits = computed(()=>{
      return $store.state.guiRuntime.sourceUnits
    })

    const fromWavelengthInSourceUnits = computed({
      get: () => toUnits($store.state.simulationSetup.gui.fromWL, sourceUnits.value),
      set: val => {
        if (!isAlmostSame($store.state.simulationSetup.gui.fromWL, fromUnits(sourceUnits.value, val))) {
          $store.commit('simulationSetup/setFromWL', fromUnits(sourceUnits.value, val))
        }
      }
    })

    const toWavelengthInSourceUnits = computed({
      get: () => toUnits($store.state.simulationSetup.gui.toWL, sourceUnits.value),
      set: val => {
        if (!isAlmostSame($store.state.simulationSetup.gui.toWL, fromUnits(sourceUnits.value, val))) {
          $store.commit('simulationSetup/setToWL', fromUnits(sourceUnits.value, val))
        }
      }
    })

    const rowTitle = computed(()=>{
      if (sourceUnits.value.endsWith('Hz')) {return 'Frequency'}
      if (sourceUnits.value.endsWith('eV')) {return 'Energy'}
      if (sourceUnits.value.endsWith('s')) {return 'Period'}
      return 'Wavelength'
    })

    const directOrder = computed(()=>{
      if (['Frequency','Energy'].includes(rowTitle.value)) return false
      return true
    })

    const fromSource = computed({
      get: () => {
        if (directOrder.value) return fromWavelengthInSourceUnits.value
        return toWavelengthInSourceUnits.value
      },
      set: val => {
        if (directOrder.value) fromWavelengthInSourceUnits.value = val
        else toWavelengthInSourceUnits.value = val
      }
    })

    const toSource = computed({
      get: () => {
        if (!directOrder.value) return fromWavelengthInSourceUnits.value
        return toWavelengthInSourceUnits.value
      },
      set: val => {
        if (!directOrder.value) fromWavelengthInSourceUnits.value = val
        else toWavelengthInSourceUnits.value = val
      }
    })

    const pointsSource = computed({
      get: () => {
        if (isPointsToggle.value) return $store.state.simulationSetup.gui.pointsWL
        return (toSource.value-fromSource.value)/($store.state.simulationSetup.gui.pointsWL-1)
      },
      set: val => {
        if (isPointsToggle.value) $store.commit('simulationSetup/setPointsWL', val)
        else {
          $store.commit('simulationSetup/setPointsWL', 1+(toSource.value-fromSource.value)/val)
        }
      }
    })

    return { fromSource, toSource, pointsSource, isShowingHelpForInputWithUnits,
      flexRowTitleStyle, rowTitle,
      sourceUnits, isPointsToggle,  }
  },
})
</script>
