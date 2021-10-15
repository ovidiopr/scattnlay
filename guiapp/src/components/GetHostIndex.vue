<template>
  <div class="row items-baseline">
    <div class="col-xs-12 col-sm-auto text-weight-bold text-center q-px-md q-py-sm">
      Host media
    </div>
    <div class="col-xs-auto col-sm">
      <div class="row justify-start items-center">

        <div class="col-auto" ><input-with-units
            v-model:input-result="hostIndex"
            v-model:is-showing-help="isShowingHelpForInputWithUnits"
            :initial-expression="hostIndex.toString()"
            title="Re(n)"
            units=""
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

export default defineComponent({

  name: 'GetHostIndex',
  components: {InputWithUnits,},

  setup() {
    const $store = useStore()

    const isShowingHelpForInputWithUnits = computed({
      get: () => $store.state.plotRuntime.isShowingHelpForInputWithUnits,
      set: val => $store.commit('plotRuntime/setIsShowingHelpForInputWithUnits', val)
    })

    const hostIndex = computed({
      get: () => $store.state.simulationSetup.gui.hostIndex,
      set: val => $store.commit('simulationSetup/setHostIndex', val)
    })

    return { hostIndex, isShowingHelpForInputWithUnits }
  },
})
</script>
