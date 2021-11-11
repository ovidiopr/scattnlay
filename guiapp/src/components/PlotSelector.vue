<template>
  <div class="row justify-xs-center justify-sm-start items-baseline">
    <div class="col-auto">
    </div>
    <q-toggle v-model="isQscaToggle" dense class="q-ml-md">Qsca</q-toggle>
    <q-toggle v-model="isQabsToggle" dense class="q-ml-md">Qabs</q-toggle>
    <q-toggle v-model="isQextToggle" dense class="q-ml-md">Qext</q-toggle>

  </div>
</template>

<script lang="ts">
import {
  defineComponent,
  computed,
  } from 'vue'
import { useStore } from 'src/store'
import { flexRowTitleStyle } from 'components/utils'

export default defineComponent({
  name: 'GetPlotSettings',

  setup() {
    const $store = useStore()
    const isQscaToggle = computed({
      get: () => $store.state.plotRuntime.isQscaToggle,
      set: val => $store.commit('plotRuntime/setQscaToggle', val)
    })
    const isQabsToggle = computed({
      get: () => $store.state.plotRuntime.isQabsToggle,
      set: val => $store.commit('plotRuntime/setQabsToggle', val)
    })
    const isQextToggle = computed({
      get: () => $store.state.plotRuntime.isQextToggle,
      set: val => $store.commit('plotRuntime/setQextToggle', val)
    })

    return { isQscaToggle, isQabsToggle, isQextToggle
      }
  },
})
</script>
