<template>
  <span>
  <q-tooltip
      v-if=" spectrumRangeStart > fromWavelengthStore ||
             spectrumRangeEnd   <   toWavelengthStore"
      anchor="top middle" self="bottom middle"
      class="bg-red">
    Mismatch with spectrum simulation
  </q-tooltip>
  <span :class="spectrumRangeStart > fromWavelengthStore?'text-red':'text-black'">
              {{ Math.ceil(spectrumRangeStart) }}
            </span>
  &ndash;
  <span :class="spectrumRangeEnd < toWavelengthStore?'text-red':'text-black'">
              {{ Math.floor(spectrumRangeEnd) }}
            </span>
  &NonBreakingSpace;nm
    </span>
</template>

<script lang="ts">
import {
  computed,
  defineComponent,
} from 'vue'
import { useStore } from 'src/store'

export default defineComponent({
  name: 'ShowSpectrumRange',
  props: {
    spectrumRangeStart: {
      type: Number,
      required: true,
    },
    spectrumRangeEnd: {
      type: Number,
      required: true,
    },
  },
  setup() {
    const $store = useStore()
    const fromWavelengthStore = computed(()=>$store.state.simulationSetup.gui.fromWL)
    const toWavelengthStore = computed(()=>$store.state.simulationSetup.gui.toWL)

    return {
      fromWavelengthStore, toWavelengthStore
    }
  }
})
</script>
