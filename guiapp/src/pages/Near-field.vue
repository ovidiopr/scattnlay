<template>
  <q-page class="column q-px-md">
    <div class="q-ma-md"/>
    <GetWlFromPlot/>
    <div class="col-auto">
      Input result: {{$store.state.simulationSetup.gui.nearFieldWL}}
    </div>
    <!-- place QPageScroller at end of page -->
    <q-page-scroller position="bottom-right" :scroll-offset="150" :offset="[18, 18]">
      <q-btn size="xs" fab icon="keyboard_arrow_up" color="primary" />
    </q-page-scroller>
  </q-page>
</template>

<script lang='ts'>
import {
  defineComponent
} from 'vue'
import GetWlFromPlot from 'components/GetWlFromPlot.vue'
// import { useStore } from 'src/store'


export default defineComponent({
  name: 'NearField',
  components: {GetWlFromPlot
  },
  setup() {

    // const $store = useStore()

    return {
    }
  }
})
</script>
