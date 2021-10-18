<template>
  <q-page class="column q-px-md">
    <div class="q-ma-sm"/>
    <GetUnits/>
    <div class="q-ma-xs"/>
    <GetHostIndex/>
    <div class="q-ma-xs"/>
    <GetHostIndex/>
    <div class="q-ma-xs"/>
    <GetSourceParameters/>
    <div class="q-ma-xs"/>
    <GetSourceParameters/>
    <div class="col-auto">
      Input result: {{$store.state.simulationSetup.gui.fromWL}}
    </div>
  </q-page>
</template>

<script lang='ts'>
import {
  defineComponent
} from 'vue'
import GetUnits from 'components/GetUnits.vue'
import GetHostIndex from 'components/GetHostIndex.vue'
import GetSourceParameters from 'components/GetSourceParameters.vue'
// import { useStore } from 'src/store'


export default defineComponent({
  name: 'PageIndex',
  components: {GetUnits, GetHostIndex, GetSourceParameters },
  setup() {
    // const $store = useStore()
    return {}
  }
})
</script>
