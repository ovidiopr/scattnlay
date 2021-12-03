<template>
  <q-page class="column q-px-md">
    Under construction!
    <div class="q-ma-sm"/>
    <GetSourceParameters/>
    <GetHostIndex/>
    <div class="col-auto">
      Input result: {{$store.state.simulationSetup.gui.fromWL}}
    </div>
  </q-page>
</template>

<script lang='ts'>
import {
  defineComponent
} from 'vue'
import GetHostIndex from 'components/GetHostIndex.vue'
import GetSourceParameters from 'components/GetSourceParameters.vue'
// import GetParticleParameters from 'components/GetParticleParameters.vue'
// import { useStore } from 'src/store'


export default defineComponent({
  name: 'PageIndex',
  components: {GetHostIndex, GetSourceParameters,
    // GetParticleParameters
  },
  setup() {
    // const $store = useStore()
    return {}
  }
})
</script>
