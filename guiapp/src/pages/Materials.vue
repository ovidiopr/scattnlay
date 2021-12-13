<template>
  <q-page class="col-12 q-px-md">
    <div class="q-ma-sm"/>
    <MaterialsSelector/>
    <div class="q-ma-xs"/>
    <router-link to="/spectrum#GetSourceParametersHref"> <GetSourceParameters :is-info-mode="true"/> </router-link>
    <MaterialsActivated/>
    <div class="q-ma-xs"/>
  </q-page>
</template>

<script lang='ts'>
import {
  defineComponent,
} from 'vue'

import MaterialsSelector from 'components/MaterialsSelector.vue'
import MaterialsActivated from 'components/MaterialsActivated.vue'
import GetSourceParameters from 'components/GetSourceParameters.vue'
// import { useStore } from 'src/store'


export default defineComponent({
  name: 'PageIndex',
  components: {
    MaterialsSelector, GetSourceParameters,
    MaterialsActivated
  },
})
</script>
