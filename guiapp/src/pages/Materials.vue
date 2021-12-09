<template>
  <q-page class="column q-px-md">
    <div class="q-ma-sm"/>
    <MaterialsSelector/>
    <div class="q-ma-xs"/>
    <router-link to="/spectrum#GetSourceParametersHref"> <GetSourceParameters :is-info-mode="true"/> </router-link>
  </q-page>
</template>

<script lang='ts'>
import {
  defineComponent,
} from 'vue'

import MaterialsSelector from 'components/MaterialsSelector.vue'
import GetSourceParameters from 'components/GetSourceParameters.vue'
// import { useStore } from 'src/store'


export default defineComponent({
  name: 'PageIndex',
  components: {
    MaterialsSelector, GetSourceParameters
  },
})
</script>
