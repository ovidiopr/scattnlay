<template>
  <q-page class="column q-px-md">
    <div class="q-ma-md"/>
    <router-link to="/spectrum#GetSourceParametersHref"> <GetSourceParameters :is-info-mode="true"/> </router-link>
    <div class="q-ma-xs"/>
    <MaterialsSelector/>
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
