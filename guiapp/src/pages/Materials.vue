<template>
  <q-page class="column q-px-md">
    <div class="q-ma-md"/>
    <MaterialsSelector/>
  </q-page>
</template>

<script lang='ts'>
import {
  defineComponent,
} from 'vue'

import MaterialsSelector from 'components/MaterialsSelector.vue'
// import { useStore } from 'src/store'


export default defineComponent({
  name: 'PageIndex',
  components: {
    MaterialsSelector
  },
})
</script>
