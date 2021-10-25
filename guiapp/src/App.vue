<template>
  <router-view />
</template>
<script lang="ts">
import { defineComponent } from 'vue';
import { useStore } from 'src/store'

export default defineComponent({
  name: 'App',
  setup(){
    const $store = useStore()
    $store.dispatch('simulationSetup/loadScattnlay')

  }
})
</script>
