<template>
  <q-page class="col-12 q-px-md">
    <div class="q-ma-sm" />
    <MaterialsSelector />
    <div class="q-ma-xs" />
    <router-link to="/spectrum">
      <GetSourceParameters :is-info-mode="true" />
    </router-link>
    <MaterialsActivated />
    <div class="q-ma-xs" />
    <PlotMaterials />
    <div class="q-ma-xs" />
    <!-- place QPageScroller at end of page -->
    <div class="q-ma-lg" />
    <q-page-scroller
      position="bottom-right"
      :scroll-offset="150"
      :offset="[18, 18]"
    >
      <q-btn size="xs" fab icon="keyboard_arrow_up" color="primary" />
    </q-page-scroller>
  </q-page>
</template>

<script lang="ts">
import { defineComponent } from 'vue';

import MaterialsSelector from 'components/materials/MaterialsSelector.vue';
import MaterialsActivated from 'components/materials/MaterialsActivated.vue';
import GetSourceParameters from 'components/spectrum/GetSourceParameters.vue';
import PlotMaterials from 'components/materials/PlotMaterials.vue';
// import { useStore } from 'src/store'

export default defineComponent({
  name: 'MaterialsPage',
  components: {
    MaterialsSelector,
    GetSourceParameters,
    MaterialsActivated,
    PlotMaterials,
  },
});
</script>
