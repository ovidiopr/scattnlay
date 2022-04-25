<template>
  <router-view />
</template>
<script lang="ts">
import { defineComponent, watch } from 'vue';
import { useStore } from 'src/store';

export default defineComponent({
  name: 'App',
  setup() {
    const $store = useStore();
    void (async () => $store.dispatch('simulationSetup/loadScattnlay'))();
    void $store.dispatch('guiRuntime/activateMaterial', 'main/Ag/McPeak.yml');
    void $store.dispatch('guiRuntime/activateMaterial', 'main/Au/McPeak.yml');
    void $store.dispatch('guiRuntime/activateMaterial', 'main/Al/McPeak.yml');
    void $store.dispatch('guiRuntime/activateMaterial', 'main/Cu/McPeak.yml');
    void $store.dispatch(
      'guiRuntime/activateMaterial',
      'main/Si/Green-2008.yml'
    );
    void $store.dispatch('guiRuntime/activateMaterial', 'main/SiO2/Gao.yml');
    void $store.dispatch('guiRuntime/activateMaterial', 'main/TiO2/Sarkar.yml');

    let isPlotInitialToggle = false;
    watch($store.state.guiRuntime.activatedMaterials, () => {
      if (!isPlotInitialToggle) {
        const indexToToggle =
          $store.state.guiRuntime.activatedMaterials.findIndex(
            (val) => val.name == 'Ag_McPeak'
          );
        // Materials are activated in async actions, so toggle Ag_McPeak to be plotted as soon as it is loaded.
        if (indexToToggle != -1) {
          $store.commit(
            'guiRuntime/toggleIsPlot',
            $store.state.guiRuntime.activatedMaterials[indexToToggle].name
          );
          isPlotInitialToggle = true;
        }
      }
    });
  },
});
</script>
