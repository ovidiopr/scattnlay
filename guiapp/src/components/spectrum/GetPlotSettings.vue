<template>
  <div class="row items-baseline">
    <div class="col-xs-12 col-sm-auto text-center q-px-md q-py-sm">
      <div :style="flexRowTitleStyle">
        <span class="text-weight-bold">Plot label</span> (optional)
      </div>
    </div>
    <div class="col-xs-grow col-sm q-px-xs">
      <div class="row justify-xs-center justify-sm-start items-center">
        <div class="col-auto">
          <q-input
            v-model="plotLabel"
            :shadow-text="
              plotLabel == '' ? shadowText + ' (Tab to complete)' : ''
            "
            outlined
            dense
            class="q-py-xs"
            :style="'width: ' + basicWidthStyle"
            @keydown="processInputFill"
            @focus="processInputFill"
          />
        </div>
      </div>
    </div>
  </div>
</template>

<script lang="ts">
import { defineComponent, computed } from 'vue';
import { event } from 'quasar';
const { stopAndPrevent } = event;

import { useStore } from 'src/store';
import {
  flexRowTitleStyle,
  basicSelectorWidthStyle,
  basicWidthStyle,
} from 'components/config';

export default defineComponent({
  name: 'GetPlotSettings',

  setup() {
    const $store = useStore();

    const plotLabel = computed({
      get: () => $store.state.simulationSetup.gui.plotLabel,
      set: (val) => $store.commit('simulationSetup/setPlotLabel', val),
    });

    const shadowText = computed(() => {
      const numberOfLayers = $store.state.simulationSetup.gui.layers.length;
      let particleType = 'bulk';
      if (numberOfLayers == 2) return 'core-shell';
      if (numberOfLayers > 2) return 'multilayer';
      const R =
        $store.state.simulationSetup.gui.layers[0].layerWidth.toString();
      const units = $store.state.guiRuntime.units;
      return particleType + ' R=' + R + units;
    });

    return {
      flexRowTitleStyle,
      basicSelectorWidthStyle,
      basicWidthStyle,
      plotLabel,
      shadowText,

      processInputFill(e: KeyboardEvent) {
        if (e === void 0) {
          return;
        }
        if (e.code === 'Tab') {
          if (shadowText.value.length > 0 && plotLabel.value == '') {
            stopAndPrevent(e);
            plotLabel.value = shadowText.value;
          }
        }
      },
    };
  },
});
</script>
