<template>
  <div>
    <div class="row items-baseline">
      <div
        class="col-xs-12 col-sm-auto text-weight-bold text-center q-pr-md q-py-sm"
      >
        <div :style="flexRowTitleStyle">Source plane wave</div>
      </div>
      <div class="col-xs-grow col-sm">
        <div class="row justify-xs-center justify-sm-start items-baseline">
          <div class="col-auto">
            <input-with-units
              v-model:input-result="currentWavelengthInSourceUnits"
              v-model:is-showing-help="isShowingHelpForInputWithUnits"
              :initial-expression="currentWavelengthInSourceUnits.toString()"
              :units="sourceUnits"
              title="at"
            />
          </div>
          <div class="col-auto q-px-sm">
            or <span class="text-bold">click on plot</span> below to select a
            data point
          </div>
        </div>
      </div>
    </div>
    <div class="q-ma-xs" />
    <ReactiveChart
      :chart="chartContent"
      :window-height-share="0.4"
      @plotCreated="manageID($event)"
    />
  </div>
</template>

<script lang="ts">
import ReactiveChart from 'components/ReactiveChart.vue';
import { useStore } from 'src/store';
import { defineComponent, computed } from 'vue';
import { fromUnits, toUnits } from 'components/utils';
import InputWithUnits from 'components/InputWithUnits.vue';
import { flexRowTitleStyle } from 'components/config';
import { PlotlyHTMLElement } from 'plotly.js-dist-min';
import { cloneDeep } from 'lodash';

export default defineComponent({
  name: 'GetWlFromPlot',
  components: {
    ReactiveChart,
    InputWithUnits,
  },
  setup() {
    const $store = useStore();

    const sourceUnits = computed(() => $store.state.guiRuntime.sourceUnits);
    const isShowingHelpForInputWithUnits = computed({
      get: () => $store.state.guiRuntime.isShowingHelpForInputWithUnits,
      set: (val) =>
        $store.commit('guiRuntime/setIsShowingHelpForInputWithUnits', val),
    });

    const currentWavelengthInSourceUnits = computed({
      get: () =>
        toUnits(
          $store.state.simulationSetup.gui.nearFieldSetup.atWL,
          sourceUnits.value
        ),
      set: (val) =>
        $store.commit(
          'simulationSetup/setNearFieldWL',
          fromUnits(sourceUnits.value, val)
        ),
    });

    const chartContent = computed(() => {
      let content = cloneDeep($store.state.plotRuntime.spectrumPlots);
      if (content.layout.yaxis) content.layout.yaxis.title = 'Cross-section';
      content.layout['shapes'] = [
        {
          type: 'line',
          x0: currentWavelengthInSourceUnits.value,
          y0: 0,
          x1: currentWavelengthInSourceUnits.value,
          yref: 'paper',
          y1: 1,
          line: {
            color: 'red',
            width: 1.5,
            dash: 'dot',
          },
        },
      ];
      return content;
    });
    function manageID(chartID: string) {
      const myPlot = document.getElementById(chartID) as PlotlyHTMLElement;
      myPlot.on('plotly_click', function (data) {
        for (let i = 0; i < data.points.length; i++) {
          let val = data.points[i].x;
          if (val)
            currentWavelengthInSourceUnits.value = parseFloat(val.toString());
        }
      });
    }
    return {
      currentWavelengthInSourceUnits,
      sourceUnits,
      chartContent,
      flexRowTitleStyle,
      isShowingHelpForInputWithUnits,
      manageID,
    };
  },
});
</script>
