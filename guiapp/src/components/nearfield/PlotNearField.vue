<template>
  <div>
    <ReactiveChart :chart="nearFieldPlot" />
  </div>
</template>

<script lang="ts">
import ReactiveChart from 'components/ReactiveChart.vue';
import { useStore } from 'src/store';
import {
  defineComponent,
  // ref,
  reactive,
  computed,
  watch,
} from 'vue';
// import { flexRowTitleStyle } from 'components/config'
import {
  toUnits,
  getMaxFromHeatmap,
  getMinFromHeatmap,
  limitMap,
} from 'components/utils';
import { plotlyChart } from 'src/store/plot-runtime/state';
import { PlotData, DataTitle } from 'plotly.js-dist-min';
import { nearFieldPlane } from 'src/store/simulation-setup/state';

export default defineComponent({
  name: 'PlotNearField',
  components: {
    ReactiveChart,
  },
  setup() {
    const $store = useStore();
    const nearFieldPlotInit: plotlyChart = {
      data: [],
      layout: {
        shapes: [],
        margin: {
          l: 0,
          r: 40,
          b: 50,
          t: 30,
        },
        xaxis: {
          title: '',
        },
        yaxis: {
          scaleanchor: 'x',
          title: '',
        },
        showlegend: false,
      },
      config: {
        responsive: true,
        // showEditInChartStudio: true,
        displaylogo: false,
      },
    };

    const nearFieldPlot = reactive(nearFieldPlotInit);
    const crossSection = computed(
      () => $store.state.simulationSetup.current.nearFieldSetup.crossSection
    );
    const relativePlotSize = computed(
      () => $store.state.simulationSetup.current.nearFieldSetup.relativePlotSize
    );
    const plotYSideResolution = computed(
      () =>
        $store.state.simulationSetup.current.nearFieldSetup.plotYSideResolution
    );
    const plotXSideResolution = computed(
      () =>
        $store.state.simulationSetup.current.nearFieldSetup.plotXSideResolution
    );
    const layerWidths = computed(() =>
      $store.state.simulationSetup.current.layers.map((x) => x.layerWidth)
    );
    const totalR = computed(() => layerWidths.value.reduce((a, b) => a + b, 0));
    const x0 = computed(() => totalR.value * relativePlotSize.value);
    const dx = computed(
      () => (x0.value * 2.0) / (plotXSideResolution.value - 1)
    );
    const units = computed(() => $store.state.guiRuntime.units);
    const at_x = computed(() => {
      if (
        crossSection.value == nearFieldPlane.Ek ||
        crossSection.value == nearFieldPlane.Hk
      ) {
        return $store.state.simulationSetup.current.nearFieldSetup.atRelativeZ0;
      }
      return $store.state.simulationSetup.current.nearFieldSetup.atRelativeY0;
    });

    const at_y = computed(() => {
      if (crossSection.value == nearFieldPlane.Hk) {
        return $store.state.simulationSetup.current.nearFieldSetup.atRelativeY0;
      }
      return $store.state.simulationSetup.current.nearFieldSetup.atRelativeX0;
    });

    const xy = computed(() => {
      let x: number[] = [];
      let y: number[] = [];
      const xi =
        at_x.value * totalR.value -
        (dx.value * (plotXSideResolution.value - 1)) / 2;
      const yi =
        at_y.value * totalR.value -
        (dx.value * (plotYSideResolution.value - 1)) / 2;

      for (let j = 0; j < plotYSideResolution.value; ++j) {
        for (let i = 0; i < plotXSideResolution.value; ++i) {
          x.push(toUnits(xi + i * dx.value, units.value));
          y.push(toUnits(yi + j * dx.value, units.value));
        }
      }
      $store.commit('plotRuntime/setNearFieldCoords', { x: x, y: y });
      return { x: x, y: y };
    });

    const limitFrom = computed(
      () => $store.state.plotRuntime.nearFieldLimitFrom
    );
    const limitTo = computed(() => $store.state.plotRuntime.nearFieldLimitTo);

    const nearFieldStore = computed(() => {
      let nearFieldStoreLocal = $store.state.plotRuntime.nearFieldEabs;
      $store.commit(
        'plotRuntime/setNearFieldDataTo',
        getMaxFromHeatmap(nearFieldStoreLocal)
      );
      $store.commit(
        'plotRuntime/setNearFieldDataFrom',
        getMinFromHeatmap(nearFieldStoreLocal)
      );
      $store.commit(
        'plotRuntime/setNearFieldLimitTo',
        $store.state.plotRuntime.nearFieldDataTo
      );
      $store.commit(
        'plotRuntime/setNearFieldLimitFrom',
        $store.state.plotRuntime.nearFieldDataFrom
      );
      return nearFieldStoreLocal;
    });
    const nearFieldProc = computed(() => {
      if (!nearFieldStore.value) return nearFieldStore.value;
      return limitMap(nearFieldStore.value, limitFrom.value, limitTo.value);
      // return nearFieldStore.map(x=>x>limitTo.value?limitTo.value:x)
    });
    watch([nearFieldProc, xy], () => {
      nearFieldPlot.data.length = 0;
      const heatMapSettings: Partial<PlotData> = {
        type: 'heatmap',
        colorscale: 'Jet',
        colorbar: { title: '|𝐸|∕|𝐸𝜊|' },
        z: nearFieldProc.value,
      };
      nearFieldPlot.data.push({ ...xy.value, ...heatMapSettings });

      if (nearFieldPlot.layout.shapes) {
        nearFieldPlot.layout.shapes.length = 0;
        let r = 0;
        for (let widthLayer of layerWidths.value) {
          r += widthLayer;
          nearFieldPlot.layout.shapes.push({
            type: 'circle',
            xref: 'x',
            yref: 'y',
            x0: toUnits(-r, units.value),
            y0: toUnits(-r, units.value),
            x1: toUnits(r, units.value),
            y1: toUnits(r, units.value),
            line: {
              color: 'rgba(235, 235, 235, 0.8)',
              width: 3,
            },
          });
        }
      }
    });

    const xaxisTitle = computed(() => {
      let title: string | Partial<DataTitle> = 'x [' + units.value + ']';
      return title;
    });
    if (nearFieldPlot.layout.xaxis)
      nearFieldPlot.layout.xaxis.title = xaxisTitle.value;
    watch(xaxisTitle, () => {
      if (nearFieldPlot.layout.xaxis)
        nearFieldPlot.layout.xaxis.title = xaxisTitle.value;
    });

    return {
      nearFieldPlot,
      totalR,
    };
  },
});
</script>
