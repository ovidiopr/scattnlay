<template>
  <q-btn
    color="primary"
    no-caps
    :disable="data ? false : true"
    @click="saveSpectrumSimulation"
    >Save</q-btn
  >
  <q-checkbox v-model="isSaveWithPythonScript" size="sm"
    >Include *.py plotting script</q-checkbox
  >
</template>

<script lang="ts">
import {
  computed,
  defineComponent,
  // ref,
  // watch
  // nextTick
} from 'vue';
import { useStore } from 'src/store';
import { saveAs } from 'file-saver';

export default defineComponent({
  name: 'SaveSimulationNearField',

  setup() {
    const $store = useStore();

    const isSaveWithPythonScript = computed({
      get: () => $store.state.guiRuntime.isSaveWithPythonScript,
      set: (val) => $store.commit('guiRuntime/setIsSaveWithPythonScript', val),
    });

    const coordX = computed(() => $store.state.plotRuntime.nearFieldCoordX);
    const coordY = computed(() => $store.state.plotRuntime.nearFieldCoordY);
    const data = computed(() => $store.state.plotRuntime.nearFieldEabs);
    const layerWidths = computed(() =>
      $store.state.simulationSetup.current.layers.map((x) => x.layerWidth)
    );
    const crossSection = computed(() => {
      let val =
        $store.state.simulationSetup.current.nearFieldSetup.crossSection;
      if (val == 0) return 'Ek';
      if (val == 1) return 'Hk';
      return 'EH';
    });
    const radii = computed(() => {
      let r = 0;
      let radiiLocal = '';
      for (let widthLayer of layerWidths.value) {
        r += widthLayer;
        radiiLocal += r.toString() + ',';
      }
      return radiiLocal;
    });
    const lengthUnits = computed(() => $store.state.guiRuntime.units);
    return {
      data,
      isSaveWithPythonScript,
      saveSpectrumSimulation() {
        if (!data.value) return;
        const pythonScript =
          '# You can open and plot this file using Python\n' +
          '# (without manually removing this header, it will be skipped), see example below.\n' +
          'import numpy as np\n' +
          'from matplotlib import pyplot as plt\n' +
          "data = np.genfromtxt('scattnlay-near-field.txt', names=True, delimiter=', ')\n" +
          'x = data[data.dtype.names[0]] # x-axis has units\n' +
          'y = data[data.dtype.names[1]] # x-axis has units\n' +
          'x_size = len(np.unique(x))\n' +
          'y_size = len(np.unique(y))\n' +
          'x_min = np.min(x)\n' +
          'y_min = np.min(y)\n' +
          'x_max = np.max(x)\n' +
          'y_max = np.max(y)\n' +
          'dx = (x_max-x_min)/(x_size-1)\n' +
          'dy = (y_max-y_min)/(y_size-1)\n' +
          "Eabs = data['Eabs'].reshape((y_size,x_size))\n" +
          '\n' +
          'plt.figure()\n' +
          "plt.imshow(Eabs, cmap='jet', origin='lower',\n" +
          '           extent=(x_min-dx/2, x_max+dx/2, y_min-dy/2, y_max+dy/2))\n' +
          'cbar = plt.colorbar(shrink=0.7)\n' +
          "cbar.ax.set_title(r'$|E|/|E_0|$')\n" +
          "plt.xlabel(data.dtype.names[0].replace('_',', '))\n" +
          "plt.ylabel(data.dtype.names[1].replace('_',', '))\n" +
          'fig = plt.gcf()\n' +
          'ax = fig.gca()\n' +
          'radii = [' +
          radii.value +
          ']\n' +
          'for r in radii:\n' +
          "    ax.add_patch(plt.Circle((0, 0), r, color='w', alpha=0.7, lw=3, fill=False))\n" +
          "plt.title('" +
          crossSection.value +
          " cross-section')\n" +
          'plt.show()\n\n';

        let columnNames =
          '# X [' +
          lengthUnits.value.toString() +
          '] , Y [' +
          lengthUnits.value.toString() +
          '], Eabs\n';
        let body = '';
        for (let i = 0; i < data.value.length; ++i) {
          let row =
            coordX.value[i].toString() +
            ', ' +
            coordY.value[i].toString() +
            ', ' +
            data.value[i].toString();
          row += '\n';
          body += row;
        }

        const scattnlayNearField = new Blob(
          [columnNames + body],
          { type: 'text/plain;charset=utf-8', endings: 'native' } //TODO test if newline is correctly written in Windows, MacOS
        );
        saveAs(scattnlayNearField, 'scattnlay-near-field.txt');

        const scattnlayNearFieldScript = new Blob(
          [pythonScript],
          { type: 'text/plain;charset=utf-8', endings: 'native' } //TODO test if newline is correctly written in Windows, MacOS
        );
        saveAs(scattnlayNearFieldScript, 'scattnlay-near-field-plot.py');
      },
    };
  },
});
</script>
