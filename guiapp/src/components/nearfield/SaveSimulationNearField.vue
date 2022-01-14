<template>
  <q-btn
      color="primary"
      no-caps
      :disable="data ? false : true"
      @click="saveSpectrumSimulation"
  >Save</q-btn>
</template>

<script lang="ts">
import {
  computed,
  defineComponent,
    // ref,
  // watch
  // nextTick
} from 'vue'
import {useStore} from 'src/store'
// import {getModeName, range, rangeInt} from 'components/utils'
// import {cloneDeep} from 'lodash'
import {saveAs} from 'file-saver'
// import {nearFieldPlane} from 'src/store/simulation-setup/state';


export default defineComponent({
  name: 'SaveSimulationNearField',

  setup() {
    const $store = useStore()


    const coordX = computed(() => $store.state.plotRuntime.nearFieldCoordX)
    const coordY = computed(() => $store.state.plotRuntime.nearFieldCoordY)
    const data = computed(() => $store.state.plotRuntime.nearFieldEabs)

    return {
      data,
      saveSpectrumSimulation(){
        if (!data.value) return
        const fileHeader = '# # You can open and plot this file using Python\n' +
            '# # (without manually removing this header, it will be skipped), see example below.\n' +
            '# import numpy as np\n' +
            '# from matplotlib import pyplot as plt\n' +
            '# data = np.genfromtxt(\'scattnlay-spectra.txt\', skip_header=21, names=True, delimiter=\', \')\n' +
            '# x = data[data.dtype.names[0]] # x-axis has units\n' +
            '# # Possible labels for spectrum data: Qsca, Qabs, Qext,\n' +
            '# # Qsca_E_dipole, etc. (see last comment before spectra data)\n' +
            '# a = data[\'Qsca\']\n' +
            '# b = data[\'Qsca_E_dipole\']\n' +
            '# c = data[\'Qsca_H_dipole\']\n' +
            '# \n' +
            '# plt.figure()\n' +
            '# plt.plot(x, a, label=\'Qsca\')\n' +
            '# plt.plot(x, b, label=\'Qsca E dipole\')\n' +
            '# plt.plot(x, c, label=\'Qsca H dipole\')\n' +
            '# plt.legend()\n' +
            '# plt.xlabel(data.dtype.names[0].replace(\'_\', \', \'))\n' +
            '# plt.ylabel(\'Normalized cross-sections\')\n' +
            '# plt.show()\n\n'
        let columnNames = '# X, Y, Eabs\n'
        let body = ''
        for (let i = 0; i < data.value.length; ++i) {
          let row = coordX.value[i].toString() + ', '
              + coordY.value[i].toString() + ', '
              + data.value[i].toString()
          row += '\n'
          body += row
        }

        const scattnlayNearField = new Blob([fileHeader+columnNames+body
            ],
            {type: 'text/plain;charset=utf-8',
              endings: 'native'}  //TODO test if newline is correctly written in Windows, MacOS
        )
        saveAs(scattnlayNearField, 'scattnlay-near-field.txt');
      }
    }
  },
})
</script>
