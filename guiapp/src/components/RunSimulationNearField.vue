<template>
  <div class="row items-baseline">
    <div class="col-xs-12 col-sm-auto text-weight-bold text-center q-px-md q-py-sm">
      <q-tooltip
          v-if=" $store.state.guiRuntime.safeFromWL > $store.state.simulationSetup.gui.nearFieldSetup.atWL ||
                 $store.state.guiRuntime.safeToWL < $store.state.simulationSetup.gui.nearFieldSetup.atWL "
          anchor="top middle" self="center middle"
          class="bg-amber-4 text-black shadow-4">
        Will use materials<br> spectrum range.
      </q-tooltip>
        <q-btn :loading="isRunning"
               :disable="isRunning||!isNmieLoaded"
               color="primary"
               no-caps
               :label="isNmieLoaded ? 'Run simulation' : 'Loading...'"
               @click="runNearFieldSimulation">
          <template #loading>
            <q-spinner-gears />
          </template>
        </q-btn>
    </div>
    <div class="col-xs-grow col-sm q-px-xs">
      <div class="row justify-xs-center justify-sm-start items-baseline">

        <div class="col-auto">
          <q-btn
              color="primary"
              no-caps
              @click="saveSpectrumSimulation"
          >Save</q-btn>
        </div>

      </div>
    </div>
  </div>
</template>

<script lang="ts">
import {
  defineComponent,
  computed,
  // watch,
    onActivated,
    nextTick
} from 'vue'
import { useStore } from 'src/store'
import { getModeName, range, rangeInt} from 'components/utils'
import { cloneDeep } from 'lodash'
import { saveAs } from 'file-saver'

export default defineComponent({
  name: 'RunSimulationNearField',

  setup() {
    const $store = useStore()

    const isRunning = computed({
      get: ()=> $store.state.simulationSetup.nmies.nearField.isNmieRunning,
      set: val => {
        val ? $store.commit('simulationSetup/markNmieNearFieldAsStarted') : $store.commit('simulationSetup/markNmieNearFieldAsFinished')
      }
    })

    const isNmieLoaded = computed(()=>{ return $store.state.simulationSetup.nmies.nearField.instance })

    const atWL = computed(
        () => $store.state.simulationSetup.current.nearFieldSetup.atWL)
    //-------------------------------------------------------------------------//
    //---------------------  Main  --------------------------------------------//
    //-------------------------------------------------------------------------//
    function runNearFieldSimulation() {
      if (isRunning.value) {
        console.log('Some Nmie is already running!')
        return
      }
      isRunning.value = true
      void nextTick(()=> {
        $store.commit('simulationSetup/copySetupFromGuiToCurrent')

        const host = $store.state.simulationSetup.current.hostIndex

        try {
          if (!$store.state.simulationSetup.nmies.nearField.instance) throw 'ERROR! Scattnlay module was not loaded'
          const nmie = $store.state.simulationSetup.nmies.nearField.instance
          const layers = cloneDeep($store.state.simulationSetup.current.layers)
          const nmieStartedTime = performance.now()

          nmie.SetWavelength(atWL.value)
          nmie.ClearTarget()
          for (const layer of layers) {
            if (layer.material.nSpline) layer.n = layer.material.nSpline.at(atWL.value)
            if (layer.material.kSpline) layer.k = layer.material.kSpline.at(atWL.value)
            nmie.AddTargetLayerReIm(layer.layerWidth * host, layer.n / host, layer.k / host)
          }

          nmie.SetModeNmaxAndType(-1, -1)

          nmie.RunFieldCalculationCartesian(
              $store.state.simulationSetup.current.nearFieldSetup.plotSideResolution,
              $store.state.simulationSetup.current.nearFieldSetup.relativePlotSize,
              $store.state.simulationSetup.current.nearFieldSetup.crossSection,
              0, 0, 0, 0
          )


          const nmieTotalRunTime = (performance.now()-nmieStartedTime)/1000
          // console.log('Total simulation time:', nmieTotalRunTime, 's')
          $store.commit('simulationSetup/setNmieNearFieldTotalRunTime', nmieTotalRunTime)

        } catch (e) {
          console.log('Some error:', e)
        }
        isRunning.value = false
      })
    }

    onActivated(()=>{
      if (isNmieLoaded.value) runNearFieldSimulation()
    })



    return { isRunning, isNmieLoaded,
      runNearFieldSimulation,
      saveSpectrumSimulation(){
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
        let xTitle = 'x'
        if ( $store.state.plotRuntime.spectrumPlots.layout.xaxis ) {
          xTitle = String($store.state.plotRuntime.spectrumPlots.layout.xaxis.title)
        }

        let columnNames = '# ' + xTitle + ', Qsca, Qabs, Qext, '
        const mode_n = rangeInt($store.state.simulationSetup.current.numberOfModesToPlot, 1);
        const mode_types = range(0, 1);
        for (const n of mode_n) {
          for (const mode_type of mode_types) {
            const modeTypeName = mode_type == 0 ? 'E' : 'H'
            columnNames += 'Qsca_' + modeTypeName + '_' +getModeName(n)+', '
            columnNames += 'Qabs_' + modeTypeName + '_' +getModeName(n)+', '
            columnNames += 'Qext_' + modeTypeName + '_' +getModeName(n)+', '
          }
        }
        columnNames = columnNames.slice(0, -2)
        columnNames += '\n'
        let body = ''
        const WLs = $store.state.plotRuntime.WLsInUnits
        const Qsca = $store.state.plotRuntime.Qsca
        const Qabs = $store.state.plotRuntime.Qabs
        const Qext = $store.state.plotRuntime.Qext
        const Qsca_n = $store.state.plotRuntime.Qsca_n
        const Qabs_n = $store.state.plotRuntime.Qabs_n
        const Qext_n = $store.state.plotRuntime.Qext_n
        for (let i = 0; i < WLs.length; ++i) {
          let row = WLs[i].toString() + ', '
              + Qsca[i].toString() + ', '
              + Qabs[i].toString() + ', '
              + Qext[i].toString() + ', '
          for (const n of mode_n) {
            for (const mode_type of mode_types) {
              row += Qsca_n[mode_type][n - 1][i].toString() + ', '
              row += Qabs_n[mode_type][n - 1][i].toString() + ', '
              row += Qext_n[mode_type][n - 1][i].toString() + ', '
            }
          }
          row = row.slice(0, -2)
          row += '\n'
          body += row
        }

        const scattnlaySpectra = new Blob([fileHeader+columnNames+body],
            {type: 'text/plain;charset=utf-8',
              endings: 'native'}  //TODO test if newline is correctly written in Windows, MacOS
        )
        saveAs(scattnlaySpectra, 'scattnlay-spectra.txt');
      }
    }
  },
})
</script>
