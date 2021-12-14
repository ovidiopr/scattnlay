<template>
  <div class="row items-baseline">
    <div class="col-xs-12 col-sm-auto text-weight-bold text-center q-px-md q-py-sm">
        <q-btn :loading="isRunning"
               :disable="isRunning||!isNmieLoaded"
               color="primary"
               no-caps
               :label="isNmieLoaded ? 'Run simulation' : 'Loading...'"
               @click="runSpectrumSimulation">
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
  computed, watch,
} from 'vue'
import { useStore } from 'src/store'
import { getModeName, range, rangeInt} from 'components/utils'
import { cloneDeep } from 'lodash'
import { saveAs } from 'file-saver'

export default defineComponent({
  name: 'RunSimulationSpectrum',

  setup() {
    const $store = useStore()

    const isRunning = computed({
      get: ()=> $store.state.simulationSetup.isNmieRunning,
      set: val => {
        val ? $store.commit('simulationSetup/markNmieAsStarted') : $store.commit('simulationSetup/markNmieAsFinished')
      }
    })

    const isNmieLoaded = computed(()=>{ return $store.state.simulationSetup.isNmieLoaded })

    const sourceUnits = computed( ()=>$store.state.guiRuntime.sourceUnits)

    function getWLs(){
      const fromWL = $store.state.simulationSetup.current.fromWL
      const toWL = $store.state.simulationSetup.current.toWL
      const pointsWL = $store.state.simulationSetup.current.pointsWL
      if (sourceUnits.value.endsWith('Hz') || sourceUnits.value.endsWith('eV')) {
        const fromF = 1./fromWL
        const toF = 1./toWL
        const stepF = (fromF-toF)/(pointsWL-1)
        const Fs = range(toF, fromF, stepF);
        let WLs = []
        for (const f of Fs) WLs.push(1./f)
        return WLs
      }
      const stepWL = (toWL-fromWL)/(pointsWL-1)
      const WLs = range(fromWL, toWL, stepWL);
      return WLs
    }

    function initQ(mode_n:number[], mode_types:number[]) {
      let Qsca:number[] = [], Qabs:number[] = [], Qext:number[] = []
      let Qsca_n:number[][][] = [[], []]
      let Qabs_n:number[][][] = [[], []]
      let Qext_n:number[][][] = [[], []]
      mode_types.forEach(function (mode_type) {
        mode_n.forEach(function () {
          Qsca_n[mode_type].push([])
          Qabs_n[mode_type].push([])
          Qext_n[mode_type].push([])
        })
      })
      return {Qsca, Qabs, Qext, Qsca_n, Qabs_n, Qext_n}
    }

    //-------------------------------------------------------------------------//
    //---------------------  Main  --------------------------------------------//
    //-------------------------------------------------------------------------//
    function runSpectrumSimulation() {
      if (isRunning.value) {
        console.log('Some Nmie is already running!')
        return
      }
      isRunning.value = true
      setTimeout(()=> {

        $store.commit('simulationSetup/copySetupFromGuiToCurrent')

        const host = $store.state.simulationSetup.current.hostIndex
        const WLs = getWLs()
        const mode_n = rangeInt($store.state.simulationSetup.current.numberOfModesToPlot, 1);
        const mode_types = range(0, 1);
        let {Qsca, Qabs, Qext, Qsca_n, Qabs_n, Qext_n} = initQ(mode_n, mode_types)

        try {
          if (!$store.state.simulationSetup.nmie) throw 'ERROR! Scattnlay module was not loaded'
          const nmie = $store.state.simulationSetup.nmie
          const layers = cloneDeep($store.state.simulationSetup.current.layers)
          const nmieStartedTime = performance.now()
          for (const WL of WLs) {
            nmie.SetWavelength(WL)

            nmie.ClearTarget()
            for (const layer of layers) {
              if (layer.material.nSpline) layer.n = layer.material.nSpline.at(WL)
              if (layer.material.kSpline) layer.k = layer.material.kSpline.at(WL)
              nmie.AddTargetLayerReIm(layer.layerWidth * host, layer.n / host, layer.k / host)
            }

            nmie.SetModeNmaxAndType(-1, -1)
            nmie.RunMieCalculation()
            Qsca.push(nmie.GetQsca())
            Qabs.push(nmie.GetQabs())
            Qext.push(nmie.GetQext())

            for (const mode_type of mode_types) {
              for (const n of mode_n) {
                nmie.SetModeNmaxAndType(n, mode_type)
                nmie.RunMieCalculation()
                Qsca_n[mode_type][n - 1].push(nmie.GetQsca())
                Qabs_n[mode_type][n - 1].push(nmie.GetQabs())
                Qext_n[mode_type][n - 1].push(nmie.GetQext())
              }
            }
          }
          const nmieTotalRunTime = (performance.now()-nmieStartedTime)/1000
          // console.log('Total simulation time:', nmieTotalRunTime, 's')
          $store.commit('simulationSetup/setNmieTotalRunTime', nmieTotalRunTime)

          $store.commit('plotRuntime/setQ', {WLs, Qsca, Qabs, Qext, Qsca_n, Qabs_n, Qext_n})
          $store.commit('plotRuntime/setWLsInUnits', sourceUnits.value)

          $store.commit('plotRuntime/updateNumberOfPlotsFromPreviousSimulations')
          $store.commit('plotRuntime/setCommonLabel', $store.state.simulationSetup.current.plotLabel)
          $store.commit('simulationSetup/setPlotLabel', '')
          $store.commit('plotRuntime/updateSpectraPlot')
        } catch (e) {
          console.log(e)
        }
        isRunning.value = false
      },200)
    }

    watch(isNmieLoaded, ()=>{
      if (isNmieLoaded.value) runSpectrumSimulation()
    })

    return { isRunning, isNmieLoaded,
      runSpectrumSimulation,
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
        if ( $store.state.plotRuntime.spectraPlot.layout.xaxis ) {
          xTitle = String($store.state.plotRuntime.spectraPlot.layout.xaxis.title)
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
