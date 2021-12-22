import { MutationTree } from 'vuex';
import { simulationSetupStateInterface as sssi, simulationSetup, layer, nearFieldType } from './state';
import { cloneDeep } from 'lodash'
import { markRaw} from 'vue'

const mutation: MutationTree<sssi> = {
  setNmie (state: sssi,
           newVal: import('src/nmiejs').nmie_class) {
    state.nmie = markRaw(newVal)
  },
  markNmieAsLoaded  (state: sssi) {state.isNmieLoaded  = true },
  markNmieAsStarted (state: sssi) {
    state.isNmieRunning = true
  },
  markNmieAsFinished(state: sssi) {
    state.isNmieRunning = false
  },
  setNmieTotalRunTime(state: sssi, val:number) {state.nmieTotalRunTime = val},

  setCurrentState (state: sssi,
                   newVal: simulationSetup) {
    state.current = cloneDeep(newVal)
  },

  copySetupFromGuiToCurrent (state: sssi) {
    state.current = cloneDeep(state.gui)
  },

  // Mutations for simulation setup as represented\set in GUI
  setGuiState (state: sssi,
               newVal: simulationSetup) {
    state.gui = cloneDeep(newVal)
  },
  setLayers    (state: sssi,
                newVal: layer[]) {
    state.gui.layers = cloneDeep(newVal)
  },
  setHostIndex (state: sssi, val: number) {state.gui.hostIndex = val},
  setFromWL    (state: sssi, val: number) {state.gui.fromWL    = val},
  setToWL      (state: sssi, val: number) {state.gui.toWL      = val},
  setPointsWL  (state: sssi, val: number) {state.gui.pointsWL  = val},
  setPlotLabel (state: sssi, val: string) {state.gui.plotLabel = val},
  setNumberOfModesToPlot  (state: sssi, val: number) {state.gui.numberOfModesToPlot  = val},

  setNearFieldWL                 (state: sssi, val: number)        {state.gui.nearFieldSetup.atWL               = val},
  setNearFieldRelativePlotSize   (state: sssi, val: number)        {state.gui.nearFieldSetup.relativePlotSize   = val},
  setNearFieldPlotSideResolution (state: sssi, val: number)        {state.gui.nearFieldSetup.plotSideResolution = val},
  setNearFieldCrossSection       (state: sssi, val: nearFieldType) {state.gui.nearFieldSetup.crossSection       = val},
  setNearFieldMaxComputeTime     (state: sssi, val: number)        {state.gui.nearFieldSetup.maxComputeTime     = val},

  setFarFieldWL    (state: sssi, val: number) { state.gui.farFieldWL  = val},

};

export default mutation;
