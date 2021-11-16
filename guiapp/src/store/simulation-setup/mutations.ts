import { MutationTree } from 'vuex';
import { simulationSetupStateInterface as sssi, simulationSetup, layer } from './state';
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
  setCurrentWL (state: sssi, val: number) {state.gui.currentWL = val},
  setNumberOfModesToPlot  (state: sssi, val: number) {state.gui.numberOfModesToPlot  = val},
  setPlotLabel (state: sssi, val: string) {state.gui.plotLabel = val},

};

export default mutation;
