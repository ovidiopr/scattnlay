import { MutationTree } from 'vuex';
import { simulationSetupStateInterface as sssi, simulationSetup, layer } from './state';
import { cloneDeep } from 'lodash'

const mutation: MutationTree<sssi> = {
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
  setNumberOfModesToPlot  (state: sssi, val: number) {state.gui.numberOfModesToPlot  = val},

};

export default mutation;
