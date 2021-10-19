import { MutationTree } from 'vuex';
import { simulationSetupStateInterface as sssi, simulationSetup, layer } from './state';
import { cloneDeep } from 'lodash'

const mutation: MutationTree<sssi> = {
  setGuiState (state: sssi,
               newVal: simulationSetup) {
    state.gui = cloneDeep(newVal)
    // // Possible usage in component
    // let simulationSetupGui = reactive(cloneDeep($store.state.simulationSetup.gui))
    // const unsubscribe = $store.subscribe((mutation, /*state*/) => {
    //   if (mutation.type === 'simulationSetup/setGuiState') {
    //     let key: keyof typeof simulationSetupGui
    //     for (key in simulationSetupGui) {
    //       simulationSetupGui[key] = $store.state.simulationSetup.gui[key]
    //     }
    //   }
    // })
    // onBeforeUnmount(()=>unsubscribe())
    // watch(simulationSetupGui, () => $store.commit('simulationSetup/setGuiState',simulationSetupGui))
  },
  setHostIndex (state: sssi, val: number) {state.gui.hostIndex = val},
  setFromWL    (state: sssi, val: number) {state.gui.fromWL    = val},
  setToWL      (state: sssi, val: number) {state.gui.toWL      = val},
  setPointsWL  (state: sssi, val: number) {state.gui.pointsWL  = val},

  setLayers    (state: sssi, val: layer[]) {state.gui.layers = cloneDeep(val) }
};

export default mutation;
