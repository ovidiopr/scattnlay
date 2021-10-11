import { MutationTree } from 'vuex';
import { simulationSetupStateInterface, simulationSetup } from './state';
import { cloneDeep } from 'lodash'

const mutation: MutationTree<simulationSetupStateInterface> = {
  setGuiState (state: simulationSetupStateInterface,
               newVal: simulationSetup) {
    // your code
    state.gui = cloneDeep(newVal)
  }
};

export default mutation;
