import { MutationTree } from 'vuex';
import { simulationSetupStateInterface } from './state';

const mutation: MutationTree<simulationSetupStateInterface> = {
  setIsShowingHelpForInputWithUnits (state: simulationSetupStateInterface, newVal: boolean) {
    // your code
    state.isShowingHelpForInputWithUnits = newVal
  }
};

export default mutation;
