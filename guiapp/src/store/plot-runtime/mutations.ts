import { MutationTree } from 'vuex';
import { plotRuntimeStateInterface } from './state';

const mutation: MutationTree<plotRuntimeStateInterface> = {
  setIsShowingHelpForInputWithUnits (state: plotRuntimeStateInterface, newVal: boolean) {
    // your code
    state.isShowingHelpForInputWithUnits = newVal
  }
};

export default mutation;
