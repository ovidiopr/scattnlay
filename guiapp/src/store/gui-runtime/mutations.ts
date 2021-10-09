import { MutationTree } from 'vuex';
import { guiRuntimeStateInterface } from './state';

const mutation: MutationTree<guiRuntimeStateInterface> = {
  setIsShowingHelpForInputWithUnits (state: guiRuntimeStateInterface, newVal: boolean) {
    // your code
    state.isShowingHelpForInputWithUnits = newVal
  }
};

export default mutation;
