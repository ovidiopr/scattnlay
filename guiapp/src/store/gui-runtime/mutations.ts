import { MutationTree } from 'vuex';
import { guiRuntimeStateInterface as grsi} from './state';

const mutation: MutationTree<grsi> = {
  setIsShowingHelpForInputWithUnits (state: grsi, val: boolean) {state.isShowingHelpForInputWithUnits = val},
  setUnits             (state: grsi, val: string ) {state.units             = val},
  setSourceUnits       (state: grsi, val: string ) {state.sourceUnits       = val},
  setIsSourceSameUnits (state: grsi, val: boolean) {state.isSourceSameUnits = val},
};

export default mutation;
