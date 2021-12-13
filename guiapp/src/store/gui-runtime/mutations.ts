import { MutationTree } from 'vuex';
import { guiRuntimeStateInterface as grsi} from './state';
import { composeLabelFromPageData } from 'components/utils'

const mutation: MutationTree<grsi> = {
  setIsShowingHelpForInputWithUnits (state: grsi, val: boolean) {state.isShowingHelpForInputWithUnits = val},
  setUnits             (state: grsi, val: string ) {state.units             = val},
  setSourceUnits       (state: grsi, val: string ) {state.sourceUnits       = val},
  setIsSourceSameUnits (state: grsi, val: boolean) {state.isSourceSameUnits = val},

  activateMaterial(state: grsi, filepath: string) {
    console.log(composeLabelFromPageData(filepath))
    return {state, filepath}
  },

  deleteMaterial(state: grsi, label: string) {
    const indexToDelete = state.activatedMaterials.findIndex(val => val.name==label)
    state.activatedMaterials.splice(indexToDelete,1)
  },

};

export default mutation;
