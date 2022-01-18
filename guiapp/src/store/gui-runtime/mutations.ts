import { MutationTree } from 'vuex';
import { guiRuntimeStateInterface as grsi } from './state';
import { material } from 'src/store/simulation-setup/state';
import { cloneDeep } from 'lodash';

function compare(a: material, b: material) {
  if (a.name == 'link') return -1;
  if (b.name == 'link') return 1;
  if (a.name == 'nk-constant') return -1;
  if (b.name == 'nk-constant') return 1;

  if (a.name < b.name) return -1;
  if (a.name > b.name) return 1;
  return 0;
}

const mutation: MutationTree<grsi> = {
  setIsSaveWithPythonScript(state: grsi, val: boolean) {
    state.isSaveWithPythonScript = val;
  },
  setIsAutoRefineNearField(state: grsi, val: boolean) {
    state.isAutoRefineNearField = val;
  },
  setIsShowingHelpForInputWithUnits(state: grsi, val: boolean) {
    state.isShowingHelpForInputWithUnits = val;
  },
  setUnits(state: grsi, val: string) {
    state.units = val;
  },
  setSourceUnits(state: grsi, val: string) {
    state.sourceUnits = val;
  },
  setIsSourceSameUnits(state: grsi, val: boolean) {
    state.isSourceSameUnits = val;
  },

  setNearFieldZoom(
    state: grsi,
    val: { fromX: number; toX: number; fromY: number; toY: number }
  ) {
    state.nearFieldZoom = cloneDeep(val);
  },

  setSafeWL(state: grsi, val: { safeFromWL: number; safeToWL: number }) {
    state.safeFromWL = val.safeFromWL;
    state.safeToWL = val.safeToWL;
  },

  addMaterial(state: grsi, material: material) {
    state.activatedMaterials.push(material);
    state.activatedMaterials.sort(compare);
  },

  deleteMaterial(state: grsi, label: string) {
    const indexToDelete = state.activatedMaterials.findIndex(
      (val) => val.name == label
    );
    state.activatedMaterials.splice(indexToDelete, 1);
  },
  // @click="$store.commit('guiRuntime/toggleIsPlot',
  // props.row.name

  toggleIsPlot(state: grsi, label: string) {
    const indexToToggle = state.activatedMaterials.findIndex(
      (val) => val.name == label
    );
    state.activatedMaterials[indexToToggle].isPlot =
      !state.activatedMaterials[indexToToggle].isPlot;
  },
};

export default mutation;
