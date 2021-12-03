import { ActionTree } from 'vuex';
import { StateInterface } from '../index';
import { simulationSetupStateInterface } from './state';
import nmiejs from 'src/nmiejs.js';


const actions: ActionTree<simulationSetupStateInterface, StateInterface> = {
  async loadScattnlay ({commit,state}/* context */) {
    const module = await nmiejs()
    const nmie = new module.nmie()
    commit('setNmie', nmie)
    // // Test nmiejs if working
    // if (state.nmie && !state.isNmieRunning) {
    //   commit('markNmieAsStarted')
    //   state.nmie.ClearTarget()
    //   const R = 100.0
    //   const reN = 4.0
    //   const imN = 0.01
    //   state.nmie.AddTargetLayerReIm(R, reN, imN)
    //   state.nmie.SetModeNmaxAndType(-1, -1)
    //   const WL = 800
    //   state.nmie.SetWavelength(WL)
    //   state.nmie.RunMieCalculation()
    //   console.log(state.nmie.GetQsca())
    //   // outer_arc_points, radius_points, from_Rho, to_Rho,
    //   // from_Theta, to_Theta, from_Phi, to_Phi, isIgnoreAvailableNmax
    //   state.nmie.RunFieldCalculationPolar(2, 2,
    //       0.1, 1.5, 0, 3.1415, 0, 3.1415,
    //       0)
    //   console.log('Field Eabs:', state.nmie.GetFieldEabs())
    //   commit('markNmieAsFinished')
    // }
    if (state.nmie) {
      commit('markNmieAsLoaded')
    }
  }
};

export default actions;
