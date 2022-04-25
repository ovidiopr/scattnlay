import { ActionTree } from 'vuex';
import { StateInterface } from '../index';
import { simulationSetupStateInterface } from './state';
import nmiejs from 'src/nmiejs.js';
import { useQuasar } from 'quasar';

const actions: ActionTree<simulationSetupStateInterface, StateInterface> = {
  async loadScattnlay({ commit /*state*/ } /* context */) {
    const $q = useQuasar();
    $q.loading.show({
      message: 'Loading Mie calculator. Please wait...',
      boxClass: 'bg-grey-2 text-grey-9',
      spinnerColor: 'primary',
    });
    const module = await nmiejs();
    const nmies = {
      spectrum: new module.nmie(),
      nearField: new module.nmie(),
      farField: new module.nmie(),
    };
    commit('setNmies', nmies);

    // // Test nmiejs if working
    // if (state.nmies.spectrum.instance && !state.nmies.spectrum.isNmieRunning) {
    //   commit('markNmieAsStarted')
    //   state.nmies.spectrum.instance.ClearTarget()
    //   const R = 100.0
    //   const reN = 4.0
    //   const imN = 0.01
    //   state.nmies.spectrum.instance.AddTargetLayerReIm(R, reN, imN)
    //   state.nmies.spectrum.instance.SetModeNmaxAndType(-1, -1)
    //   const WL = 800
    //   state.nmies.spectrum.instance.SetWavelength(WL)
    //   state.nmies.spectrum.instance.RunMieCalculation()
    //   console.log(state.nmies.spectrum.instance.GetQsca())
    //   // outer_arc_points, radius_points, from_Rho, to_Rho,
    //   // from_Theta, to_Theta, from_Phi, to_Phi, isIgnoreAvailableNmax
    //   state.nmies.spectrum.instance.RunFieldCalculationPolar(2, 2,
    //       0.1, 1.5, 0, 3.1415, 0, 3.1415,
    //       0)
    //   console.log('Field Eabs:', state.nmies.spectrum.instance.GetFieldEabs())
    //   commit('markNmieAsFinished')
    // }

    $q.loading.hide();
  },
};

export default actions;
