import { ActionTree } from 'vuex';
import { StateInterface } from '../index';
import { simulationSetupStateInterface } from './state';

const actions: ActionTree<simulationSetupStateInterface, StateInterface> = {
  someAction (/* context */) {
    // your code
  }
};

export default actions;
