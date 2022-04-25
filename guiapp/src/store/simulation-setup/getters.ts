import { GetterTree } from 'vuex';
import { StateInterface } from '../index';
import { simulationSetupStateInterface } from './state';

const getters: GetterTree<simulationSetupStateInterface, StateInterface> = {
  someAction (/* context */) {
    // your code
  }
};

export default getters;
