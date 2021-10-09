import { ActionTree } from 'vuex';
import { StateInterface } from '../index';
import { plotRuntimeStateInterface } from './state';

const actions: ActionTree<plotRuntimeStateInterface, StateInterface> = {
  someAction (/* context */) {
    // your code
  }
};

export default actions;
