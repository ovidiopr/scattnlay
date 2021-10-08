import { ActionTree } from 'vuex';
import { StateInterface } from '../index';
import { guiRuntimeStateInterface } from './state';

const actions: ActionTree<guiRuntimeStateInterface, StateInterface> = {
  someAction (/* context */) {
    // your code
  }
};

export default actions;
