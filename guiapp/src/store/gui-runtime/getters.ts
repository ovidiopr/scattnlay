import { GetterTree } from 'vuex';
import { StateInterface } from '../index';
import { guiRuntimeStateInterface } from './state';

const getters: GetterTree<guiRuntimeStateInterface, StateInterface> = {
  someAction (/* context */) {
    // your code
  }
};

export default getters;
