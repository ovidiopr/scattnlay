import { GetterTree } from 'vuex';
import { StateInterface } from '../index';
import { plotRuntimeStateInterface } from './state';

const getters: GetterTree<plotRuntimeStateInterface, StateInterface> = {
  someAction (/* context */) {
    // your code
  }
};

export default getters;
