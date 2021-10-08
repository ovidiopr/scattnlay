import { Module } from 'vuex';
import { StateInterface } from '../index';
import state, { guiRuntimeStateInterface } from './state';
import actions from './actions';
import getters from './getters';
import mutations from './mutations';

const guiRuntimeModule: Module<guiRuntimeStateInterface, StateInterface> = {
  namespaced: true,
  actions,
  getters,
  mutations,
  state
};

export default guiRuntimeModule;
