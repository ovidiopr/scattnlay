import { Module } from 'vuex';
import { StateInterface } from '../index';
import state, { plotRuntimeStateInterface } from './state';
import actions from './actions';
import getters from './getters';
import mutations from './mutations';

const plotRuntimeModule: Module<plotRuntimeStateInterface, StateInterface> = {
  namespaced: true,
  actions,
  getters,
  mutations,
  state
};

export default plotRuntimeModule;
