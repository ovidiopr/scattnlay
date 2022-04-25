import { Module } from 'vuex';
import { StateInterface } from '../index';
import state, { simulationSetupStateInterface } from './state';
import actions from './actions';
import getters from './getters';
import mutations from './mutations';

const simulationSetupModule: Module<simulationSetupStateInterface, StateInterface> = {
  namespaced: true,
  actions,
  getters,
  mutations,
  state
};

export default simulationSetupModule;
