import { store } from 'quasar/wrappers';
import { InjectionKey } from 'vue';
import { Router } from 'vue-router'
import {
  createStore,
  Store as VuexStore,
  useStore as vuexUseStore,
} from 'vuex';

import guiRuntime from './gui-runtime';
import { guiRuntimeStateInterface } from './gui-runtime/state';

import plotRuntime from './plot-runtime';
import { plotRuntimeStateInterface } from './plot-runtime/state';

import simulationSetup from './simulation-setup';
import { simulationSetupStateInterface } from './simulation-setup/state';

/*
 * If not building with SSR mode, you can
 * directly export the Store instantiation;
 *
 * The function below can be async too; either use
 * async/await or return a Promise which resolves
 * with the Store instance.
 */

export interface StateInterface {
  // Define your own store structure, using submodules if needed
  guiRuntime: guiRuntimeStateInterface;
  plotRuntime: plotRuntimeStateInterface;
  simulationSetup: simulationSetupStateInterface;
  // Declared as unknown to avoid linting issue. Best to strongly type as per the line above.
  // example: unknown
}

// provide typings for `this.$store`
declare module '@vue/runtime-core' {
  interface ComponentCustomProperties {
    $store: VuexStore<StateInterface>;
  }
}

// provide typings for `useStore` helper
export const storeKey: InjectionKey<VuexStore<StateInterface>> =
 Symbol('vuex-key');

// Provide typings for `this.$router` inside Vuex stores
declare module 'vuex' {
  // eslint-disable-next-line @typescript-eslint/no-unused-vars
  export interface Store<S> {
    readonly $router: Router;
  }
}

export default store(function (/* { ssrContext } */) {
  const Store = createStore<StateInterface>({
    modules: {
      guiRuntime,
      plotRuntime,
      simulationSetup,
    },

    // enable strict mode (adds overhead!)
    // for dev mode and --debug builds only
    strict: !!process.env.DEBUGGING,
  });

  return Store;
});

export function useStore() {
  return vuexUseStore(storeKey);
}
