<template>
  <div class="row items-baseline">
    <div
      class="col-xs-12 col-sm-auto text-weight-bold text-center q-px-md q-py-sm"
    >
      <q-tooltip
        v-if="
          $store.state.guiRuntime.safeFromWL >
            $store.state.simulationSetup.gui.nearFieldSetup.atWL ||
          $store.state.guiRuntime.safeToWL <
            $store.state.simulationSetup.gui.nearFieldSetup.atWL
        "
        anchor="top middle"
        self="center middle"
        class="bg-amber-4 text-black shadow-4"
      >
        Will use materials<br />
        spectrum range.
      </q-tooltip>
      <q-btn
        :loading="isRunning"
        :disable="isRunning || !isNmieLoaded"
        color="primary"
        no-caps
        :label="isNmieLoaded ? 'Run simulation' : 'Loading...'"
        @click="runNearFieldSimulation"
      >
        <template #loading>
          <q-spinner-gears />
        </template>
      </q-btn>
    </div>
    <div class="col-xs-grow col-sm q-px-xs">
      <div class="row justify-xs-center justify-sm-start items-baseline">
        <div class="col-auto">
          <SaveSimulationNearField />
        </div>
      </div>
    </div>
  </div>
</template>

<script lang="ts">
import { computed, defineComponent, nextTick } from 'vue';
import { useStore } from 'src/store';
import { cloneDeep } from 'lodash';
import SaveSimulationNearField from 'components/nearfield/SaveSimulationNearField.vue';

export default defineComponent({
  name: 'RunSimulationNearField',
  components: { SaveSimulationNearField },
  setup() {
    const $store = useStore();

    const isRunning = computed({
      get: () => $store.state.simulationSetup.nmies.nearField.isNmieRunning,
      set: (val) => {
        val
          ? $store.commit('simulationSetup/markNmieNearFieldAsStarted')
          : $store.commit('simulationSetup/markNmieNearFieldAsFinished');
      },
    });

    const isNmieLoaded = computed(() => {
      return $store.state.simulationSetup.nmies.nearField.instance;
    });

    const atWL = computed(
      () => $store.state.simulationSetup.current.nearFieldSetup.atWL
    );
    //-------------------------------------------------------------------------//
    //---------------------  Main  --------------------------------------------//
    //-------------------------------------------------------------------------//
    function runNearFieldSimulation() {
      if (isRunning.value) {
        console.log('Some Nmie is already running!');
        return;
      }
      isRunning.value = true;
      void setTimeout(() => {
        void nextTick(() => {
          $store.commit('simulationSetup/copySetupFromGuiToCurrent');

          const host = $store.state.simulationSetup.current.hostIndex;
          const plotXSideResolution =
            $store.state.simulationSetup.current.nearFieldSetup
              .plotXSideResolution;
          const plotYSideResolution =
            $store.state.simulationSetup.current.nearFieldSetup
              .plotYSideResolution;
          const relativePlotSize =
            $store.state.simulationSetup.current.nearFieldSetup
              .relativePlotSize;
          const crossSection =
            $store.state.simulationSetup.current.nearFieldSetup.crossSection;
          const atX =
            $store.state.simulationSetup.current.nearFieldSetup.atRelativeX0;
          const atY =
            $store.state.simulationSetup.current.nearFieldSetup.atRelativeY0;
          const atZ =
            $store.state.simulationSetup.current.nearFieldSetup.atRelativeZ0;
          const isMarkUnconvergedAsNaN = 1; // 0 - do not mark, else do mark
          try {
            if (!$store.state.simulationSetup.nmies.nearField.instance)
              throw 'ERROR! Scattnlay module was not loaded';
            const nmie = $store.state.simulationSetup.nmies.nearField.instance;
            const layers = cloneDeep(
              $store.state.simulationSetup.current.layers
            );
            const nmieStartedTime = performance.now();

            nmie.SetWavelength(atWL.value);
            nmie.ClearTarget();
            for (const layer of layers) {
              if (layer.material.nSpline)
                layer.n = layer.material.nSpline.at(atWL.value);
              if (layer.material.kSpline)
                layer.k = layer.material.kSpline.at(atWL.value);
              nmie.AddTargetLayerReIm(
                layer.layerWidth * host,
                layer.n / host,
                layer.k / host
              );
            }

            nmie.SetModeNmaxAndType(-1, -1);

            nmie.RunFieldCalculationCartesian(
              plotYSideResolution,
              plotXSideResolution, // in simulation z-axis resolution for Ek and Hk cross-sections
              relativePlotSize,
              crossSection,
              atX,
              atY,
              atZ,
              isMarkUnconvergedAsNaN
            );
            const Eabs = nmie.GetFieldEabs();
            $store.commit('plotRuntime/setNearField', Eabs);
            const nmieTotalRunTime =
              (performance.now() - nmieStartedTime) / 1000;
            // console.log('Total simulation time:', nmieTotalRunTime, 's')
            $store.commit(
              'simulationSetup/setNmieNearFieldTotalRunTime',
              nmieTotalRunTime
            );
          } catch (e) {
            console.log('Some error:', e);
          }
          isRunning.value = false;
        });
      }, 100);
    }

    runNearFieldSimulation();

    return { isRunning, isNmieLoaded, runNearFieldSimulation };
  },
});
</script>
