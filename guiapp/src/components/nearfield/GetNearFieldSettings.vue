<template>
  <div class="row items-baseline">
    <div
      class="col-xs-12 col-sm-auto text-weight-bold text-center q-pr-md q-py-sm"
    >
      <div :style="flexRowTitleStyle">Plot</div>
    </div>
    <div class="col-xs-grow col-sm">
      <div class="row justify-center items-baseline">
        <div class="col-auto text-center q-py-xs q-pr-md">
          <div :style="flexRowTitleStyle">x-side resolution</div>
        </div>
        <div class="col-xs-grow col-sm">
          <div class="row justify-xs-center justify-sm-start items-center">
            <div class="col-auto">
              <input-with-units
                v-model:input-result="plotXSideResolution"
                v-model:is-showing-help="isShowingHelpForInputWithUnits"
                :initial-expression="plotXSideResolution.toString()"
                title="points"
                units=""
              />
            </div>
          </div>
        </div>
      </div>
      <div class="q-ma-xs" />
      <div class="row justify-center items-center">
        <div class="col-auto text-center q-py-xs q-pr-md">
          <div :style="flexRowTitleStyle">ratio</div>
        </div>
        <div class="col-xs-grow col-sm">
          <div class="row justify-xs-center justify-sm-start items-center">
            <div class="q-gutter-md q-py-sm q-px-xs">
              <q-radio
                v-model="plotRatioLabel"
                dense
                size="sm"
                :val="'any'"
                label="any"
              />
              <q-radio
                v-model="plotRatioLabel"
                dense
                size="sm"
                :val="'fixed'"
                label="fixed"
              />
              <q-radio
                v-model="plotRatioLabel"
                dense
                size="sm"
                :val="'1:1'"
                label="1:1"
              />
              <q-radio
                v-model="plotRatioLabel"
                dense
                size="sm"
                :val="'3:2'"
                label="3:2"
              />
              <q-radio
                v-model="plotRatioLabel"
                dense
                size="sm"
                :val="'2:1'"
                label="2:1"
              />
            </div>
          </div>
        </div>
      </div>

      <div
        v-if="plotRatioLabel == 'any'"
        class="row justify-center items-baseline"
      >
        <div class="col-auto text-center q-py-xs q-pr-md">
          <div :style="flexRowTitleStyle">y-side resolution</div>
        </div>
        <div class="col-xs-grow col-sm">
          <div class="row justify-xs-center justify-sm-start items-center">
            <input-with-units
              v-model:input-result="plotYSideResolution"
              v-model:is-showing-help="isShowingHelpForInputWithUnits"
              :initial-expression="plotYSideResolution.toString()"
              :is-info-mode="plotRatioLabel != 'any'"
              title="points"
              units=""
            />
          </div>
        </div>
      </div>

      <div class="row justify-center items-baseline">
        <div class="col-auto text-center q-py-xs q-pr-md">
          <div :style="flexRowTitleStyle">relative side length</div>
        </div>
        <div class="col-xs-grow col-sm">
          <div class="row justify-xs-center justify-sm-start items-center">
            <input-with-units
              v-model:input-result="relativePlotSize"
              v-model:is-showing-help="isShowingHelpForInputWithUnits"
              :initial-expression="relativePlotSize.toString()"
              title="𝐿&thinsp;/&hairsp;𝟐𝑅"
              units=""
            />
          </div>
        </div>
      </div>
      <div class="q-ma-xs" />
      <div class="row justify-center items-center">
        <div class="col-auto text-center q-py-xs q-pr-md">
          <div :style="flexRowTitleStyle">cross-section</div>
        </div>
        <div class="col-xs-grow col-sm">
          <div class="row justify-xs-center justify-sm-start items-center">
            <div class="q-gutter-md q-py-sm q-px-xs">
              <q-radio
                v-model="crossSection"
                dense
                size="sm"
                :val="nearFieldPlane.Ek"
                label="Ek"
              />
              <q-radio
                v-model="crossSection"
                dense
                size="sm"
                :val="nearFieldPlane.Hk"
                label="Hk"
              />
              <q-radio
                v-model="crossSection"
                dense
                size="sm"
                :val="nearFieldPlane.EH"
                label="EH"
              />
            </div>
          </div>
        </div>
      </div>
    </div>
  </div>
</template>

<script lang="ts">
import { defineComponent, computed, watch } from 'vue';
import { useStore } from 'src/store';
import InputWithUnits from 'components/InputWithUnits.vue';
import { flexRowTitleStyle } from 'components/config';
import { nearFieldPlane } from 'src/store/simulation-setup/state';

export default defineComponent({
  name: 'GetNearFieldSettings',
  components: { InputWithUnits },

  setup() {
    const $store = useStore();

    const isShowingHelpForInputWithUnits = computed({
      get: () => $store.state.guiRuntime.isShowingHelpForInputWithUnits,
      set: (val) =>
        $store.commit('guiRuntime/setIsShowingHelpForInputWithUnits', val),
    });

    const crossSection = computed({
      get: () => $store.state.simulationSetup.gui.nearFieldSetup.crossSection,
      set: (val) =>
        $store.commit('simulationSetup/setNearFieldCrossSection', val),
    });

    const relativePlotSize = computed({
      get: () =>
        $store.state.simulationSetup.gui.nearFieldSetup.relativePlotSize,
      set: (val) =>
        $store.commit('simulationSetup/setNearFieldRelativePlotSize', val),
    });

    // const maxComputeTime = computed({
    //   get: () => $store.state.simulationSetup.gui.nearFieldSetup.maxComputeTime,
    //   set: val => $store.commit('simulationSetup/setNearFieldMaxComputeTime', val)
    // })

    const plotXSideResolution = computed({
      get: () =>
        $store.state.simulationSetup.gui.nearFieldSetup.plotXSideResolution,
      // TODO: make InputWithUnits to handle integer input, so no need to use floor() in the next line.
      set: (val) =>
        $store.commit(
          'simulationSetup/setNearFieldPlotXSideResolution',
          Math.floor(val)
        ),
    });

    const plotYSideResolution = computed({
      get: () =>
        $store.state.simulationSetup.gui.nearFieldSetup.plotYSideResolution,
      set: (val) =>
        $store.commit(
          'simulationSetup/setNearFieldPlotYSideResolution',
          Math.floor(val)
        ),
    });

    const plotRatioLabel = computed({
      get: () => $store.state.guiRuntime.plotRatioLabel,
      set: (val) => $store.commit('guiRuntime/setPlotRatioLabel', val),
    });
    const plotRatio = computed({
      get: () => $store.state.guiRuntime.plotRatio,
      set: (val) => $store.commit('guiRuntime/setPlotRatio', val),
    });

    watch(plotRatioLabel, () => {
      switch (plotRatioLabel.value) {
        case '1:1':
          plotRatio.value = 1;
          break;
        case '3:2':
          plotRatio.value = 2 / 3;
          break;
        case '2:1':
          plotRatio.value = 1 / 2;
          break;
        case 'fixed':
          plotRatio.value =
            plotYSideResolution.value / plotXSideResolution.value;
          break;
        default:
          break;
      }
    });
    watch([plotXSideResolution, plotRatioLabel, plotRatio], () => {
      if (plotRatioLabel.value != 'any')
        plotYSideResolution.value = plotRatio.value * plotXSideResolution.value;
    });
    return {
      crossSection,
      plotRatioLabel,
      isShowingHelpForInputWithUnits,
      flexRowTitleStyle,
      relativePlotSize,
      // maxComputeTime,
      plotXSideResolution,
      plotYSideResolution,
      nearFieldPlane,
    };
  },
});
</script>
