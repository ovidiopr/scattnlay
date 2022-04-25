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
          <div :style="flexRowTitleStyle" @click="hover = !hover">ratio</div>
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

      <!--      <transition-->
      <!--        name="scale-y"-->
      <!--        @before-enter="beforeScale"-->
      <!--        @before-leave="beforeScale"-->
      <!--      >-->
      <div
        class="collapsible-wrapper"
        v-bind:class="{ collapsed: plotRatioLabel != 'any' || hover }"
      >
        <div class="collapsible">
          <!--        <div v-show="plotRatioLabel == 'any' || hover">-->
          <div class="row justify-center items-baseline row--slide-transition">
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
        </div>
      </div>
      <!--      </transition>-->
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
import { defineComponent, computed, watch, ref } from 'vue';
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
    const hover = ref(false);
    return {
      hover,
      crossSection,
      plotRatioLabel,
      isShowingHelpForInputWithUnits,
      flexRowTitleStyle,
      relativePlotSize,
      // maxComputeTime,
      plotXSideResolution,
      plotYSideResolution,
      nearFieldPlane,
      beforeScale(el: Element) {
        let hel = el as HTMLElement;
        console.log(hel.offsetHeight);
        // for (let child in hel.childNodes) {
        //   console.log(child);
        // }

        // Array.prototype.reduce.call(
        //   hel.childNodes,
        //   function (p: number, c: HTMLElement) {
        //     return p + (c.offsetHeight || 0);
        //   },
        //   0
        // )
      },
    };
  },
});
</script>

<style scoped>
.scale-y-enter-from,
.scale-y-leave-to {
  max-height: 0;
  transform: scaleY(0);
  transform-origin: top;
}
.scale-y-enter-to,
.scale-y-leave-from {
  max-height: 5000px;
  transform: scaleY(1);
  transform-origin: top;
}
.scale-y-enter-active {
  transition: max-height 0.5s cubic-bezier(1, 0.01, 1, 0.01);
}
.scale-y-leave-active {
  transition: max-height 0.5s cubic-bezier(0.01, 1, 0.01, 1);
}

/* based on https://stackoverflow.com/a/43965099 */
.collapsible-wrapper {
  display: flex;
  overflow: hidden;
}
.collapsible-wrapper:after {
  content: '';
  height: 50px;
  transition: height 0.3s linear, max-height 0s 0.3s linear;
  max-height: 0px;
  opacity: 0;
}
.collapsible {
  transition: margin-bottom 0.3s cubic-bezier(0, 0, 0, 1),
    opacity 0.3s cubic-bezier(0.9, 0, 1, 1);
  margin-bottom: 0;
  max-height: 1000000px;
  opacity: 1;
}
.collapsible-wrapper.collapsed > .collapsible {
  margin-bottom: -2000px;
  transition: margin-bottom 0.3s cubic-bezier(1, 0, 1, 1), visibility 0s 0.3s,
    max-height 0s 0.3s, opacity 0.3s;
  visibility: hidden;
  max-height: 0;
  opacity: 0;
}
.collapsible-wrapper.collapsed:after {
  height: 0;
  opacity: 0;
  transition: height 0.3s linear;
  max-height: 50px;
}
</style>
