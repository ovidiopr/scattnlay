<template>
  <div class="row items-baseline">
    <div class="col-xs-12 col-sm-auto text-center q-px-md q-py-sm">
      <div :style="flexRowTitleStyle">
        |&thinsp;𝐸&thinsp;|&emsp13;∕&emsp13;|&thinsp;𝐸𝜊&thinsp;|
      </div>
    </div>
    <div class="col-xs-grow col-sm">
      <div class="row justify-xs-center justify-sm-start items-center">
        <div class="col-auto">
          <input-with-units
            v-model:input-result="limitFrom"
            v-model:is-showing-help="isShowingHelpForInputWithUnits"
            :initial-expression="limitFrom.toString()"
            title="from"
            units=""
          />
        </div>
        <div class="col-auto">
          <input-with-units
            v-model:input-result="limitTo"
            v-model:is-showing-help="isShowingHelpForInputWithUnits"
            :initial-expression="limitTo.toString()"
            title="to"
            units=""
          />
        </div>
        <div class="col-auto q-pa-xs">
          <q-btn
            no-caps
            flat
            icon="restart_alt"
            color="primary"
            label="reset"
            @click="resetLimits()"
          />
        </div>
      </div>
    </div>
  </div>
  <div class="q-ma-xs" />
  <div class="row items-baseline">
    <div class="col-xs-12 col-sm-auto text-center q-px-md q-py-sm">
      <div :style="flexRowTitleStyle"></div>
    </div>
    <div class="col-xs-grow col-sm">
      <div class="row justify-xs-center justify-sm-start items-center">
        <div class="col-auto">
          <q-checkbox v-model="isLogColorscale" size="sm">
            use log scale
          </q-checkbox>
        </div>
        <div class="col-auto q-mx-sm">
          <q-select
            v-model="colorscale"
            :options="colorscaleOptions"
            options-dense
            outlined
            dense
          />
        </div>
        <div class="col-auto">
          <q-checkbox v-model="moreColors" size="sm">
            use more colors
          </q-checkbox>
          <span class="q-px-sm">
            <q-tooltip anchor="top middle" self="center middle">
              more on color preception
              <q-icon name="launch" />
            </q-tooltip>
            <a
              href="https://www.semanticscholar.org/paper/Why-We-Use-Bad-Color-Maps-and-What-You-Can-Do-About-Moreland/028423bb486b2a963ee6a330ea5cceab467a5349"
            >
              <q-icon name="o_info" size="sm" />
            </a>
          </span>
        </div>
      </div>
    </div>
  </div>
</template>

<script lang="ts">
import { defineComponent, watch, computed, ref } from 'vue';
import { useStore } from 'src/store';
import InputWithUnits from 'components/InputWithUnits.vue';
import { flexRowTitleStyle } from 'components/config';

export default defineComponent({
  name: 'GetNearFieldColorScale',
  components: { InputWithUnits },

  setup() {
    const $store = useStore();

    const colorscale = computed({
      get: () => $store.state.guiRuntime.colorscale,
      set: (val) => $store.commit('guiRuntime/setColorscale', val),
    });

    const isLogColorscale = computed({
      get: () => $store.state.guiRuntime.isLogColorscale,
      set: (val) => $store.commit('guiRuntime/setIsLogColorscale', val),
    });

    const moreColors = ref(false);
    const basicColors = ['Jet', 'Portland', 'Hot', 'Greys', 'Greens'];
    const additionalColors = [
      'YlOrRd',
      'YlGnBu',
      'RdBu',
      'Electric',
      'Earth',
      'Picnic',
      'Bluered',
      'Blackbody',
    ];

    const colorscaleOptions = computed(() => {
      if (moreColors.value) return [...basicColors, ...additionalColors];
      return basicColors;
    });

    const isShowingHelpForInputWithUnits = computed({
      get: () => $store.state.guiRuntime.isShowingHelpForInputWithUnits,
      set: (val) =>
        $store.commit('guiRuntime/setIsShowingHelpForInputWithUnits', val),
    });

    const dataFrom = computed(() => $store.state.plotRuntime.nearFieldDataFrom);
    const dataTo = computed(() => $store.state.plotRuntime.nearFieldDataTo);
    const limitFrom = computed({
      get: () => $store.state.plotRuntime.nearFieldLimitFrom,
      set: (val) => $store.commit('plotRuntime/setNearFieldLimitFrom', val),
    });
    const limitTo = computed({
      get: () => $store.state.plotRuntime.nearFieldLimitTo,
      set: (val) => $store.commit('plotRuntime/setNearFieldLimitTo', val),
    });

    limitFrom.value = dataFrom.value;
    limitTo.value = dataTo.value;

    watch(dataFrom, () => {
      if (dataFrom.value > limitFrom.value) limitFrom.value = dataFrom.value;
    });

    watch(dataTo, () => {
      if (dataTo.value < limitTo.value) limitTo.value = dataTo.value;
    });

    return {
      isShowingHelpForInputWithUnits,
      flexRowTitleStyle,
      limitFrom,
      limitTo,
      resetLimits() {
        limitTo.value = dataTo.value;
        limitFrom.value = dataFrom.value;
      },
      colorscale,
      moreColors,
      colorscaleOptions,
      isLogColorscale,
    };
  },
});
</script>
