<template>
  <q-page class="column q-px-md">
    <!--    <img src="https://www.pikpng.com/pngl/m/8-80773_new-content-coming-soon-website-under-construction-banner.png" alt="New Content Coming Soon - Website Under Construction">-->
    <q-img :src="img_path" />
    <div class="q-ma-sm" />
    Under construction!

    <div class="q-ma-sm" />
    <GetSourceParameters />
    <GetHostIndex />
    <div class="col-auto">
      Input result: {{ $store.state.simulationSetup.gui.fromWL }}
    </div>
    <!-- place QPageScroller at end of page -->
    <div class="q-ma-lg" />
    <q-page-scroller
      position="bottom-right"
      :scroll-offset="150"
      :offset="[18, 18]"
    >
      <q-btn size="xs" fab icon="keyboard_arrow_up" color="primary" />
    </q-page-scroller>
  </q-page>
</template>

<script lang="ts">
import { defineComponent } from 'vue';
import GetHostIndex from 'components/spectrum/GetHostIndex.vue';
import GetSourceParameters from 'components/spectrum/GetSourceParameters.vue';
// import GetParticleParameters from 'components/GetParticleParameters.vue'
// import { useStore } from 'src/store'

export default defineComponent({
  name: 'FarFieldPage',
  components: {
    GetHostIndex,
    GetSourceParameters,
    // GetParticleParameters
  },
  setup() {
    let img_path = 'underconstruction.jpg';
    if (process.env.publicPath)
      img_path = process.env.publicPath.toString() + 'underconstruction.jpg';
    // const $store = useStore()
    return { img_path };
  },
});
</script>
