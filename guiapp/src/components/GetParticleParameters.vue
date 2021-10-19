<template>
  <div class="row items-baseline">
    <div class="col-xs-12 col-sm-auto text-weight-bold text-center q-px-md q-py-sm">
      <div :style="flexRowTitleStyle">
        Spherical particle
      </div>
    </div>
    <div class="col-xs-grow col-sm">
      <div class="row justify-xs-center justify-sm-start items-baseline">

        <div class="col-auto" >
          <div class="q-gutter-x-md">
            <q-radio v-model="particleType" dense size='sm' val="bulk" label="bulk" />
            <q-radio v-model="particleType" dense size='sm' val="core-shell" label="core-shell" />
          </div>

        </div>

        <div class="col-auto q-px-md">
          <div class="row justify-xs-center justify-sm-start items-baseline">
            <div class="col-auto" >
              <q-radio v-model="particleType" dense size='sm' val="multilayer" label="multilayer" />
            </div>
            <div class="col-auto" >
              <q-input
                  v-model.number="numberOfLayers"
                  :disable="particleType!='multilayer'"
                  :outlined="particleType=='multilayer'"
                  :filled="particleType!='multilayer'"
                  type="number"
                  class="q-px-sm"
                  min=1
                  max=10
                  step=1
                  dense
                  style="width: 6em"
              />
            </div>
          </div>
        </div>

        <div class="col-auto" >
          numOfLayers: {{numberOfLayers}}
          <br> layers: {{layers}}
        </div>
      </div>
    </div>
  </div>
</template>

<script lang="ts">
import {
  defineComponent,
    ref,
    reactive,
  // computed,
  watch

  } from 'vue'
import { useStore } from 'src/store'
import { flexRowTitleStyle } from 'components/utils'
import {cloneDeep} from 'lodash'


export default defineComponent({

  name: 'GetHostIndex',
  components: {},

  setup() {
    const $store = useStore()

    // const isSourceSameUnits = computed({
    //   get: () => $store.state.guiRuntime.isSourceSameUnits,
    //   set: val => $store.commit('guiRuntime/setIsSourceSameUnits', val)
    // })

    const particleType=ref('bulk')
    const numberOfLayers=ref(1)

    // const layers = computed({
    //   get: () => cloneDeep($store.state.simulationSetup.gui.layers),
    //   set: /*val*/ () => {
    //     // // Doesn't work?
    //     // $store.commit('simulationSetup/setLayers', val)
    //   }
    // })

    let layers = reactive(cloneDeep($store.state.simulationSetup.gui.layers))

    watch($store.state.simulationSetup.gui.layers, ()=>{
      layers = reactive(cloneDeep($store.state.simulationSetup.gui.layers))
    })

    watch(layers, ()=>{
      $store.commit('simulationSetup/setLayers', layers)
    })

    watch(particleType, ()=>{
      if (particleType.value=='bulk') numberOfLayers.value = 1
      if (particleType.value=='core-shell') numberOfLayers.value = 2
      if (particleType.value=='multilayer') numberOfLayers.value = 3
    })

    watch(numberOfLayers, ()=>{
      numberOfLayers.value = parseInt(numberOfLayers.value.toString())
      if (isNaN(numberOfLayers.value)) numberOfLayers.value = 3
      if (numberOfLayers.value < 1) numberOfLayers.value = 1
      if (numberOfLayers.value > 10) numberOfLayers.value = 10

      while (numberOfLayers.value < layers.length) {
        layers.pop();
      }
      let coreR = layers[0].layerWidth;
      while (numberOfLayers.value > layers.length) {
        // r_prev = r_prev*1.1;
        layers.push(
            {
              layerWidth: coreR*0.1,
              materialName: 'nk',
              n: 4.0,
              k: 0.01,
              nSpline: undefined,
              kSpline: undefined,
            }
        );
      }

    })

    return { flexRowTitleStyle, particleType, numberOfLayers, layers
      }
  },
})
</script>
