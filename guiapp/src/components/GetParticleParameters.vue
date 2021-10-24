<template>
  <div class="row items-baseline">
    <div class="col-xs-12 col-sm-auto text-weight-bold text-center q-px-md q-py-sm">
      <div :style="flexRowTitleStyle">
        Spherical particle
      </div>
    </div>
    <div class="col-xs-grow col-sm  q-px-xs">
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

      </div>
    </div>
  </div>
  <div v-for="(layer, index) in layers" :key="index" class="row items-baseline q-py-xs">
    <div class="col-xs-12 col-sm-auto text-center q-px-md q-py-sm">
      <div :style="flexRowTitleStyle">
        {{getLayerTitle(particleType, index)}}
      </div>
    </div>
    <div class="col-xs-grow col-sm">
      <div class="row justify-xs-center justify-sm-start items-baseline">

        <div class="col-auto">
          <q-tooltip
              class="bg-primary shadow-4"
              anchor="center left"
              self="center middle">
            {{getTooltipText(index)}}
          </q-tooltip>
          <input-with-units
            v-model:is-showing-help="isShowingHelpForInputWithUnits"
            :initial-expression="toUnits(layer.layerWidth, units).toString()"
            :units="units"
            :title=" index==0 ? 'R' : 'h' "
            :input-result="toUnits(layer.layerWidth, units)"
            @update:input-result="layer.layerWidth = fromUnits(units,$event)"
          />
        </div>

        <div class="col-auto">
          <q-select
              v-model="layer.materialName"
              :options="activatedMaterials"
              style="width: 17.7em"
              class="q-px-xs"
              dense
              options-dense
              outlined
              options-value="value"
              options-label="label"
          />
        </div>

      </div>
      <div v-if="layer.materialName=='nk-constant'"
               class="row justify-xs-center justify-sm-start items-baseline">

        <div class="col-auto"> <input-with-units
            v-model:input-result="layer.n"
            :initial-expression="layer.n.toString()"
            title="Re(n)"
        /></div>
        <div class="col-auto"> <input-with-units
            v-model:input-result="layer.k"
            :initial-expression="layer.k.toString()"
            title="Im(n)"

        /></div>
      </div>
    </div>
  </div>
</template>

<script lang="ts">
import {
  defineComponent,
  ref,
  reactive,
  computed,
  watch
  } from 'vue'
import { useStore } from 'src/store'
import { flexRowTitleStyle, fromUnits, toUnits } from 'components/utils'
import {cloneDeep} from 'lodash'
import InputWithUnits from 'components/InputWithUnits.vue';

export default defineComponent({
  name: 'GetParticleParameters',
  components: {InputWithUnits},

  setup() {
    const $store = useStore()
    const particleType=ref('bulk')
    const numberOfLayers=ref(1)

    const isShowingHelpForInputWithUnits = computed({
      get: () => $store.state.plotRuntime.isShowingHelpForInputWithUnits,
      set: val => $store.commit('plotRuntime/setIsShowingHelpForInputWithUnits', val)
    })

    const units = computed({
      get: () => $store.state.guiRuntime.units,
      set: val => $store.commit('guiRuntime/setUnits', val)
    })

    const activatedMaterials = computed(() => $store.state.guiRuntime.activatedMaterials)

    function getReactiveLayers() {
      return reactive( cloneDeep($store.state.simulationSetup.gui.layers) )
    }

    let layers = getReactiveLayers()

    watch($store.state.simulationSetup.gui.layers, ()=>{
      layers = getReactiveLayers()
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
              materialName: 'nk-constant',
              n: 4.0,
              k: 0.01,
              nSpline: undefined,
              kSpline: undefined,
            }
        );
      }
    })

    function getLayerTitle (particleType:string, index:number) {
      if (particleType=='core-shell'  && index == 0) return 'core'
      if (particleType=='core-shell'  && index == 1) return 'shell'
      if (particleType=='multilayer') return 'layer '+ (index+1).toString()
      return 'bulk'
    }

    function getTooltipText(index:number) {
      if (index == 0) return 'radius'
      return 'thickness'

    }

    return { flexRowTitleStyle, particleType,
      numberOfLayers, layers, getLayerTitle, getTooltipText,
      units, toUnits, fromUnits, isShowingHelpForInputWithUnits,
      activatedMaterials
      }
  },
})
</script>
