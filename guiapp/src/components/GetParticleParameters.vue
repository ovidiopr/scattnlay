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
                  :max="maxNumberOfLayers"
                  step=1
                  dense
                  :style="'width: '+basicSelectorWidthStyle"
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
              v-model="layer.material"
              :options="activatedMaterials"
              :style="'width: '+basicWidthStyle"
              class="q-px-xs"
              dense
              options-dense
              outlined
              options-value="name"
              options-label="name"
              :display-value="layer.material.name"
          >
            <template #option="scope">
              <span v-if="scope.opt.name =='link'">
                <q-item clickable dense
                        to="/materials"
                >
                  <q-item-section side>
                    <q-icon size="xs" name="settings"/>
                  </q-item-section>
                  <q-item-section>
                    <q-item-label >Manage materials...</q-item-label>
                  </q-item-section>
                </q-item>
                <q-separator inset />
              </span>
              <span v-else>
                <q-item
                    v-bind="scope.itemProps"
                    :disable="(!scope.opt.kSpline || !scope.opt.kSpline)
                     && !(scope.opt.name=='nk-constant' || scope.opt.name=='PEC')"
                >
                  <q-item-section>
                    <q-item-label>{{ scope.opt.name }}</q-item-label>
                    <!--                  <q-item-label caption>{{ scope.opt.description }}</q-item-label>-->
                  </q-item-section>
                </q-item>
              </span>
            </template>
          </q-select>
        </div>

      </div>
      <div v-if="layer.material.name=='nk-constant'"
               class="row justify-xs-center justify-sm-start items-baseline">

        <div class="col-auto"> <input-with-units
            v-model:is-showing-help="isShowingHelpForInputWithUnits"
            v-model:input-result="layer.n"
            :initial-expression="layer.n.toString()"
            title="Re(n)"
        /></div>
        <div class="col-auto"> <input-with-units
            v-model:is-showing-help="isShowingHelpForInputWithUnits"
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
import { fromUnits, isAlmostSame, toUnits } from 'components/utils'
import { flexRowTitleStyle, basicWidthStyle, basicSelectorWidthStyle, maxNumberOfLayers } from 'components/config'
import { cloneDeep } from 'lodash'
import InputWithUnits from 'components/InputWithUnits.vue';

export default defineComponent({
  name: 'GetParticleParameters',
  components: {InputWithUnits},

  setup() {
    const $store = useStore()

    const numberOfLayers=ref($store.state.simulationSetup.gui.layers.length)
    const particleType=ref('bulk')
    if (numberOfLayers.value == 1) particleType.value='bulk'
    if (numberOfLayers.value == 2) particleType.value='core-shell'
    if (numberOfLayers.value == 3) particleType.value='multilayer'
    // Initially we set numberOfLayers and particleType to match layers in store.
    // Later on changes to numberOfLayers and particleType lead to update of layers in store.

    const isShowingHelpForInputWithUnits = computed({
      get: () => $store.state.guiRuntime.isShowingHelpForInputWithUnits,
      set: val => $store.commit('guiRuntime/setIsShowingHelpForInputWithUnits', val)
    })

    const units = computed({
      get: () => $store.state.guiRuntime.units,
      set: val => $store.commit('guiRuntime/setUnits', val)
    })

    const activatedMaterials = computed(() => $store.state.guiRuntime.activatedMaterials)
    // const activatedMaterialsOptions = computed(() => {
    //   let list = activatedMaterials.value.map(x=>{
    //     return {
    //       name: x.name, value: x.name, label: x.name,
    //       spectrumRangeStart:x.spectrumRangeStart, spectrumRangeEnd:x.spectrumRangeEnd,
    //       nSpline: x.nSpline, kSpline: x.kSpline
    //     }
    //   })
    //   console.log(list)
    //   return list
    // })

    let layers = reactive( cloneDeep($store.state.simulationSetup.gui.layers) )

    watch($store.state.simulationSetup.gui.layers,
        ()=>{
          layers.splice(0, layers.length,
              ...($store.state.simulationSetup.gui.layers.map(layer=>cloneDeep(layer)))
          )
          numberOfLayers.value = $store.state.simulationSetup.gui.layers.length
        },
        { deep: true }
    )

    watch(layers,
        ()=>{
          const storeLayers = $store.state.simulationSetup.gui.layers
          for (let i = 0; i<layers.length && i<storeLayers.length; i++) {
            if (isAlmostSame(storeLayers[i].layerWidth, layers[i].layerWidth)) {
              layers[i].layerWidth = storeLayers[i].layerWidth
            }
            // if (layers[i].material.name == 'link') {
            //   layers[i].material.name = storeLayers[i].material.name
            // }
          }
          $store.commit('simulationSetup/setLayers', layers)
        },
        { deep: true }
    )

    watch(particleType, ()=>{
      if (particleType.value=='bulk') numberOfLayers.value = 1
      if (particleType.value=='core-shell') numberOfLayers.value = 2
      if (particleType.value=='multilayer') numberOfLayers.value = 3
    })

    watch(numberOfLayers, ()=>{
      numberOfLayers.value = parseInt(numberOfLayers.value.toString())
      if (isNaN(numberOfLayers.value)) return
      if (numberOfLayers.value < 1) numberOfLayers.value = 1
      if (numberOfLayers.value > maxNumberOfLayers) numberOfLayers.value = maxNumberOfLayers

      while (layers.length > numberOfLayers.value) {
        layers.pop();
      }
      let coreR = layers[0].layerWidth;
      while (layers.length < numberOfLayers.value) {
        // r_prev = r_prev*1.1;
        layers.push(
            {
              layerWidth: coreR*0.1,
              material: {name:'nk-constant', fileFullPath:undefined,
                spectrumRangeStart:undefined, spectrumRangeEnd: undefined,
                nSpline: undefined, kSpline: undefined},
              n: 4.0,
              k: 0.01,
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

    return { flexRowTitleStyle, basicWidthStyle,
      basicSelectorWidthStyle, maxNumberOfLayers,
      particleType,
      numberOfLayers, layers, getLayerTitle, getTooltipText,
      units, toUnits, fromUnits, isShowingHelpForInputWithUnits,
      activatedMaterials
      }
  },
})
</script>
