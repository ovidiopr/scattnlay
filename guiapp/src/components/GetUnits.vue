<template>
  <div class="row items-baseline">
    <div class="col-xs-12 col-sm-auto text-weight-bold text-center q-px-md q-py-sm">
      <div :style="flexRowTitleStyle">
        Units
      </div>
    </div>
    <div class="col-xs-grow col-sm  q-px-xs">
      <div class="row justify-xs-center justify-sm-start items-center">

        <div class="col-auto" >
          <div class="row items-center justify-center">
            <div class="text-right"> length </div>
          <q-select
              v-model="units"
              :options="unitsOptions"
              class="q-pa-xs"
              outlined
              dense
              options-dense
              :style="'width: '+basicSelectorWidthStyle"
              behavior="menu"
          >
          </q-select>
            <div/>
          </div>
        </div>

        <div class="col-auto" >
          <div class="row items-center justify-center">
            <div class="text-right q-pl-md">  plane wave </div>

            <q-select
                v-model="sourceUnits"
                :options="sourceUnitsOptions"
                :style="'width: '+basicSelectorWidthStyle"
                class="q-pa-xs"
                :disable="isSourceSameUnits"
                :outlined="!isSourceSameUnits"
                :filled="isSourceSameUnits"
                dense
                options-dense
                option-value="label"
                option-label="label"
                behavior="menu"
            >
              <template #option="scope">
                <q-item :label="scope.opt.title" dense>
                  <q-item-section class="text-weight-bold">{{ scope.opt.title }}</q-item-section>
                </q-item>
                <template v-for="child in scope.opt.children" :key="child.label">
                  <q-item v-close-popup dense clickable
                          @click="sourceUnits = child"
                  >
                    <q-item-section> <q-item-label class="q-ml-md"> {{child.label}} </q-item-label> </q-item-section>
                  </q-item>
                </template>
              </template>
            </q-select>


            <div class="row q-pa-xs items-center">
              <q-checkbox v-model="isSourceSameUnits" size="sm"/>
              <div>same</div>
            </div>
          </div>
        </div>

<!--        <div class="col-auto" >-->
<!--          units: {{sourceUnits}}-->
<!--          <br> store: {{$store.state.guiRuntime}}-->
<!--        </div>-->
      </div>
    </div>
  </div>
</template>

<script lang="ts">
import {
  defineComponent,
  computed,
  watch
  } from 'vue'
import { useStore } from 'src/store'
import { flexRowTitleStyle, basicSelectorWidthStyle } from 'components/config'

export default defineComponent({

  name: 'GetUnits',
  components: {},

  setup() {
    const unitsOptions = [ 'nm', 'mkm', 'mm', 'cm', 'm']
    const sourceUnitsOptions = [
      { title: 'Frequency',
        children: [{label: 'MHz'},{label: 'GHz'},{label: 'THz'}]},
      { title: 'Energy',
        children: [{label: 'meV'}, {label: 'eV'}]},
      { title: 'Period duration',
        children: [{label: 'ps'}, {label: 'fs'}]}
    ]
    const $store = useStore()

    const isSourceSameUnits = computed({
      get: () => $store.state.guiRuntime.isSourceSameUnits,
      set: val => $store.commit('guiRuntime/setIsSourceSameUnits', val)
    })


    const units = computed({
      get: () => $store.state.guiRuntime.units,
      set: val => $store.commit('guiRuntime/setUnits', val)
    })

    const sourceUnits = computed({
      get: () => {return {label:$store.state.guiRuntime.sourceUnits}},
      set: val => $store.commit('guiRuntime/setSourceUnits', val.label)
    })

    let prevCustomSourceUnits = 'THz'
    watch(isSourceSameUnits, ()=>{
      if (isSourceSameUnits.value) {
        prevCustomSourceUnits = $store.state.guiRuntime.sourceUnits
        $store.commit('guiRuntime/setSourceUnits',units.value)
      }
      else $store.commit('guiRuntime/setSourceUnits',prevCustomSourceUnits)
    })

    watch(units, ()=>{
      if (isSourceSameUnits.value) $store.commit('guiRuntime/setSourceUnits',units.value)
    })


    return { isSourceSameUnits, flexRowTitleStyle, basicSelectorWidthStyle,
      unitsOptions, units,
      sourceUnits, sourceUnitsOptions }
  },
})
</script>
