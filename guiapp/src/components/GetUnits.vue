<template>
  <div class="row items-baseline">
    <div class="col-xs-12 col-sm-auto text-weight-bold text-center q-px-md q-py-sm">
      <div :style="flexRowTitleStyle">
        Units
      </div>
    </div>
    <div class="col-xs-grow col-sm">
      <div class="row justify-xs-center justify-sm-start items-center">

        <div class="col-auto" >
          <q-select
              class="q-pa-xs"
              outlined
              dense
              options-dense
              style="width: 7em"
              behavior="menu"
              v-model="localUnits"
              :options="unitsOptions"
          >
          </q-select>
        </div>

        <div class="col-auto" >
          <div class="row q-py-xs q-px-md items-center justify-center">
          <div class="q-pr-xs"> Source units </div>

            <q-select
                v-model="localSourceUnits"
                :options="sourceUnitsOptions"
                outlined
                dense
                options-dense
                option-value="label"
                option-label="label"
            >
              <template v-slot:option="scope">
                <q-item :label="scope.opt.title" dense>
                  <q-item-section class="text-weight-bold">{{ scope.opt.title }}</q-item-section>
                </q-item>
                <template v-for="child in scope.opt.children" :key="child.label">
                  <q-item clickable v-close-popup dense
                          @click="localSourceUnits = child"
                  >
                    <q-item-section> <q-item-label v-html="child.label" class="q-ml-md" /> </q-item-section>
                  </q-item>
                </template>
              </template>
            </q-select>


            <div class="row q-py-xs items-center">
              <q-checkbox v-model="isSourceSameUnits" size="sm"/>
              <div>same</div>
            </div>
          </div>
        </div>

        <div class="col-auto" >
          units: {{localSourceUnits}}
          <br> isSource: {{isSourceSameUnits}}

        </div>
      </div>
    </div>
  </div>
</template>

<script lang="ts">
import {
  defineComponent,
  computed,
  ref
  } from 'vue'
import { useStore } from 'src/store'
import { flexRowTitleStyle } from 'components/utils'

export default defineComponent({

  name: 'GetHostIndex',
  components: {},

  setup() {
    const unitsOptions = [ 'nm', 'mkm', 'mm', 'cm', 'm', 'km' ]
    let localUnits = ref('nm')
    const sourceUnitsOptions = [
      { title: 'Frequency',
        children: [{label: 'Hz'},{label: 'kHz'},{label: 'MHz'},{label: 'GHz'},{label: 'THz'}]},
      { title: 'Energy',
        children: [{label: 'meV'}, {label: 'eV'}]},
      { title: 'Period duration',
        children: [{label: 'ps'}, {label: 'fs'}]}
    ]
    let localSourceUnits = ref({ "label": "THz" } )
    const $store = useStore()

    // let isSourceOtherUnits = ref(true)
    const isSourceSameUnits = computed({
      get: () => $store.state.guiRuntime.isSourceSameUnits,
      set: val => $store.commit('guiRuntime/setIsSourceSameUnits', val)
    })


    return { isSourceSameUnits, flexRowTitleStyle,
      unitsOptions, localUnits,
      localSourceUnits, sourceUnitsOptions }
  },
})
</script>
