<template>
  <div class="col-12 q-pa-md">
    <q-table
          :rows="activatedMaterials"
          :columns="columns"
          :rows-per-page-options="[0]"
          hide-pagination
          dense
          title="Available materials"
          title-class="text-h6"
          row-key="name"
      >

      <template #header="props">
        <q-tr :props="props">
          <q-th auto-width/>
          <q-th auto-width class="text-left"> Label </q-th>
          <q-th auto-width class="text-left"> Range </q-th>
          <q-th class="text-left"> Interpolator </q-th>
        </q-tr>
      </template>

      <template #body="props">
        <q-tr
            v-if="props.row.name!='link' && props.row.name!='nk-constant' && props.row.name!='PEC'"
            :props="props"
        >
          <q-td auto-width>
            <q-tooltip anchor="top start" self="bottom start" >
              Delete from simulation</q-tooltip>
            <q-btn size="sm" padding="5px" color="primary" round dense icon="delete"
                   @click="deleteFromSimulation(props.row.name)"/>
          </q-td>


          <q-td auto-width class="">
            {{composeLabelFromPageData(props.row.name)}}
          </q-td>

          <q-td auto-width>
            <ShowSpectrumRange
                :spectrum-range-start="props.row.spectrumRangeStart"
                :spectrum-range-end="props.row.spectrumRangeEnd"
            />
          </q-td>

          <q-td class="">
            <span v-if="props.row.nSpline && props.row.kSpline">
              <q-icon size='sm' color="green" name="done" />
            </span>
            <span v-else>
              <q-icon size='xs' color="red" name="do_not_disturb" />
            </span>
          </q-td>

        </q-tr>
      </template>

      </q-table>
  </div>
</template>

<script lang="ts">
import {
  computed,
  defineComponent,
  reactive,
  watch
} from 'vue'
import { useStore } from 'src/store'
import { composeLabelFromPageData } from 'components/utils'
import { cloneDeep } from 'lodash'
import ShowSpectrumRange from 'components/ShowSpectrumRange.vue'

export default defineComponent({

  name: 'MaterialsActivated',
  components: {ShowSpectrumRange},

  setup: function () {
    const $store = useStore()

    const fromWavelengthStore = computed(()=>$store.state.simulationSetup.gui.fromWL)
    const toWavelengthStore = computed(()=>$store.state.simulationSetup.gui.toWL)

    const columns = [
      // do not change the order without updating the template
      {name: 'spectrumRangeStart', label: 'RangeStart', field: 'spectrumRangeStart'},
      {name: 'spectrumRangeEnd', label: 'RangeEnd', field: 'spectrumRangeEnd'},
      {name: 'name', label: 'name', field: 'name'},
      {name: 'nSpline', nSpline: 'nSpline', field: 'nSpline'},
      {name: 'kSpline', kSpline: 'kSpline', field: 'kSpline'},
    ]

    function getReactiveActivatedMaterials() {
      return reactive( cloneDeep($store.state.guiRuntime.activatedMaterials) )
    }
    let activatedMaterials = getReactiveActivatedMaterials()
    watch($store.state.guiRuntime.activatedMaterials, ()=>{
      // to keep reactivity for local activatedMaterials array: remove all old value, add updated values
      activatedMaterials.splice(0, activatedMaterials.length, ...$store.state.guiRuntime.activatedMaterials)
    })

    return {columns, activatedMaterials,
      fromWavelengthStore, toWavelengthStore,
      composeLabelFromPageData,
      deleteFromSimulation(name:string) {
        $store.commit('guiRuntime/deleteMaterial', name)
      }
    }
  },
})
</script>
