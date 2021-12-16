<template>
  <div class="col-12 q-pa-md">
    <q-table
          :rows="$store.state.guiRuntime.activatedMaterials"
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
          <q-th auto-width> Plot </q-th>
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

          <q-td auto-width>
            <q-checkbox :model-value="props.row.isPlot" size="sm" color="primary" dense
                        @click="$store.commit('guiRuntime/toggleIsPlot', props.row.name )"
            />
          </q-td>


          <q-td auto-width class=""
                @click="$store.commit('guiRuntime/toggleIsPlot', props.row.name )"
          >
            {{composeLabelFromPageData(props.row.name)}}
          </q-td>

          <q-td auto-width
                @click="$store.commit('guiRuntime/toggleIsPlot', props.row.name )"
          >
            <ShowSpectrumRange
                :spectrum-range-start="props.row.spectrumRangeStart"
                :spectrum-range-end="props.row.spectrumRangeEnd"
            />
          </q-td>

          <q-td class=""
                @click="$store.commit('guiRuntime/toggleIsPlot', props.row.name )"
          >
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
} from 'vue'
import { useStore } from 'src/store'
import { composeLabelFromPageData } from 'components/utils'
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

    return {columns,
      fromWavelengthStore, toWavelengthStore,
      composeLabelFromPageData,
      deleteFromSimulation(name:string) {
        $store.commit('guiRuntime/deleteMaterial', name)
      }
    }
  },
})
</script>
