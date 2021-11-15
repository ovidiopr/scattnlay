<template>
  <div class="col-auto">
    <div class="row justify-xs-center justify-sm-start items-baseline">
      <q-toggle v-model="isPlotQsca" dense class="q-ml-md">Qsca</q-toggle>
      <q-toggle v-model="isPlotQabs" dense class="q-ml-md">Qabs</q-toggle>
      <q-toggle v-model="isPlotQext" dense class="q-ml-md">Qext</q-toggle>
    </div>
    <div class="q-ma-sm"/>
    <q-toggle
        v-model="rows[0]['1']"
    />

    <div class="row justify-xs-center justify-sm-start items-baseline">
      <q-table
          :rows="rows"
          :columns="columns"
          hide-bottom
          dense
          flat
          row-key="name"
      >
        <template #body="props">
          <q-tr :props="props">
            <q-th key="name" :props="props">
              {{ props.row.name }}
            </q-th>
            <q-td v-for="(val, index) in props.row" :key="index" :props="props">
              <div v-if="index!='name'" >
                <q-toggle
                    v-model="props.row[index.toString()]"
                />
                <!--                    @update:input-result="layer.layerWidth = fromUnits(units,$event)"-->
                {{ props.row[index.toString()]}}
              </div>
            </q-td>
          </q-tr>
        </template>
      </q-table>
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
import { cloneDeep } from 'lodash'

export default defineComponent({
  name: 'GetPlotSettings',

  setup() {
    const $store = useStore()
    const isPlotQsca = computed({
      get: () => $store.state.plotRuntime.isPlotQsca,
      set: val => $store.commit('plotRuntime/setQscaPlotToggle', val)
    })
    const isPlotQabs = computed({
      get: () => $store.state.plotRuntime.isPlotQabs,
      set: val => $store.commit('plotRuntime/setQabsPlotToggle', val)
    })
    const isPlotQext = computed({
      get: () => $store.state.plotRuntime.isPlotQext,
      set: val => $store.commit('plotRuntime/setQextPlotToggle', val)
    })

    const guiNumberOfModes = computed(()=> $store.state.simulationSetup.gui.numberOfModesToPlot)
    const simulatedNumberOfModes = computed(()=> $store.state.simulationSetup.current.numberOfModesToPlot)

    $store.commit('plotRuntime/resizeIsPlotMode', guiNumberOfModes)


    const columns = computed(()=> {
      let columns = [{ name: 'name', label: '',  align: 'left', field: '' }]
      for (let i=1; i<=guiNumberOfModes.value; ++i) {
        let label_computed = Math.pow(2, i).toString()
        if (i == 1) label_computed = 'dipole'
        if (i == 2) label_computed = 'quadrupole'
        if (i == 3) label_computed = 'octupole'
        columns.push({
          name: i.toString(),
          label: label_computed,
          align:'right',
          field: i.toString()
        })
      }
      return columns
    })

    const rows_store = computed({
      get: ()=> {
        let rows = []
        let row = {name: 'E'}
        for (let i = 1; i <= guiNumberOfModes.value; ++i) {
          // eslint-disable-next-line
          (row as any)[i.toString()] = false
        }
        rows.push(cloneDeep(row))
        row.name = 'H'
        rows.push(cloneDeep(row))
        console.log(rows)
        return rows
      },
      set: val => { console.log(val)}
    })
    const rows = ref(rows_store.value)

    return { isPlotQsca, isPlotQabs, isPlotQext,
      guiNumberOfModes, simulatedNumberOfModes,
      columns, rows
      }
  },
})
</script>
