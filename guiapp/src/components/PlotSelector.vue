<template>
  <div class="col-auto">
    <div class="row justify-xs-center justify-sm-start items-baseline">
      <q-toggle v-model="isPlotQsca" dense class="q-ml-md">Qsca</q-toggle>
      <q-toggle v-model="isPlotQabs" dense class="q-ml-md">Qabs</q-toggle>
      <q-toggle v-model="isPlotQext" dense class="q-ml-md">Qext</q-toggle>
    </div>
    <div class="q-ma-sm"/>
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
            <template v-for="(val, index) in props.row" :key="index" :props="props">
              <q-td v-if="index!='name'">
                <q-checkbox
                    v-model="props.row[index.toString()]"
                    size="sm"
                    dense
                    :disable="index>simulatedNumberOfModes"
                />
              </q-td>
            </template>
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
  ref,
  watch
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

    const columns = computed(()=> {
      let columns = [{ name: 'name', label: '',  align: 'left', field: '', headerStyle:''}]
      for (let i=1; i<=guiNumberOfModes.value; ++i) {
        let label_computed = Math.pow(2, i).toString()
        if (i == 1) label_computed = 'dipole'
        if (i == 2) label_computed = 'quadrupole'
        if (i == 3) label_computed = 'octupole'
        let text_color = ''
        if (i > simulatedNumberOfModes.value ) text_color='color:LightGray'
        columns.push({
          name: i.toString(),
          label: label_computed,
          align:'left',
          field: i.toString(),
          headerStyle: text_color
        })
      }
      return columns
    })

    $store.commit('plotRuntime/resizeIsPlotMode', guiNumberOfModes.value)

    const rows_store = computed({
      get: ()=> {
        let rows = []
        let rowE = {name: 'E'}
        let rowH = {name: 'H'}
        for (let i = 0; i < guiNumberOfModes.value; ++i) {
          // eslint-disable-next-line
          (rowE as any)[(i+1).toString()] = $store.state.plotRuntime.isPlotModeE[i];
          // eslint-disable-next-line
          (rowH as any)[(i+1).toString()] = $store.state.plotRuntime.isPlotModeH[i]
        }
        rows.push(cloneDeep(rowE))
        rows.push(cloneDeep(rowH))
        return rows
      },
      set: val => {
        let isPlotModeE: boolean[] = []
        let isPlotModeH: boolean[] = []
        const rowE = val[0]
        const rowH = val[1]
        for (let i = 0; i < guiNumberOfModes.value; ++i) {
          // eslint-disable-next-line
          isPlotModeE.push((rowE as any)[(i+1).toString()])
          // eslint-disable-next-line
          isPlotModeH.push((rowH as any)[(i+1).toString()])
        }
        $store.commit('plotRuntime/setIsPlotModeE', isPlotModeE)
        $store.commit('plotRuntime/setIsPlotModeH', isPlotModeH)
        }
    })

    const rows = ref(rows_store.value)
    watch(rows_store, ()=>{
      rows.value = rows_store.value
        },
        { deep: true }
    )
    watch(rows, ()=>{
      rows_store.value = rows.value
        },
        { deep: true }
    )

    return { isPlotQsca, isPlotQabs, isPlotQext,
      guiNumberOfModes, simulatedNumberOfModes,
      columns, rows
      }
  },
})
</script>
