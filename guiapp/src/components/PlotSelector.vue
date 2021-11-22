<template>
  <div class="row items-baseline">
    <div class="col-auto q-py-xs q-px-sm">
      <q-checkbox v-model="isRemovePlots" size="sm">
        remove previous spectra
      </q-checkbox>
    </div>
  </div>
  <div class="row items-baseline">
    <div class="col-xs-grow col-sm-auto q-px-sm">
      <q-table
              title="Values to plot"
              title-class="text-h6"
              :rows="rowsQ"
              :columns="columnsQ"
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
                        v-model="props.row[index]"
                        dense
                    />
                  </q-td>
                </template>
              </q-tr>
            </template>
          </q-table>
    </div>
    <div v-if="isPlotQabs||isPlotQsca||isPlotQext"  class="col-xs-grow col-sm-auto q-px-sm">
      <q-table
          title="Modes to plot"
          title-class="text-h6"
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
import { getModeName} from 'components/utils'
import { flexRowTitleStyle,
} from 'components/config'


export default defineComponent({
  name: 'GetPlotSettings',

  setup() {
    const $store = useStore()

    const isRemovePlots = computed({
      get: ()=> $store.state.plotRuntime.isRemovePlots,
      set: val => {
        $store.commit('plotRuntime/setIsRemovePlots', val)
        $store.commit('plotRuntime/updateSpectraPlot')
      }
    })

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

    const isPlotQscaTotal = computed({
      get: () => $store.state.plotRuntime.isPlotQscaTotal,
      set: val => $store.commit('plotRuntime/setQscaTotalPlotToggle', val)
    })
    const isPlotQabsTotal = computed({
      get: () => $store.state.plotRuntime.isPlotQabsTotal,
      set: val => $store.commit('plotRuntime/setQabsTotalPlotToggle', val)
    })
    const isPlotQextTotal = computed({
      get: () => $store.state.plotRuntime.isPlotQextTotal,
      set: val => $store.commit('plotRuntime/setQextTotalPlotToggle', val)
    })

    const columnsQ = computed(()=>{
      return [
        { name: 'name', label: '',  align: 'left', field: 'name', headerStyle:''},
        { name: 'Qsca', label: 'Qsca',  align: 'left', field: 'Qsca', headerStyle:''},
        { name: 'Qabs', label: 'Qabs',  align: 'left', field: 'Qabs', headerStyle:''},
        { name: 'Qext', label: 'Qext',  align: 'left', field: 'Qext', headerStyle:''},
      ]
    })

    const rowsQ_store = computed({
    get: ()=>[
        { name: 'total', Qsca: isPlotQscaTotal.value, Qabs: isPlotQabsTotal.value, Qext: isPlotQextTotal.value},
        { name: 'modes', Qsca: isPlotQsca.value, Qabs: isPlotQabs.value, Qext: isPlotQext.value},
    ],
      set: val =>{
        $store.commit('plotRuntime/setQscaTotalPlotToggle', val[0].Qsca)
        $store.commit('plotRuntime/setQabsTotalPlotToggle', val[0].Qabs)
        $store.commit('plotRuntime/setQextTotalPlotToggle', val[0].Qext)
        $store.commit('plotRuntime/setQscaPlotToggle', val[1].Qsca)
        $store.commit('plotRuntime/setQabsPlotToggle', val[1].Qabs)
        $store.commit('plotRuntime/setQextPlotToggle', val[1].Qext)
        // console.log(val)
      }
    })
    const rowsQ = ref(rowsQ_store.value)
    watch(rowsQ_store, ()=>{
          rowsQ.value = rowsQ_store.value
        },
        { deep: true }
    )
    watch(rowsQ, ()=>{
          rowsQ_store.value = rowsQ.value
        },
        { deep: true }
    )

    const guiNumberOfModes = computed(()=> $store.state.simulationSetup.gui.numberOfModesToPlot)
    const simulatedNumberOfModes = computed(()=> $store.state.simulationSetup.current.numberOfModesToPlot)

    const columns = computed(()=> {
      let columns = [{ name: 'name', label: '',  align: 'left', field: '', headerStyle:''}]
      for (let i=1; i<=guiNumberOfModes.value; ++i) {
        let label_computed = getModeName(i)
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

    $store.commit('plotRuntime/resizeSelectorIsPlotMode', guiNumberOfModes.value)

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

    return {flexRowTitleStyle,
      isRemovePlots,
      isPlotQsca, isPlotQabs, isPlotQext,
      guiNumberOfModes, simulatedNumberOfModes,
      columns, rows,
      columnsQ, rowsQ
      }
  },
})
</script>
