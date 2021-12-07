<template>
  <div class="row items-baseline">
    <div class="col-xs-grow col-sm-auto q-px-sm">
      <q-table
          title="Activated materials"
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
  </div>
</template>

<script lang="ts">
import {
  defineComponent,
  computed,
  } from 'vue'
import { useStore } from 'src/store'

export default defineComponent({

  name: 'MaterialsSelector',
  components: {},

  setup() {
    const $store = useStore()

    const activatedMaterials = computed(() => $store.state.guiRuntime.activatedMaterials)
    const columnsQ = computed(()=>{
      return [
        { name: 'name', label: '',  align: 'left', field: 'name', headerStyle:''},
        { name: 'Qsca', label: 'Qsca',  align: 'left', field: 'Qsca', headerStyle:''},
        { name: 'Qabs', label: 'Qabs',  align: 'left', field: 'Qabs', headerStyle:''},
        { name: 'Qext', label: 'Qext',  align: 'left', field: 'Qext', headerStyle:''},
      ]
    })


    return { }
  },
})
</script>
