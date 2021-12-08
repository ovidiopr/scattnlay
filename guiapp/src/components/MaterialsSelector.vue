<template>
  <div class="q-pa-md">

    <q-table
          :rows="rows"
          :columns="columns"
          :filter="filter"
          :loading="loading"
          :rows-per-page-options="[3, 12, 24, 0]"
          dense
          title="Available materials"
          title-class="text-h6"
          row-key="id"
      >

        <template #top>
          <div>
            <q-tooltip anchor="top end" self="center middle" >
              Using a copy of RefractiveIndex.info. Analytical models are not implemented.
            </q-tooltip>

            <q-icon size='sm' name="o_info" />
          </div>
          <q-input v-model="filter" dense  debounce="300" color="primary" >
            <template #append>
              <q-icon name="search" />
            </template>
            <template v-if="filter" #after>
              <q-btn  flat round dense icon="cancel" @click="filter=''"/>
            </template>

          </q-input>
          <q-space />
        </template>

      <template #header="props">
        <q-tr :props="props">
          <q-th auto-width />
          <q-th class="text-left"> {{ props.cols[1].label }} </q-th>
          <q-th class="text-left"> {{ props.cols[2].label }} </q-th>
          <q-th class="text-left"> {{ props.cols[3].label }} </q-th>
        </q-tr>
      </template>

      <template #body="props">
        <q-tr :props="props">
          <q-td auto-width>
            <q-tooltip anchor="top end" self="center middle" >
              Add to simulation</q-tooltip>
            <q-btn size="sm" color="accent" round dense icon="add" @click="addToSimulation(props.row.id)"/>
<!--            <q-btn size="sm" color="accent" round dense @click="props.expand = !props.expand" :icon="props.expand ? 'remove' : 'add'" />-->
          </q-td>
          <q-td auto-width> {{ props.row.bookDivider }} </q-td>
          <q-td auto-width>
            <!-- eslint-disable-next-line vue/no-v-html -->
            <div v-html="props.row.bookName"/>
            <q-tooltip anchor="top middle" self="center middle" >
              {{ props.row.shelfDivider }}</q-tooltip>
          </q-td>
          <q-td> {{ props.row.pageName }} </q-td>
        </q-tr>
      </template>

      <template #pagination="scope">
        <q-btn
            v-if="scope.pagesNumber > 2"
            icon="first_page"
            color="grey-8"
            round
            dense
            flat
            :disable="scope.isFirstPage"
            @click="scope.firstPage"
        />

        <q-btn
            icon="chevron_left"
            color="grey-8"
            round
            dense
            flat
            :disable="scope.isFirstPage"
            @click="scope.prevPage"
        />
        {{scope.pagination.page}} of {{scope.pagesNumber}}

        <q-btn
            icon="chevron_right"
            color="grey-8"
            round
            dense
            flat
            :disable="scope.isLastPage"
            @click="scope.nextPage"
        />

        <q-btn
            v-if="scope.pagesNumber > 2"
            icon="last_page"
            color="grey-8"
            round
            dense
            flat
            :disable="scope.isLastPage"
            @click="scope.lastPage"
        />
      </template>

      </q-table>
  </div>
</template>

<script lang="ts">
import {
  // computed,
  defineComponent,
    reactive,
    ref
} from 'vue'
// import { useStore } from 'src/store'
import { load } from 'js-yaml'

export default defineComponent({

  name: 'MaterialsSelector',
  components: {},

  setup: function () {
    // const $store = useStore()
    const loading = ref(true)

    const columns = [
      // do not change the order without updating the template
      {name: 'id', label: '', field: 'id'},
      {name: 'bookDivider', label: 'State', field: 'bookDivider'},
      {name: 'bookName', label: 'Material', field: 'bookName'},
      {name: 'pageName', label: 'Authors', field: 'pageName'},
      {name: 'shelfDivider', label: 'Group', field: 'shelfDivider'},
      {name: 'pageData', label: 'file', field: 'pageData'},
    ]

    async function GetPagesFromLib() { /* eslint-disable */
      let rows = []
      // lib has an irregular structure
      const response = await fetch('refractiveindex.info-database/database/library.yml')
      const data = await response.text()
      const lib = await load(data) as any
      let i = 1
      for (const shelf of lib) {
        let shelfDivider = ''
        let bookName = ''
        for (const bookOrDivider of shelf.content) {
          if (bookOrDivider.DIVIDER) {
            shelfDivider = bookOrDivider.DIVIDER
            continue
          } else if (bookOrDivider.name) {

            bookName = bookOrDivider.name
            let bookDivider = ''
            let pageName = ''
            let pageData = ''
            for (const pageOrDivider of bookOrDivider.content) {
              if (pageOrDivider.DIVIDER) {
                bookDivider = pageOrDivider.DIVIDER
                continue
              } else if (pageOrDivider.name) {
                pageName = pageOrDivider.name
                pageData = pageOrDivider.data
                if (bookDivider.includes('Model') || bookDivider.includes('model')
                    || pageName.includes('Model') || pageName.includes('model')) continue

                rows.push({ id: i,
                  shelfDivider: shelfDivider, bookName: bookName,
                  bookDivider: bookDivider, pageName: pageName, pageData: pageData
                })
                i = i + 1
                // console.log (shelfDivider, '==', bookName, '==',bookDivider, '==', pageName, '==', pageData)
                // console.log (shelfDivider, '==', bookName, '==',bookDivider, '==', pageName, '==', pageData)
              } else {
                console.log('###################### Unknown type in pageOrDivider', pageOrDivider)
              }

            }
          } else {
            console.log('###################### Unknown type in bookOrDivider', bookOrDivider)
          }
        }
      }
      return rows
    }


    let rowsProto:{
      shelfDivider: string, bookName: string,
      bookDivider: string, pageName: string, pageData: string
    }[] = []
    const rows= reactive(rowsProto)
    GetPagesFromLib().then(val => {
      rows.push(...val);
      loading.value = false
    })

    // const activatedMaterials = computed(() => $store.state.guiRuntime.activatedMaterials)
    // const columnsQ = computed(()=>{
    //   return [
    //     { name: 'name', label: '',  align: 'left', field: 'name', headerStyle:''},
    //     { name: 'Qsca', label: 'Qsca',  align: 'left', field: 'Qsca', headerStyle:''},
    //     { name: 'Qabs', label: 'Qabs',  align: 'left', field: 'Qabs', headerStyle:''},
    //     { name: 'Qext', label: 'Qext',  align: 'left', field: 'Qext', headerStyle:''},
    //   ]
    // })


    const filter = ref('')
    return {columns, rows, loading, filter,
      addToSimulation(val:number) {
        console.log(rows[val-1])
    }}
  },
})
</script>
e
