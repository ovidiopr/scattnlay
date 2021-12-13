<template>
  <div class="col-12 q-pa-md">
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
          <div class="q-mr-md">
            <q-tooltip anchor="top end" self="center middle" >
              Using a copy of RefractiveIndex.info website.<br> Analytical models are not implemented.
            </q-tooltip>
            <q-icon size='sm' name="o_info" />
          </div>

          <q-input v-model="filter" dense  debounce="200" color="primary" >
          <q-tooltip  v-if="!filter" anchor="top middle" self="center middle">filter by any string</q-tooltip>
            <template #append>
              <q-icon name="filter_alt" />
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
          <q-th auto-width />
          <q-th auto-width class="text-left"> Label </q-th>
          <q-th auto-width class="text-left"> Range </q-th>
          <q-th auto-width class="text-left"> Material </q-th>
          <q-th class="text-left"> Details </q-th>
        </q-tr>
      </template>

      <template #body="props">
        <q-tr :props="props">

          <q-td auto-width>
            <q-tooltip anchor="top end" self="center middle" >
              Add to simulation</q-tooltip>
            <q-btn size="sm" color="primary" round dense icon="add" @click="addToSimulation(props.row.id)"/>
          </q-td>

          <q-td auto-width>
            <q-tooltip anchor="top end" self="center middle" >
              Download *.yml file</q-tooltip>
            <q-btn flat
                   size="md"
                   color="primary"
                   icon="download"
                   padding="xs 2px"
                   @click="downloadPageData(props.row.pageData)"
            />
          </q-td>

          <q-td class="">
            {{composeLabelFromPageData(props.row.pageData)}}
          </q-td>

          <q-td auto-width>
            <q-tooltip
                v-if="props.row.spectrumRangeStart>=fromWavelengthStore/1000 ||
                 props.row.spectrumRangeEnd<=toWavelengthStore/1000"
                anchor="top middle" self="center middle"
                class="bg-red">
              Mismatch with spectrum simulation
            </q-tooltip>
            <span :class="props.row.spectrumRangeStart>=fromWavelengthStore/1000?'text-red':'text-black'">
              {{ props.row.spectrumRangeStart }}
            </span>
            <span v-if="props.row.spectrumRangeStart">&ndash;</span>
            <span :class="props.row.spectrumRangeEnd<=toWavelengthStore/1000?'text-red':'text-black'">
              {{ props.row.spectrumRangeEnd }}
            </span>
            <span v-if="props.row.spectrumRangeStart">&NonBreakingSpace;mkm</span>
          </q-td>

          <q-td auto-width>
            <!-- eslint-disable-next-line vue/no-v-html -->
            <div v-html="props.row.bookName"/>
            <q-tooltip anchor="top middle" self="center middle" >
              {{ props.row.shelfDivider }}</q-tooltip>
          </q-td>

          <q-td> {{ props.row.bookDivider}}<span v-if="props.row.bookDivider">;&nbsp;</span>{{props.row.pageName }} </q-td>

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
  ref
} from 'vue'
import { useStore } from 'src/store'
import { load } from 'js-yaml'
import { composeLabelFromPageData } from 'components/utils'
import { saveAs } from 'file-saver'


export default defineComponent({

  name: 'MaterialsSelector',
  components: {},

  setup: function () {
    const $store = useStore()
    const loading = ref(true)

    const fromWavelengthStore = computed(()=>$store.state.simulationSetup.gui.fromWL)
    const toWavelengthStore = computed(()=>$store.state.simulationSetup.gui.toWL)

    const columns = [
      // do not change the order without updating the template
      {name: 'id', label: '', field: 'id'},
      {name: 'bookDivider', label: 'State', field: 'bookDivider'},
      {name: 'bookName', label: 'Material', field: 'bookName'},
      {name: 'spectrumRangeStart', label: 'RangeStart', field: 'spectrumRangeStart'},
      {name: 'spectrumRangeEnd', label: 'RangeEnd', field: 'spectrumRangeEnd'},
      {name: 'pageName', label: 'Details', field: 'pageName'},
      {name: 'shelfDivider', label: 'Group', field: 'shelfDivider'},
      {name: 'pageData', label: 'file', field: 'pageData'},
    ]

    function GetRange(val:string) {
      if (val.lastIndexOf('µm') == -1) return [val,'']
      const rangeStartPosition = val.slice(0,val.lastIndexOf('µm')-1).lastIndexOf(' ')
      const spectrumRange = val.slice(rangeStartPosition+1,val.lastIndexOf('µm')+2)
      const newPageName = val.replace(spectrumRange,'')
      const spectrumRangeStart = spectrumRange.slice(0, spectrumRange.lastIndexOf('-'))
      const spectrumRangeEnd = spectrumRange.slice(spectrumRange.lastIndexOf('-')+1,
          spectrumRange.lastIndexOf('µm')-1)
      return [newPageName,spectrumRangeStart, spectrumRangeEnd]
    }

    function splitBookName(val:string) {
      if (val.lastIndexOf('(') == -1) return val
      const splitPosition = val.lastIndexOf('(')
      const mainPart = val.slice(0,splitPosition-1)
      const otherPart = val.slice(splitPosition+1, val.length-1)
      return mainPart+'<br><span style="font-size: 0.5em;">'+otherPart+'</span>'
    }

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

            bookName = splitBookName(bookOrDivider.name)
            let bookDivider = ''
            let pageName = ''
            let pageData = ''
            for (const pageOrDivider of bookOrDivider.content) {
              if (pageOrDivider.DIVIDER) {

                bookDivider = pageOrDivider.DIVIDER.replace('Experimental data:','').replace('Experimental data','')
                continue

              } else if (pageOrDivider.name) {

                pageName = pageOrDivider.name
                pageData = pageOrDivider.data
                if (bookDivider.includes('Model') || bookDivider.includes('model')
                    || pageName.includes('Model') || pageName.includes('model')) continue
                const pageNameSplit = GetRange(pageName)
                rows.push({ id: i, shelfDivider: shelfDivider, bookName: bookName, bookDivider: bookDivider,
                  pageName: pageNameSplit[0], spectrumRangeStart: pageNameSplit[1],
                  spectrumRangeEnd: pageNameSplit[2], pageData: pageData })
                i = i + 1

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

    const activatedMaterials = computed(() => $store.state.guiRuntime.activatedMaterials)
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
      fromWavelengthStore, toWavelengthStore,
      composeLabelFromPageData,
      addToSimulation(val:number) {
        console.log(rows[val-1].pageData)
        $store.commit('guiRuntime/activateMaterial', rows[val-1].pageData)
      },
      async downloadPageData(filepath:string) {
        const response = await fetch('refractiveindex.info-database/database/data/'+filepath)
        const data = await response.text()

        const scattnlaySpectra = new Blob([data],
            {type: 'text/plain;charset=utf-8',
              endings: 'native'}  //TODO test if newline is correctly written in Windows, MacOS
        )
        saveAs(scattnlaySpectra, composeLabelFromPageData(filepath));
      }
    }
  },
})
</script>
