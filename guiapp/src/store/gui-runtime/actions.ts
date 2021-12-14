import { ActionTree } from 'vuex'
import { StateInterface } from '../index'
import { guiRuntimeStateInterface } from './state'
import { composeLabelFromPageData } from 'components/utils'
import { load } from 'js-yaml'

async function loadMaterialData(filename:string){ /* eslint-disable */
  // TODO enable eslint, which is disabled now due to unknown result type of load() from js-yaml
  // TODO implement formulas
  // TODO implement multiple '- type:' in one file (e.g. it can be 'tabulated n' and 'tabulated k'

  let Ag_data

  try {
    const response = await fetch('refractiveindex.info-database/database/data/'+filename)
    const Ag_data = await response.text()

    const doc = await load(Ag_data) as any
    console.log(doc)
    if (doc.DATA[0].type == 'tabulated nk' || doc.DATA[0].type == 'tabulated n') {
      const csv = doc.DATA[0].data
      const rows:string[] = csv.split("\n")
      let data =  rows.map(row =>row.split(" "))
      console.log(data)
      data.pop()
      let data_num = data.map(elem => elem.map(elem2 => parseFloat(elem2)))

      function transpose(array:number[][]) {
        return array[0].map((col, i) => array.map(row => row[i]));
      }

      let data_columns = transpose(data_num)
      // Convert from default refractiveindex.info mkm to nm

      for (let i=0; i<data_columns[0].length; i++)
        data_columns[0][i] *= 1000
      console.log(data_columns)
      return data_columns
    }
  } catch (e) {
    console.log(e)
  }
}


const actions: ActionTree<guiRuntimeStateInterface, StateInterface> = {
  activateMaterial({commit,state}/* context */, filepath:string) {
    // console.log(composeLabelFromPageData(filepath))
    let data_columns = loadMaterialData(filepath)
    console.log(data_columns)
    return {state, filepath}
  },

}

export default actions
