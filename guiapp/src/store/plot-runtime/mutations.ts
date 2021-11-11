import { MutationTree } from 'vuex'
import { plotRuntimeStateInterface as prsi, spectraData } from './state'
import {cloneDeep} from 'lodash'
import { Data } from 'plotly.js-dist-min'

const mutation: MutationTree<prsi> = {
  setQ (state: prsi, val: spectraData) {
    state.WLs    = cloneDeep(val.WLs)
    state.Qsca   = cloneDeep(val.Qsca)
    state.Qabs   = cloneDeep(val.Qabs)
    state.Qext   = cloneDeep(val.Qext)
    state.Qsca_n = cloneDeep(val.Qsca_n)
    state.Qabs_n = cloneDeep(val.Qabs_n)
    state.Qext_n = cloneDeep(val.Qext_n)
  },

  updateSpectraPlot (state: prsi) {
    const traceQsca:Partial<Data> = {
      x: state.WLs,
      y: state.Qsca,
      type: 'scatter',
      name: 'Qsca'
    }
    state.spectraPlot.data.push(traceQsca)

  },
  setQscaToggle (state: prsi, val: boolean) {state.isQscaToggle = val},
  setQabsToggle (state: prsi, val: boolean) {state.isQabsToggle = val},
  setQextToggle (state: prsi, val: boolean) {state.isQextToggle = val},
}

export default mutation
