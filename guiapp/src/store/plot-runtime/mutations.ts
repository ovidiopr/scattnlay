import { MutationTree } from 'vuex'
import { plotRuntimeStateInterface, spectraData } from './state'
import {cloneDeep} from 'lodash'

const mutation: MutationTree<plotRuntimeStateInterface> = {
  setQ (state: plotRuntimeStateInterface, val: spectraData) {
    state.WLs    = (val.WLs)
    state.Qsca   = (val.Qsca)
    state.Qabs   = (val.Qabs)
    state.Qext   = (val.Qext)
    state.Qsca_n = (val.Qsca_n)
    state.Qabs_n = (val.Qabs_n)
    state.Qext_n = (val.Qext_n)
    // state.WLs    = cloneDeep(val.WLs)
    // state.Qsca   = cloneDeep(val.Qsca)
    // state.Qabs   = cloneDeep(val.Qabs)
    // state.Qext   = cloneDeep(val.Qext)
    // state.Qsca_n = cloneDeep(val.Qsca_n)
    // state.Qabs_n = cloneDeep(val.Qabs_n)
    // state.Qext_n = cloneDeep(val.Qext_n)
  }
}

export default mutation
