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
  setQscaPlotToggle (state: prsi, val: boolean) {state.isPlotQsca = val},
  setQabsPlotToggle (state: prsi, val: boolean) {state.isPlotQabs = val},
  setQextPlotToggle (state: prsi, val: boolean) {state.isPlotQext = val},

  resizeSelectorIsPlotMode (state:prsi, val:number) {
    while (state.isPlotModeE.length > val) state.isPlotModeE.pop();
    while (state.isPlotModeH.length > val) state.isPlotModeH.pop();
    while (state.isPlotModeE.length < val) state.isPlotModeE.push(false);
    while (state.isPlotModeH.length < val) state.isPlotModeH.push(false);
  },

  setIsPlotModeE (state:prsi, val:boolean[]) {
    for (let i = 0; i< val.length; ++i) state.isPlotModeE[i] = val[i]
  },
  setIsPlotModeH (state:prsi, val:boolean[]) {
    for (let i = 0; i< val.length; ++i) state.isPlotModeH[i] = val[i]
  },

  setIsAppendPlots (state:prsi, val:boolean) {state.isAppendPlots = val},

}

export default mutation
