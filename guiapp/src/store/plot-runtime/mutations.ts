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

  setQscaPlotToggle (state: prsi, val: boolean) {state.isPlotQsca = val},
  setQabsPlotToggle (state: prsi, val: boolean) {state.isPlotQabs = val},
  setQextPlotToggle (state: prsi, val: boolean) {state.isPlotQext = val},
  setQscaTotalPlotToggle (state: prsi, val: boolean) {state.isPlotQscaTotal = val},
  setQabsTotalPlotToggle (state: prsi, val: boolean) {state.isPlotQabsTotal = val},
  setQextTotalPlotToggle (state: prsi, val: boolean) {state.isPlotQextTotal = val},
  setIsRemovePlots (state:prsi, val:boolean) {state.isRemovePlots = val},
  setCommonLabel  (state:prsi, val:string) {state.commonLabel = val},

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

  updateNumberOfPlotsFromPreviousSimulations(state: prsi) {
    state.numberOfPlotsFromPreviousSimulations = state.spectraPlot.data.length
  },
  updateSpectraPlot (state: prsi) {
    if (state.isRemovePlots) state.numberOfPlotsFromPreviousSimulations = 0
    state.spectraPlot.data.length = state.numberOfPlotsFromPreviousSimulations

    const label:string = state.commonLabel
    if (state.isPlotQscaTotal) {
      const traceQsca: Partial<Data> = {
        x: state.WLs,
        y: state.Qsca,
        type: 'scatter',
        name: 'Qsca '+label
      }
      state.spectraPlot.data.push(traceQsca)
    }

    if (state.isPlotQabsTotal) {
      const traceQabs: Partial<Data> = {
        x: state.WLs,
        y: state.Qabs,
        type: 'scatter',
        name: 'Qabs '+label
      }
      state.spectraPlot.data.push(traceQabs)
    }

    if (state.isPlotQextTotal) {
      const traceQext: Partial<Data> = {
        x: state.WLs,
        y: state.Qext,
        type: 'scatter',
        name: 'Qext '+label
      }
      state.spectraPlot.data.push(traceQext)
    }



  },

}

export default mutation
