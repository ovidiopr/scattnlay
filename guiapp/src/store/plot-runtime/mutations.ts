import { MutationTree } from 'vuex'
import { cloneDeep } from 'lodash'
import {AxisType, Data} from 'plotly.js-dist-min'
import { plotRuntimeStateInterface as prsi, spectraData } from './state'
import { getModeName, toUnits } from 'components/utils'
import { material } from 'src/store/simulation-setup/state'


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

  setWLsInUnits (state:prsi, sourceUnits:string) {
    const converted:number[] = []
    for (const WL of state.WLs) converted.push(toUnits(WL, sourceUnits))
    state.WLsInUnits = converted //assign it once to avoid multiple reactivity updates of the spectra plot
  },

  setQscaPlotToggle (state: prsi, val: boolean) {state.isPlotQsca = val},
  setQabsPlotToggle (state: prsi, val: boolean) {state.isPlotQabs = val},
  setQextPlotToggle (state: prsi, val: boolean) {state.isPlotQext = val},
  setQscaTotalPlotToggle (state: prsi, val: boolean) {state.isPlotQscaTotal = val},
  setQabsTotalPlotToggle (state: prsi, val: boolean) {state.isPlotQabsTotal = val},
  setQextTotalPlotToggle (state: prsi, val: boolean) {state.isPlotQextTotal = val},
  setIsRemovePlots (state:prsi, val:boolean) {state.isRemovePlots = val},
  setIsLogPlot (state:prsi, val:boolean) {state.isLogPlot = val},
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

  updateXAxisTitle (state:prsi, val:string) {
    if (state.spectrumPlots.layout.xaxis) state.spectrumPlots.layout.xaxis.title = val
    if (state.materialPlots.layout.xaxis) state.materialPlots.layout.xaxis.title = val
  },

  updateNumberOfPlotsFromPreviousSimulations(state: prsi) {
    state.numberOfPlotsFromPreviousSimulations = state.spectrumPlots.data.length
  },
  updateSpectrumPlots (state: prsi) {
    if (state.isRemovePlots) state.numberOfPlotsFromPreviousSimulations = 0
    state.spectrumPlots.data.length = state.numberOfPlotsFromPreviousSimulations
    let logState:AxisType|undefined = undefined
    if (state.isLogPlot) logState = 'log'
    if (state.spectrumPlots.layout.yaxis) state.spectrumPlots.layout.yaxis.type = logState

    const label:string = state.commonLabel
    if (state.isPlotQscaTotal) {
      const traceQsca: Partial<Data> = {
        x: state.WLsInUnits,
        y: state.Qsca,
        type: 'scatter',
        name: 'Qsca '+label
      }
      state.spectrumPlots.data.push(traceQsca)
    }

    if (state.isPlotQabsTotal) {
      const traceQabs: Partial<Data> = {
        x: state.WLsInUnits,
        y: state.Qabs,
        type: 'scatter',
        name: 'Qabs '+label
      }
      state.spectrumPlots.data.push(traceQabs)
    }

    if (state.isPlotQextTotal) {
      const traceQext: Partial<Data> = {
        x: state.WLsInUnits,
        y: state.Qext,
        type: 'scatter',
        name: 'Qext '+label
      }
      state.spectrumPlots.data.push(traceQext)
    }

    const typeNames = ['E', 'H']
    const totalEvaluatedModes = state.Qsca_n[0].length

    for (let modeType = 0; modeType < 2; ++modeType) {
      for (let mode_n = 0; mode_n < totalEvaluatedModes; ++mode_n) {
        const isPlotMode = modeType === 0 ? state.isPlotModeE : state.isPlotModeH;
        if (!isPlotMode[mode_n]) continue;
        if (state.isPlotQsca) {
          const traceQsca: Partial<Data> = {
            x: state.WLsInUnits,
            y: state.Qsca_n[modeType][mode_n],
            type: 'scatter',
            name: 'Qsca ' + typeNames[modeType] + ' ' + getModeName(mode_n + 1)+' '+label
          };
          state.spectrumPlots.data.push(traceQsca);
        }
        if (state.isPlotQabs) {
          const traceQabs: Partial<Data> = {
            x: state.WLsInUnits,
            y: state.Qabs_n[modeType][mode_n],
            type: 'scatter',
            name: 'Qabs ' + typeNames[modeType] + ' ' + getModeName(mode_n + 1)+' '+label
          };
          state.spectrumPlots.data.push(traceQabs);
        }
        if (state.isPlotQext) {
          const traceQext: Partial<Data> = {
            x: state.WLsInUnits,
            y: state.Qext_n[modeType][mode_n],
            type: 'scatter',
            name: 'Qext ' + typeNames[modeType] + ' ' + getModeName(mode_n + 1)+' '+label
          };
          state.spectrumPlots.data.push(traceQext);
        }
      }
    }
  },

  // TODO move state.materialPlots and updateMaterialPlots to PlotMaterials.vue. At the moment
  // *.vue have incomplete TypeScript support (e.g. no way to define
  //     const traceDataReN: Partial<Data>
  // )
  updateMaterialPlots(state:prsi, val:{activatedMaterials:material[], sourceUnits:string,
    fromWL:number, toWL:number, pointsWL:number, plotRange: string,
    isPlotReN: boolean, isPlotImN: boolean, isPlotInterpolation: boolean })
  {
    state.materialPlots.data.length = 0
    for (const material of val.activatedMaterials) {
      if (!material.isPlot) continue

      if (!material.nSpline) continue
      if (val.isPlotReN) {
        const traceDataReN: Partial<Data> = {
          x: material.nSpline.xs.map(x => toUnits(x, val.sourceUnits)),
          y: material.nSpline.ys,
          mode: 'markers',
          marker:{size:5},
          type: 'scatter',
          name: 'Re(n) ' + material.name + ' data'
        }
        state.materialPlots.data.push(traceDataReN)
      }

      if (!material.kSpline) continue
      if (val.isPlotImN) {
        const traceDataImN: Partial<Data> = {
          x: material.kSpline.xs.map(x => toUnits(x, val.sourceUnits)),
          y: material.kSpline.ys,
          mode: 'markers',
          marker:{size:5},
          type: 'scatter',
          name: 'Im(n) ' + material.name + ' data'
        }
        state.materialPlots.data.push(traceDataImN)
      }

      if (val.isPlotInterpolation) {
        let fromWL = val.fromWL
        let toWL = val.toWL
        let pointsWL = val.pointsWL-1
        if (state.materialPlots.layout.xaxis) state.materialPlots.layout.xaxis.range = [fromWL, toWL]
        if (val.plotRange == 'material data') {
          fromWL = material.nSpline.xs[0]
          toWL = material.nSpline.xs[material.nSpline.xs.length-1]
          pointsWL = 1000
          if (state.materialPlots.layout.xaxis) state.materialPlots.layout.xaxis.range = undefined
        }
        const stepWL = (toWL-fromWL)/pointsWL

        const WLs:number[] =[]
        const nSpline:number[] = []
        const kSpline:number[] = []
        for (let i=0; i<pointsWL; ++i) {
          WLs.push(fromWL+i*stepWL)
          nSpline.push(material.nSpline.at(fromWL+i*stepWL))
          kSpline.push(material.kSpline.at(fromWL+i*stepWL))
        }
        WLs.push(toWL)
        nSpline.push(material.nSpline.at(toWL))
        kSpline.push(material.kSpline.at(toWL))

        if (val.isPlotReN) {
          if (!material.nSpline) continue
          const traceDataReNi: Partial<Data> = {
            x: WLs.map(x => toUnits(x, val.sourceUnits)),
            y: nSpline,
            type: 'scatter',
            name: 'Re(n) ' + material.name + ' interpolate'
          }
          state.materialPlots.data.push(traceDataReNi)
        }

        if (val.isPlotImN) {
          if (!material.nSpline) continue
          const traceDataImNi: Partial<Data> = {
            x: WLs.map(x => toUnits(x, val.sourceUnits)),
            y: kSpline,
            type: 'scatter',
            name: 'Im(n) ' + material.name + ' interpolate'
          }
          state.materialPlots.data.push(traceDataImNi)
        }

      }
    }
  }

}

export default mutation
