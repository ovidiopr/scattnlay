import { MutationTree } from 'vuex';
import { cloneDeep } from 'lodash';
import { AxisType, Data } from 'plotly.js-dist-min';
import { plotRuntimeStateInterface as prsi, spectraData } from './state';
import { getModeName, toUnits } from 'components/utils';

const mutation: MutationTree<prsi> = {
  setNearField(state: prsi, val: Float64Array) {
    state.nearFieldEabs = cloneDeep(val);
  },
  setNearFieldCoords(state: prsi, val: { x: number[]; y: number[] }) {
    state.nearFieldCoordX = cloneDeep(val.x);
    state.nearFieldCoordY = cloneDeep(val.y);
  },
  setNearFieldDataFrom(state: prsi, val: number) {
    state.nearFieldDataFrom = val;
  },
  setNearFieldDataTo(state: prsi, val: number) {
    state.nearFieldDataTo = val;
  },
  setNearFieldLimitFrom(state: prsi, val: number) {
    state.nearFieldLimitFrom = val;
  },
  setNearFieldLimitTo(state: prsi, val: number) {
    state.nearFieldLimitTo = val;
  },

  setQ(state: prsi, val: spectraData) {
    state.WLs = cloneDeep(val.WLs);
    state.Qsca = cloneDeep(val.Qsca);
    state.Qabs = cloneDeep(val.Qabs);
    state.Qext = cloneDeep(val.Qext);
    state.Qsca_n = cloneDeep(val.Qsca_n);
    state.Qabs_n = cloneDeep(val.Qabs_n);
    state.Qext_n = cloneDeep(val.Qext_n);
  },

  setWLsInUnits(state: prsi, sourceUnits: string) {
    const converted: number[] = [];
    for (const WL of state.WLs) converted.push(toUnits(WL, sourceUnits));
    state.WLsInUnits = converted; //assign it once to avoid multiple reactivity updates of the spectra plot
  },

  setQscaPlotToggle(state: prsi, val: boolean) {
    state.isPlotQsca = val;
  },
  setQabsPlotToggle(state: prsi, val: boolean) {
    state.isPlotQabs = val;
  },
  setQextPlotToggle(state: prsi, val: boolean) {
    state.isPlotQext = val;
  },
  setQscaTotalPlotToggle(state: prsi, val: boolean) {
    state.isPlotQscaTotal = val;
  },
  setQabsTotalPlotToggle(state: prsi, val: boolean) {
    state.isPlotQabsTotal = val;
  },
  setQextTotalPlotToggle(state: prsi, val: boolean) {
    state.isPlotQextTotal = val;
  },
  setIsRemovePlots(state: prsi, val: boolean) {
    state.isRemovePlots = val;
  },
  setIsLogPlot(state: prsi, val: boolean) {
    state.isLogPlot = val;
  },
  setCommonLabel(state: prsi, val: string) {
    state.commonLabel = val;
  },

  resizeSelectorIsPlotMode(state: prsi, val: number) {
    while (state.isPlotModeE.length > val) state.isPlotModeE.pop();
    while (state.isPlotModeH.length > val) state.isPlotModeH.pop();
    while (state.isPlotModeE.length < val) state.isPlotModeE.push(false);
    while (state.isPlotModeH.length < val) state.isPlotModeH.push(false);
  },
  setIsPlotModeE(state: prsi, val: boolean[]) {
    for (let i = 0; i < val.length; ++i) state.isPlotModeE[i] = val[i];
  },
  setIsPlotModeH(state: prsi, val: boolean[]) {
    for (let i = 0; i < val.length; ++i) state.isPlotModeH[i] = val[i];
  },

  updateXAxisTitle(state: prsi, val: string) {
    if (state.spectrumPlots.layout.xaxis)
      state.spectrumPlots.layout.xaxis.title = val;
  },

  updateNumberOfPlotsFromPreviousSimulations(state: prsi) {
    state.numberOfPlotsFromPreviousSimulations =
      state.spectrumPlots.data.length;
  },
  updateSpectrumPlots(state: prsi) {
    if (state.isRemovePlots) state.numberOfPlotsFromPreviousSimulations = 0;
    state.spectrumPlots.data.length =
      state.numberOfPlotsFromPreviousSimulations;
    let logState: AxisType | undefined = undefined;
    if (state.isLogPlot) logState = 'log';
    if (state.spectrumPlots.layout.yaxis)
      state.spectrumPlots.layout.yaxis.type = logState;

    const label: string = state.commonLabel;
    if (state.isPlotQscaTotal) {
      const traceQsca: Partial<Data> = {
        x: state.WLsInUnits,
        y: state.Qsca,
        type: 'scatter',
        name: 'Qsca ' + label,
      };
      state.spectrumPlots.data.push(traceQsca);
    }

    if (state.isPlotQabsTotal) {
      const traceQabs: Partial<Data> = {
        x: state.WLsInUnits,
        y: state.Qabs,
        type: 'scatter',
        name: 'Qabs ' + label,
      };
      state.spectrumPlots.data.push(traceQabs);
    }

    if (state.isPlotQextTotal) {
      const traceQext: Partial<Data> = {
        x: state.WLsInUnits,
        y: state.Qext,
        type: 'scatter',
        name: 'Qext ' + label,
      };
      state.spectrumPlots.data.push(traceQext);
    }

    const typeNames = ['E', 'H'];
    const totalEvaluatedModes = state.Qsca_n[0].length;

    for (let modeType = 0; modeType < 2; ++modeType) {
      for (let mode_n = 0; mode_n < totalEvaluatedModes; ++mode_n) {
        const isPlotMode =
          modeType === 0 ? state.isPlotModeE : state.isPlotModeH;
        if (!isPlotMode[mode_n]) continue;
        if (state.isPlotQsca) {
          const traceQsca: Partial<Data> = {
            x: state.WLsInUnits,
            y: state.Qsca_n[modeType][mode_n],
            type: 'scatter',
            name:
              'Qsca ' +
              typeNames[modeType] +
              ' ' +
              getModeName(mode_n + 1) +
              ' ' +
              label,
          };
          state.spectrumPlots.data.push(traceQsca);
        }
        if (state.isPlotQabs) {
          const traceQabs: Partial<Data> = {
            x: state.WLsInUnits,
            y: state.Qabs_n[modeType][mode_n],
            type: 'scatter',
            name:
              'Qabs ' +
              typeNames[modeType] +
              ' ' +
              getModeName(mode_n + 1) +
              ' ' +
              label,
          };
          state.spectrumPlots.data.push(traceQabs);
        }
        if (state.isPlotQext) {
          const traceQext: Partial<Data> = {
            x: state.WLsInUnits,
            y: state.Qext_n[modeType][mode_n],
            type: 'scatter',
            name:
              'Qext ' +
              typeNames[modeType] +
              ' ' +
              getModeName(mode_n + 1) +
              ' ' +
              label,
          };
          state.spectrumPlots.data.push(traceQext);
        }
      }
    }
  },
};

export default mutation;
