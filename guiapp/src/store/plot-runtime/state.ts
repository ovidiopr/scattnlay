import { Data, Layout, Config } from 'plotly.js-dist-min';

export interface plotlyChart {
  data: Data[];
  layout: Partial<Layout>;
  config: Partial<Config> | undefined;
}

export interface spectraData {
  WLs: number[];
  Qsca: number[];
  Qabs: number[];
  Qext: number[];
  Qsca_n: number[][][];
  Qabs_n: number[][][];
  Qext_n: number[][][];
}

export interface plotRuntimeStateInterface {
  // Near field
  nearFieldEabs: Float64Array | undefined;
  nearFieldCoordX: number[];
  nearFieldCoordY: number[];
  nearFieldDataFrom: number;
  nearFieldDataTo: number;
  nearFieldLimitFrom: number;
  nearFieldLimitTo: number;
  // Spectra plots
  WLs: number[];
  WLsInUnits: number[];
  Qsca: number[];
  Qabs: number[];
  Qext: number[];
  Qsca_n: number[][][];
  Qabs_n: number[][][];
  Qext_n: number[][][];
  spectrumPlots: plotlyChart;
  isPlotQsca: boolean;
  isPlotQabs: boolean;
  isPlotQext: boolean;
  isPlotQscaTotal: boolean;
  isPlotQabsTotal: boolean;
  isPlotQextTotal: boolean;
  isPlotModeE: boolean[];
  isPlotModeH: boolean[];
  isRemovePlots: boolean;
  isLogPlot: boolean;
  numberOfPlotsFromPreviousSimulations: number;
  commonLabel: string;
}

function state(): plotRuntimeStateInterface {
  const nearFieldEabs = undefined;
  const nearFieldCoordX: number[] = [];
  const nearFieldCoordY: number[] = [];
  const nearFieldDataFrom = 0;
  const nearFieldDataTo = 1e300;
  const nearFieldLimitFrom = 0;
  const nearFieldLimitTo = 1e300;

  const WLs: number[] = [];
  const WLsInUnits: number[] = [];
  const Qsca: number[] = [],
    Qabs: number[] = [],
    Qext: number[] = [];
  const Qsca_n: number[][][] = [[], []];
  const Qabs_n: number[][][] = [[], []];
  const Qext_n: number[][][] = [[], []];
  const isPlotQsca = true;
  const isPlotQabs = false;
  const isPlotQext = false;
  const isPlotQscaTotal = true;
  const isPlotQabsTotal = false;
  const isPlotQextTotal = false;
  const isPlotModeE: boolean[] = [true];
  const isPlotModeH: boolean[] = [true];
  const numberOfPlotsFromPreviousSimulations = 0;
  const commonLabel = '';

  const spectrumPlots: plotlyChart = {
    data: [],
    layout: {
      margin: {
        l: 0,
        r: 40,
        b: 50,
        t: 0,
      },
      // paper_bgcolor: '#7f7f7f',
      // plot_bgcolor: '#c7c7c7',
      // title: 'reactive charts',
      xaxis: {
        // will be set on mount
        title: 'Wavelength, nm',
      },
      yaxis: {
        title: 'Normalized cross-sections',
      },
      showlegend: true,
      legend: {
        orientation: 'h',
        x: -0.15,
        y: 1.12,
      },
    },
    config: {
      responsive: true,
      // showEditInChartStudio: true,
      displaylogo: false,
    },
  };

  const isRemovePlots = true;
  const isLogPlot = false;

  return {
    nearFieldEabs,
    nearFieldCoordX,
    nearFieldCoordY,
    nearFieldDataFrom,
    nearFieldDataTo,
    nearFieldLimitFrom,
    nearFieldLimitTo,
    WLs,
    WLsInUnits,
    Qsca,
    Qabs,
    Qext,
    Qsca_n,
    Qabs_n,
    Qext_n,
    spectrumPlots,
    isPlotQsca,
    isPlotQabs,
    isPlotQext,
    isPlotQscaTotal,
    isPlotQabsTotal,
    isPlotQextTotal,
    isPlotModeE,
    isPlotModeH,
    isRemovePlots,
    isLogPlot,
    numberOfPlotsFromPreviousSimulations,
    commonLabel,
  };
}

export default state;
