import { Data, Layout, Config } from 'plotly.js-dist-min'

export interface plotlyChart {
  data: Data[],
  layout: Partial<Layout>,
  config: Partial<Config>|undefined
}

export interface spectraData {
  WLs: number[]
  Qsca:number[]
  Qabs:number[]
  Qext:number[]
  Qsca_n:number[][][]
  Qabs_n:number[][][]
  Qext_n:number[][][]
}

export interface plotRuntimeStateInterface {
  WLs: number[]
  Qsca:number[]
  Qabs:number[]
  Qext:number[]
  Qsca_n:number[][][]
  Qabs_n:number[][][]
  Qext_n:number[][][]
  spectraPlot: plotlyChart
  isPlotQsca: boolean
  isPlotQabs: boolean
  isPlotQext: boolean
  isPlotModeE: boolean[]
  isPlotModeH: boolean[]
  isAppendPlots: boolean
}

function state(): plotRuntimeStateInterface {
  const WLs:number[] = []
  const Qsca:number[] = [], Qabs:number[] = [], Qext:number[] = []
  const Qsca_n:number[][][] = [[], []]
  const Qabs_n:number[][][] = [[], []]
  const Qext_n:number[][][] = [[], []]
  const isPlotQsca = true
  const isPlotQabs = true
  const isPlotQext = false
  const isPlotModeE:boolean[] = []
  const isPlotModeH:boolean[] = []

  const spectraPlot:plotlyChart = {
    data: [],
    layout: {
      margin: {
        l: 0,
        r: 40,
        b: 50,
        t: 0
      },
      // paper_bgcolor: '#7f7f7f',
      // plot_bgcolor: '#c7c7c7',
      // title: 'reactive charts',
      xaxis: {
        // will be set on mount
        title: ''
      },
      yaxis: {
        title: 'Normalized cross-sections'
      },
      showlegend: true,
      legend: {
        orientation: 'h',
        x: -.1,
        y: 1.05
      },
    },
    config: {responsive: true,
      // showEditInChartStudio: true,
      displaylogo: false}
  }
  const isAppendPlots = false

  return { WLs,
    Qsca, Qabs, Qext, Qsca_n, Qabs_n, Qext_n,
    spectraPlot, isPlotQsca, isPlotQabs, isPlotQext,
    isPlotModeE, isPlotModeH,
    isAppendPlots
  }
}

export default state
