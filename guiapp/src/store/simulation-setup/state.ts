import { cloneDeep } from 'lodash'
import Spline from 'cubic-spline-ts'

// All numbers with units (e.g. size, radius, wavelength, e.g.) are given in nanometers.

export interface material {
  name: string
  spectrumRangeStart:number|undefined
  spectrumRangeEnd:number|undefined
  nSpline: Spline|undefined
  kSpline: Spline|undefined
  isPlot: boolean|undefined
}

export interface layer {
  layerWidth: number
  material: material
  n: number
  k: number
}

export enum nearFieldType {
  Ek = 'Ek',
  Hk = 'Hk',
  both = 'both'
}

export interface nearFieldSetup {
  atWL:number
  relativePlotSize: number
  plotSideResolution: number
  crossSection: nearFieldType
  maxComputeTime: number //in seconds
}

export interface simulationSetup {
  hostIndex: number
  fromWL: number; toWL:number; pointsWL:number
  layers: layer[]
  numberOfModesToPlot: number
  plotLabel: string

  nearFieldSetup: nearFieldSetup,
  farFieldWL:number
}

export interface simulationSetupStateInterface {
  library: Map<string,simulationSetup>
  gui: simulationSetup
  current: simulationSetup
  nmie: import('src/nmiejs').nmie_class|undefined
  isNmieLoaded: boolean
  isNmieRunning: boolean
  nmieTotalRunTime: number
}

function setupFactory(hostIndex = 1,
                      fromWL = 400, toWL=1000, pointsWL=201,
                      layers = [
                        {layerWidth:100, n:4, k:0.01,
                          material: {name:'nk-constant',
                            spectrumRangeStart:undefined, spectrumRangeEnd: undefined,
                            nSpline: undefined, kSpline: undefined, isPlot:false},
                        },
                      ],
                      numberOfModesToPlot = 4,
                      plotLabel = '',

                      nearFieldSetup = {
                        atWL: 619,
                        relativePlotSize: 2,
                        plotSideResolution: 64,
                        crossSection: nearFieldType.both,
                        maxComputeTime: 5 //in seconds
                      },

                      farFieldWL = 619,

):simulationSetup {
  return {hostIndex:hostIndex,
    fromWL:fromWL, toWL:toWL, pointsWL:pointsWL,
    layers: cloneDeep(layers),
    numberOfModesToPlot: numberOfModesToPlot,
    plotLabel: plotLabel,

    nearFieldSetup:nearFieldSetup,

    farFieldWL:farFieldWL,
  }
}

function state(): simulationSetupStateInterface {

  const gui = setupFactory()
  const current = cloneDeep(gui)
  const library = new Map<string,simulationSetup>()
  library.set('default', cloneDeep(gui))
  const nmie = undefined
  const isNmieLoaded = false
  const isNmieRunning = false
  const nmieTotalRunTime = 0

  return {
    library,
    gui, // simulation setup config as shown in GUI
    current, // simulation setup used for the latest simulation
    nmie,
    isNmieLoaded, isNmieRunning, nmieTotalRunTime
  }
};

export default state;
