import { cloneDeep } from 'lodash';
import Spline from 'cubic-spline-ts';

// All numbers with units (e.g. size, radius, wavelength, e.g.) are given in nanometers.

export interface material {
  name: string;
  spectrumRangeStart: number | undefined;
  spectrumRangeEnd: number | undefined;
  nSpline: Spline | undefined;
  kSpline: Spline | undefined;
  isPlot: boolean | undefined;
}

export interface layer {
  layerWidth: number;
  material: material;
  n: number;
  k: number;
}

export enum nearFieldPlane {
  Ek = 0,
  Hk,
  EH = 2,
}

export interface nearFieldSetup {
  atWL: number;
  relativePlotSize: number;
  // X0, Y0, and Z0 refers to simulation
  atRelativeX0: number;
  atRelativeY0: number;
  atRelativeZ0: number;
  // X and Y refers to plot
  plotXSideResolution: number;
  plotYSideResolution: number;
  crossSection: nearFieldPlane;
  // maxComputeTime: number //in seconds
}

export interface simulationSetup {
  hostIndex: number;
  fromWL: number;
  toWL: number;
  pointsWL: number;
  layers: layer[];
  numberOfModesToPlot: number;
  plotLabel: string;

  nearFieldSetup: nearFieldSetup;
  farFieldWL: number;
}

export interface nmieModule {
  //nmie instance
  instance: import('src/nmiejs').nmie_class | undefined;
  isNmieRunning: boolean;
  nmieTotalRunTime: number;
}

export interface simulationSetupStateInterface {
  library: Map<string, simulationSetup>;
  gui: simulationSetup;
  current: simulationSetup;
  nmies: { spectrum: nmieModule; nearField: nmieModule; farField: nmieModule };
}

function setupFactory(
  hostIndex = 1,
  fromWL = 400,
  toWL = 1000,
  pointsWL = 201,
  layers = [
    {
      layerWidth: 100,
      n: 4,
      k: 0.01,
      material: {
        name: 'nk-constant',
        spectrumRangeStart: undefined,
        spectrumRangeEnd: undefined,
        nSpline: undefined,
        kSpline: undefined,
        isPlot: false,
      },
    },
  ],
  numberOfModesToPlot = 3,
  plotLabel = '',

  nearFieldSetup = {
    atWL: 619.3885178,
    relativePlotSize: 2,
    atRelativeX0: 0,
    atRelativeY0: 0,
    atRelativeZ0: 0,
    plotXSideResolution: 64,
    plotYSideResolution: 64,
    crossSection: nearFieldPlane.Ek,
    // maxComputeTime: 5 //in seconds
  },

  farFieldWL = 619
): simulationSetup {
  return {
    hostIndex: hostIndex,
    fromWL: fromWL,
    toWL: toWL,
    pointsWL: pointsWL,
    layers: cloneDeep(layers),
    numberOfModesToPlot: numberOfModesToPlot,
    plotLabel: plotLabel,

    nearFieldSetup: nearFieldSetup,

    farFieldWL: farFieldWL,
  };
}

function state(): simulationSetupStateInterface {
  const gui = setupFactory();
  const current = cloneDeep(gui);
  const library = new Map<string, simulationSetup>();
  library.set('default', cloneDeep(gui));
  const nmies = {
    spectrum: {
      instance: undefined,
      isNmieRunning: false,
      nmieTotalRunTime: 0,
    },
    nearField: {
      instance: undefined,
      isNmieRunning: false,
      nmieTotalRunTime: 0,
    },
    farField: {
      instance: undefined,
      isNmieRunning: false,
      nmieTotalRunTime: 0,
    },
  };
  return {
    library,
    gui, // simulation setup config as shown in GUI
    current, // simulation setup used for the latest simulation
    nmies,
  };
}

export default state;
