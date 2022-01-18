import { material } from 'src/store/simulation-setup/state';

// All numbers with units (e.g. size, radius, wavelength, e.g.) are given in nanometers.

export interface guiRuntimeStateInterface {
  isShowingHelpForInputWithUnits: boolean;
  units: string;
  sourceUnits: string;
  isSourceSameUnits: boolean;
  activatedMaterials: material[];
  safeFromWL: number;
  safeToWL: number;
  isSaveWithPythonScript: boolean;
  isAutoRefineNearField: boolean;
  isSquareNearField: boolean;
  nearFieldZoom: { fromX: number; toX: number; fromY: number; toY: number };
}

function state(): guiRuntimeStateInterface {
  return {
    isShowingHelpForInputWithUnits: true,
    units: 'nm',
    sourceUnits: 'nm',
    isSourceSameUnits: true,
    safeFromWL: 0,
    safeToWL: 1e300,
    isSaveWithPythonScript: true,
    nearFieldZoom: { fromX: 0, toX: 1, fromY: 0, toY: 1 },
    isAutoRefineNearField: true,
    isSquareNearField: true,
    activatedMaterials: [
      // 'PEC',
      {
        name: 'link',
        spectrumRangeStart: 0,
        spectrumRangeEnd: 1e300,
        nSpline: undefined,
        kSpline: undefined,
        isPlot: false,
      },
      {
        name: 'nk-constant',
        spectrumRangeStart: 0,
        spectrumRangeEnd: 1e300,
        nSpline: undefined,
        kSpline: undefined,
        isPlot: false,
      },
    ],
  };
}

export default state;
