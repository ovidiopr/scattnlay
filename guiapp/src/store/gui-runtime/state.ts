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
}

function state(): guiRuntimeStateInterface {
  return {
    isShowingHelpForInputWithUnits: true,
    units: 'nm',
    sourceUnits: 'nm',
    isSourceSameUnits: true,
    safeFromWL: 0,
    safeToWL: 1e300,
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
