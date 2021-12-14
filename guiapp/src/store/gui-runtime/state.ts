import { material } from 'src/store/simulation-setup/state'

// All numbers with units (e.g. size, radius, wavelength, e.g.) are given in nanometers.


export interface guiRuntimeStateInterface {
  isShowingHelpForInputWithUnits: boolean
  units: string
  sourceUnits: string
  isSourceSameUnits: boolean
  activatedMaterials: material[]
}

function state(): guiRuntimeStateInterface {
  return {
    isShowingHelpForInputWithUnits: true,
    units: 'nm',
    sourceUnits: 'nm',
    isSourceSameUnits: true,
    activatedMaterials: [
        // 'PEC',
      {name:'link', fileFullPath:'', spectrumRangeStart:0, spectrumRangeEnd:1e300,
        nSpline: undefined, kSpline: undefined},
      {name:'nk-constant', fileFullPath:'', spectrumRangeStart:0, spectrumRangeEnd:1e300,
        nSpline: undefined, kSpline: undefined},
      {name:'Ag_McPeak', fileFullPath:'main/Ag/McPeak.yml', spectrumRangeStart:300, spectrumRangeEnd:1700,
        nSpline: undefined, kSpline: undefined},
      {name:'Au_McPeak', fileFullPath:'main/Au/McPeak.yml', spectrumRangeStart:300, spectrumRangeEnd:1700,
        nSpline: undefined, kSpline: undefined},
    ]
  }
}

export default state;
