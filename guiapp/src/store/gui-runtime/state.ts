export interface guiRuntimeStateInterface {
  isShowingHelpForInputWithUnits: boolean
  units: string
  sourceUnits: string
  isSourceSameUnits: boolean
}

function state(): guiRuntimeStateInterface {
  return {
    isShowingHelpForInputWithUnits: true,
    units: 'nm',
    sourceUnits: 'nm',
    isSourceSameUnits: true
  }
}

export default state;
