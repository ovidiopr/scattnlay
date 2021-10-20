export interface guiRuntimeStateInterface {
  isShowingHelpForInputWithUnits: boolean
  units: string
  sourceUnits: string
  isSourceSameUnits: boolean
  activatedMaterials: string[]
}

function state(): guiRuntimeStateInterface {
  return {
    isShowingHelpForInputWithUnits: true,
    units: 'nm',
    sourceUnits: 'nm',
    isSourceSameUnits: true,
    activatedMaterials: [
        // 'PEC',
      'nk-constant']
  }
}

export default state;
