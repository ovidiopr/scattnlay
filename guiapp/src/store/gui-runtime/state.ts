export interface guiRuntimeStateInterface {
  isShowingHelpForInputWithUnits: boolean;
}

function state(): guiRuntimeStateInterface {
  return {
    isShowingHelpForInputWithUnits: true
  }
};

export default state;
