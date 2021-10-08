export interface guiRuntimeStateInterface {
  isShowingHelpForInputWithUnits: boolean;
}

function state(): guiRuntimeStateInterface {
  return {
    isShowingHelpForInputWithUnits: false
  }
};

export default state;
