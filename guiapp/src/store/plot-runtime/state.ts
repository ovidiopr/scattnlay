export interface plotRuntimeStateInterface {
  isShowingHelpForInputWithUnits: boolean;
}

function state(): plotRuntimeStateInterface {
  return {
    isShowingHelpForInputWithUnits: true
  }
};

export default state;
