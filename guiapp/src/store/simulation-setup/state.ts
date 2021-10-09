export interface simulationSetupStateInterface {
  isShowingHelpForInputWithUnits: boolean;
}

function state(): simulationSetupStateInterface {
  return {
    isShowingHelpForInputWithUnits: true
  }
};

export default state;
