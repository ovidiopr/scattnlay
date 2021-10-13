import { cloneDeep } from 'lodash'
export interface simulationSetup {
  hostIndex: number,
  fromWL: number, toWL:number, pointsWL:number
}

export interface simulationSetupStateInterface {
  library: Map<string,simulationSetup>;
  gui: simulationSetup;
  current: simulationSetup;
}

function setupFactory(hostIndex = 1,
                      fromWL = 300, toWL=1000, pointsWL=100
                     ):simulationSetup {
  return {hostIndex:hostIndex,
    fromWL:fromWL, toWL:toWL, pointsWL:pointsWL }
}

function state(): simulationSetupStateInterface {
  const gui = setupFactory()
  const current = cloneDeep(gui)
  const library = new Map<string,simulationSetup>()
  library.set('default', cloneDeep(gui))
  return {
    library,
    gui,
    current
  }
};

export default state;
