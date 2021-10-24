import { cloneDeep } from 'lodash'
import Spline from 'cubic-spline-ts'
import nmiejs from 'src/nmiejs.js';

export interface layer {
  layerWidth: number
  materialName: string
  n: number
  k: number
  nSpline: Spline|undefined
  kSpline: Spline|undefined
}

export interface simulationSetup {
  hostIndex: number
  fromWL: number; toWL:number; pointsWL:number
  layers: layer[]
  numberOfModesToPlot: number
}

export interface simulationSetupStateInterface {
  library: Map<string,simulationSetup>
  gui: simulationSetup
  current: simulationSetup
}

function setupFactory(hostIndex = 1,
                      fromWL = 300, toWL=1000, pointsWL=101,
                      layers = [
                        {layerWidth:100, n:4, k:0.01,
                          materialName:'nk-constant',
                          nSpline:undefined, kSpline:undefined
                        },
                      ],
                      numberOfModesToPlot = 4
                     ):simulationSetup {
  return {hostIndex:hostIndex,
    fromWL:fromWL, toWL:toWL, pointsWL:pointsWL,
    layers: cloneDeep(layers),
    numberOfModesToPlot: numberOfModesToPlot,
  }
}

function state(): simulationSetupStateInterface {
  // Test nmiejs if working
  void (async () => {
    const module = await nmiejs();
    const nmie = new module.nmie();
    nmie.ClearTarget();
    const R = 100.0;
    const reN = 4.0;
    const imN = 0.01;
    nmie.AddTargetLayerReIm(R, reN, imN)
    nmie.SetModeNmaxAndType(-1, -1);
    const WL = 800;
    nmie.SetWavelength(WL);
    nmie.RunMieCalculation();
    console.log(nmie.GetQsca());
    // outer_arc_points, radius_points, from_Rho, to_Rho,
    // from_Theta, to_Theta, from_Phi, to_Phi, isIgnoreAvailableNmax
    nmie.RunFieldCalculationPolar(2, 2,
        0.1, 1.5, 0, 3.1415, 0, 3.1415,
        0);
    console.log('Field Eabs:', nmie.GetFieldEabs());
  })();

  const gui = setupFactory()
  const current = cloneDeep(gui)
  const library = new Map<string,simulationSetup>()
  library.set('default', cloneDeep(gui))
  return {
    library,
    gui, // simulation setup config as shown in GUI
    current // simulation setup used for the latest simulation
  }
};

export default state;
