import { MutationTree } from 'vuex';
import {
  simulationSetupStateInterface as sssi,
  simulationSetup,
  layer,
  nearFieldPlane,
  nearFieldSetup,
} from './state';
import { cloneDeep } from 'lodash';
import { markRaw } from 'vue';

const mutation: MutationTree<sssi> = {
  setNmies(
    state: sssi,
    newVal: {
      spectrum: import('src/nmiejs').nmie_class;
      nearField: import('src/nmiejs').nmie_class;
      farField: import('src/nmiejs').nmie_class;
    }
  ) {
    state.nmies.spectrum.instance = markRaw(newVal.spectrum);
    state.nmies.nearField.instance = markRaw(newVal.nearField);
    state.nmies.farField.instance = markRaw(newVal.farField);
  },
  markNmieAsStarted(state: sssi) {
    state.nmies.spectrum.isNmieRunning = true;
  },
  markNmieAsFinished(state: sssi) {
    state.nmies.spectrum.isNmieRunning = false;
  },
  setNmieTotalRunTime(state: sssi, val: number) {
    state.nmies.spectrum.nmieTotalRunTime = val;
  },

  markNmieNearFieldAsStarted(state: sssi) {
    state.nmies.nearField.isNmieRunning = true;
  },
  markNmieNearFieldAsFinished(state: sssi) {
    state.nmies.nearField.isNmieRunning = false;
  },
  setNmieNearFieldTotalRunTime(state: sssi, val: number) {
    state.nmies.nearField.nmieTotalRunTime = val;
  },

  markNmieFarFieldAsStarted(state: sssi) {
    state.nmies.farField.isNmieRunning = true;
  },
  markNmieFarFieldAsFinished(state: sssi) {
    state.nmies.farField.isNmieRunning = false;
  },
  setNmieFarFieldTotalRunTime(state: sssi, val: number) {
    state.nmies.farField.nmieTotalRunTime = val;
  },

  setCurrentState(state: sssi, newVal: simulationSetup) {
    state.current = cloneDeep(newVal);
  },

  copySetupFromGuiToCurrent(state: sssi) {
    state.current = cloneDeep(state.gui);
  },

  // Mutations for simulation setup as represented\set in GUI
  setGuiState(state: sssi, newVal: simulationSetup) {
    state.gui = cloneDeep(newVal);
  },
  setLayers(state: sssi, newVal: layer[]) {
    state.gui.layers = cloneDeep(newVal);
  },
  setHostIndex(state: sssi, val: number) {
    state.gui.hostIndex = val;
  },
  setFromWL(state: sssi, val: number) {
    state.gui.fromWL = val;
  },
  setToWL(state: sssi, val: number) {
    state.gui.toWL = val;
  },
  setPointsWL(state: sssi, val: number) {
    state.gui.pointsWL = val;
  },
  setPlotLabel(state: sssi, val: string) {
    state.gui.plotLabel = val;
  },
  setNumberOfModesToPlot(state: sssi, val: number) {
    state.gui.numberOfModesToPlot = val;
  },

  setNearFieldSetup(state: sssi, val: nearFieldSetup) {
    state.gui.nearFieldSetup = cloneDeep(val);
  },
  setNearFieldWL(state: sssi, val: number) {
    state.gui.nearFieldSetup.atWL = val;
  },
  setNearFieldRelativePlotSize(state: sssi, val: number) {
    state.gui.nearFieldSetup.relativePlotSize = val;
  },
  setNearFieldAtRelativeX0(state: sssi, val: number) {
    if (Math.abs(val) < 1e-15) val = 0;
    state.gui.nearFieldSetup.atRelativeX0 = val;
  },
  setNearFieldAtRelativeY0(state: sssi, val: number) {
    if (Math.abs(val) < 1e-15) val = 0;
    state.gui.nearFieldSetup.atRelativeY0 = val;
  },
  setNearFieldAtRelativeZ0(state: sssi, val: number) {
    if (Math.abs(val) < 1e-15) val = 0;
    state.gui.nearFieldSetup.atRelativeZ0 = val;
  },
  setNearFieldPlotXSideResolution(state: sssi, val: number) {
    state.gui.nearFieldSetup.plotXSideResolution = val;
  },
  setNearFieldPlotYSideResolution(state: sssi, val: number) {
    state.gui.nearFieldSetup.plotYSideResolution = val;
  },
  setNearFieldCrossSection(state: sssi, val: nearFieldPlane) {
    state.gui.nearFieldSetup.crossSection = val;
  },
  // setNearFieldMaxComputeTime     (state: sssi, val: number)         {state.gui.nearFieldSetup.maxComputeTime     = val},

  setFarFieldWL(state: sssi, val: number) {
    state.gui.farFieldWL = val;
  },
};

export default mutation;
