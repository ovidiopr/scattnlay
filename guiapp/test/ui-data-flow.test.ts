import { describe, it, expect, beforeEach } from 'vitest';
import { createStore } from 'vuex';
import plotRuntime from '../src/store/plot-runtime';
import simulationSetup from '../src/store/simulation-setup';

describe('Plotly Data Flow', () => {
  let store: any;

  beforeEach(() => {
    store = createStore({ modules: { plotRuntime, simulationSetup } });
  });

  it('should format spectra results into Plotly traces', () => {
    // Mock successful simulation data
    store.commit('plotRuntime/setQ', {
      WLs: [400, 500],
      Qsca: [1.2, 1.5],
      Qabs: [0.1, 0.2],
      Qext: [1.3, 1.7],
      Qsca_n: [[[]], [[]]], // Simplified
      Qabs_n: [[[]], [[]]],
      Qext_n: [[[]], [[]]]
    });

    store.commit('plotRuntime/updateSpectrumPlots');
    
    const chartData = store.state.plotRuntime.spectrumPlots.data;
    expect(chartData.length).toBeGreaterThan(0);
    expect(chartData[0].x).toBeDefined();
    expect(chartData[0].y[0]).toBeCloseTo(1.2);
  });
});
