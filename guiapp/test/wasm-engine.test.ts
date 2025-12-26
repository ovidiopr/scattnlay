import { describe, it, expect, beforeAll } from 'vitest';
// Point directly to the generated JS for testing
import createNmiejsModule from '../public/wasm/nmiejs.js';

describe('Scattnlay WASM Engine', () => {
  let nmie: any;

  beforeAll(async () => {
    const module = await createNmiejsModule();
    nmie = new module.nmie();
  });

  it('should initialize the nmie class', () => {
    expect(nmie).toBeDefined();
    expect(typeof nmie.RunMieCalculation).toBe('function');
  });

  it('should calculate correct Qsca for a gold-like sphere', () => {
    // Standard test case: R=100nm, n=4.0, k=0.01, WL=800nm
    nmie.ClearTarget();
    nmie.SetWavelength(800);
    nmie.AddTargetLayerReIm(100, 4.0, 0.01);
    nmie.SetModeNmaxAndType(-1, -1);
    
    nmie.RunMieCalculation();
    
    const qsca = nmie.GetQsca();
    // Verify results are numbers and within physical expectations
    expect(qsca).toBeGreaterThan(0);
    expect(qsca).toBeLessThan(10);
    expect(qsca).toBeCloseTo(5.87, 1); // Validated against scattnlay-dp
  });

  it('should handle target clearing and resets', () => {
    nmie.ClearTarget();
    nmie.AddTargetLayerReIm(50, 1.5, 0);
    nmie.RunMieCalculation();
    const q1 = nmie.GetQsca();
    
    nmie.ClearTarget();
    nmie.AddTargetLayerReIm(100, 1.5, 0);
    nmie.RunMieCalculation();
    const q2 = nmie.GetQsca();
    
    expect(q1).not.toEqual(q2);
  });
});
