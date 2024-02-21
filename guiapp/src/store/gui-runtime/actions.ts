import { ActionTree } from 'vuex';
import { StateInterface } from '../index';
import { guiRuntimeStateInterface } from './state';
import { composeLabelFromPageData } from 'components/utils';
import { load } from 'js-yaml';
import Spline from 'cubic-spline-ts';

async function loadMaterialData(
  filename: string
): Promise<number[][] | undefined> {
  // TODO enable eslint, which is disabled now due to unknown result type of load() from js-yaml
  // TODO implement formulas
  // TODO implement multiple '- type:' in one file (e.g. it can be 'tabulated n' and 'tabulated k'

  try {
    const response = await fetch(
      process.env.publicPath +
        'refractiveindex.info-database/database/data/' +
        filename
    );
    const Ag_data = await response.text();

    const doc = (await load(Ag_data)) as any;
    if (
      doc.DATA[0].type == 'tabulated nk' ||
      doc.DATA[0].type == 'tabulated n'
    ) {
      const csv = doc.DATA[0].data;
      const rows: string[] = csv.split('\n');
      const data = rows.map((row) => row.split(' '));
      data.pop();
      const data_num = data.map((elem) => elem.map((elem2) => parseFloat(elem2)));

      function transpose(array: number[][]) {
        return array[0].map((col, i) => array.map((row) => row[i]));
      }

      const data_columns = transpose(data_num);
      // Convert from default refractiveindex.info mkm to nm

      for (let i = 0; i < data_columns[0].length; i++)
        data_columns[0][i] *= 1000;
      return data_columns;
    }
  } catch (e) {
    console.log(e);
  }
  return undefined;
}

const actions: ActionTree<guiRuntimeStateInterface, StateInterface> = {
  async activateMaterial({ commit, state } /* context */, filepath: string) {
    const data_columns: number[][] | undefined = await loadMaterialData(
      filepath
    );
    if (!data_columns) return;
    let xs: number[] = data_columns[0];
    let ys1: number[] = data_columns[1];
    let ys2: number[] = data_columns[1].map(() => 0);
    if (data_columns[2]) ys2 = data_columns[2];
    const maxVal = 350; // TODO move it to config.ts
    if (xs.length > maxVal) {
      const delta = Math.floor(xs.length / maxVal);
      const tmp_xs: number[] = [];
      const tmp_ys1: number[] = [];
      const tmp_ys2: number[] = [];
      for (let i = 0; i < xs.length; i = i + delta) {
        tmp_xs.push(xs[i]);
        tmp_ys1.push(ys1[i]);
        tmp_ys2.push(ys2[i]);
      }
      xs = tmp_xs;
      ys1 = tmp_ys1;
      ys2 = tmp_ys2;
    }

    // TODO use 10.1016/j.cagd.2010.10.002 or https://en.wikipedia.org/wiki/Monotone_cubic_interpolation
    const spline_n = new Spline(xs, ys1);
    const spline_k = new Spline(xs, ys2);

    const name = composeLabelFromPageData(filepath);
    const spectrumRangeStart = xs[0];
    const spectrumRangeEnd = xs[xs.length - 1];
    commit('addMaterial', {
      name: name,
      spectrumRangeStart: spectrumRangeStart,
      spectrumRangeEnd: spectrumRangeEnd,
      nSpline: spline_n,
      kSpline: spline_k,
      isPlot: false,
    });
  },
};

export default actions;
