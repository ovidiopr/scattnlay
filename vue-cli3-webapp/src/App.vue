<template>
  <div id="app">
    <b-loading :active.sync="isLoading" :can-cancel="isLoading"/>
    <ShowInfo/>
    <section>
      <GetUnits v-bind:units="units"
                v-bind:source_units="source_units"
                v-bind:isSourceOtherUnits="isSourceOtherUnits"
                @unitsData="units=$event"
                @source_unitsData="source_units=$event"
                @isSourceOtherUnitsData="isSourceOtherUnits=$event"
      />
      <GetHostIndex v-bind:hostIndex="simulationSetup.hostIndex"
                    @hostIndexData="simulationSetup.hostIndex=$event" />

      <GetSourceParameters v-bind:fromWL="simulationSetup.fromWL"
                           v-bind:toWL="simulationSetup.toWL"
                           v-bind:stepWL="simulationSetup.stepWL"
                           v-bind:source_units="source_units"
                           @fromWLData="simulationSetup.fromWL=$event"
                           @toWLData="simulationSetup.toWL=$event"
                           @stepWLData="simulationSetup.stepWL=$event"
                           @source_unitsData="source_units=$event"
      />
      <GetMaterials v-bind:materials="materials"/>
      <GetParticleParameters v-bind:layers="simulationSetup.layers"
                             v-bind:units="units"
                             v-bind:materials="materials"
      />

      <div  class="field is-horizontal">
        <div class="field-label is-normal">
          <label class="label">Modes to plot</label>
        </div>
        <div class="field-body" >
          <b-field label="">
            <b-input v-model="simulationSetup.total_mode_n" type='number' min=1 style="width:7rem"/>
          </b-field>
        </div>
      </div>
      <div>
        <b-button class="is-primary is-medium" @click="runSimulation();">
          Run simulation
        </b-button><br>
        It took {{ ttime }} s.
      </div>
      <div class="field is-grouped is-grouped-multiline">
        <b-switch v-model="plotSelector.isPlotQsca">
          Qsca
        </b-switch>
        <b-switch v-model="plotSelector.isPlotQabs">
          Qabs
        </b-switch>
      </div>
      <!--        <b-table :data="plotSelectorData" :columns="plotSelectorColumns" :mobile-cards="isShowInfo"></b-table>-->
      <table class="table is-narrow">
        <thead>
        <tr>
          <th v-for="mode_name in simulationRuntime.mode_n_names" v-bind:key="mode_name.name">
            {{ mode_name.name == 'type' ? '': mode_name.name}}
          </th>
        </tr>
        </thead>
        <tbody>
        <tr>
          <th>E</th>
          <td v-for="(mode,index) in plotSelector.isPlotModeE" v-bind:key="index">
            <template v-if="(simulationRuntime.total_mode_n_evaluated - 1) < index">
              <b-checkbox v-model="plotSelector.isPlotModeE[index]" disabled> </b-checkbox>
            </template>
            <template v-else>
              <b-checkbox v-model="plotSelector.isPlotModeE[index]"> </b-checkbox>
            </template>
          </td>
        </tr>
        <tr>
          <th>H</th>
          <td v-for="(mode,index) in plotSelector.isPlotModeH" v-bind:key="index">
            <template v-if="(simulationRuntime.total_mode_n_evaluated - 1) < index">
              <b-checkbox v-model="plotSelector.isPlotModeH[index]" disabled> </b-checkbox>
            </template>
            <template v-else>
              <b-checkbox v-model="plotSelector.isPlotModeH[index]"> </b-checkbox>
            </template>
          </td>
        </tr>
        </tbody>
      </table>

    </section>
    <div class="chart-container">
      <reactive-chart :chart="chart"/>
    </div>

  </div>
</template>

<script>
  // You should put *.wasm to public/ and *.js to src/ folder

  // To compile fibbonacci example use
  //   emcc -O3 -s WASM=1 -s EXTRA_EXPORTED_RUNTIME_METHODS='["cwrap"]' -s ALLOW_MEMORY_GROWTH=1 -s MODULARIZE=1 -s 'EXPORT_NAME="fibonacci"' -o ./fibonacci.js fibonacci.c
  // for and example from https://gist.github.com/ashleygwilliams/32c31a3f5b8c87bf2894108b3534ee4f
  // import fibonacci from './fibonacci.js';
  // const module = fibonacci({
  //   locateFile(path) {
  //     console.log(path);
  //     return path;
  //   }
  // });
  // module.onRuntimeInitialized = () => {
  //   console.log(module._fib(12));
  // };
  // // Test the size of wasm file
  // fetch('nmiejs.wasm'
  // ).then(response =>
  //         response.arrayBuffer()
  // ).then(bytes =>
  //         console.log(bytes)
  // );

  import nmiejs from './nmiejs.js';
  const module = nmiejs({
    locateFile(path) {
      // console.log(path);
      return path;
    }
  });

  // // Test nmiejs if working
  // module.onRuntimeInitialized = () => {
  //   const nmie = new module.nmie();
  //   nmie.ClearTarget();
  //   let R = 100.0;
  //   let reN = 4.0;
  //   let imN = 0.01;
  //   nmie.AddTargetLayerReIm(R, reN, imN)
  //   nmie.SetModeNmaxAndType(-1, -1);
  //   let WL = 800;
  //   nmie.SetWavelength(WL);
  //   nmie.RunMieCalculation();
  //   console.log(nmie.GetQsca());
  // }

  const range = (start, stop, step = 1) => Array(Math.ceil((stop - start) / step)).fill(start).map((x, y) => x + y * step);

  function rangeInt(size, startAt = 0) {
    return [...Array(size).keys()].map(i => i + startAt);
  }

  import ReactiveChart from "./components/ReactiveChart.vue";
  import ShowInfo from "./components/ShowInfo.vue";
  import GetMaterials from "./components/GetMaterials.vue";
  import GetHostIndex from "./components/GetHostIndex.vue";
  import GetUnits from "./components/GetUnits.vue";
  import GetSourceParameters from "./components/GetSourceParameters.vue";
  import GetParticleParameters from "./components/GetParticleParameters.vue";

  export default {
    name: 'app',
    components: {
      GetHostIndex,
      GetMaterials,
      GetSourceParameters,
      GetParticleParameters,
      GetUnits,
      ReactiveChart,
      ShowInfo
    },
    data() {
        return {
          nmie: undefined,
          ttime: "",
          window: {
            width: 0,
            height: 0
          },
          units: 'nm',
          units_prev: 'nm',
          source_units: 'nm',
          source_units_prev: 'nm',
          materials: [],
          isSourceOtherUnits: false,
          simulationSetup: {
            hostIndex: 1,
            stepWL: 0.5,
            fromWL: 300.0,
            toWL: 1000.0,
            layers: [
              {
                R: 100.0,
                material: 'nk',
                reN: 4.0,
                imN: 0.01,
                index: 0
              }
            ],
            total_mode_n: 4
          },
          simulationRuntime: {
            r_units: 'nm',
            r_source_units: 'nm',
            stepWL: 0.5,
            fromWL: 300.0,
            toWL: 1000.0,
            layers: [
              {
                R: 100.0,
                material: 'nk',
                reN: 4.0,
                imN: 0.01,
                index:0
              }
            ],
            mode_n: [],
            mode_n_names: [],
            mode_types: range(0, 2),
            mode_names: ['E', 'H'],
            WLs: [],
            Qsca: [],
            Qabs: [],
            Qsca_n: [[], []],
            Qabs_n: [[], []],
            layout: {},
            trace1: {},
            trace2: {},
            total_mode_n_evaluated: 4,
          },
          plotData: {
            pWLs: [],
            pQsca: [],
            pQabs: [],
            pQsca_n: [[], []],
            pQabs_n: [[], []],
          },
          changes: 0,
          isShowInfo: false,
          isLoading: true,
          isRunning: false,
          plotSelector: {
            isPlotModeE: [],
            isPlotModeH: [],
            isPlotQabs: true,
            isPlotQsca: true
          },
          // plotSelectorData: undefined,
          // plotSelectorColumns: undefined,
          // plotSelectorNames: undefined,
          chart: {
            uuid: "123",
            traces: [
              {
                y: [],
                line: {
                  color: "#5e9e7e",
                  width: 4,
                  shape: "line"
                }
              }
            ],
            layout: {
              // title: 'reactive charts',
              xaxis: {
                // will be set on mount
                title: ''
              },
              yaxis: {
                title: 'Normalized cross-sections'
              },
              width: 100,
              height: 100
            }
          }
        };
      },
    created() {
        window.addEventListener('resize', this.handleResize);
        this.handleResize();
      },
    destroyed() {
        window.removeEventListener('resize', this.handleResize)
      },
    mounted() {
      this.setXaxisTitle();
      module.onRuntimeInitialized = () => {
        this.nmie = new module.nmie();
        this.runMie();
        this.setIsPlotMode();
        this.setModeNNames();
        this.plotResults();
        this.isLoading = false;
      }


    },
    watch: {
      plotSelector: {
          handler: function () {
            this.plotResults();
          },
          deep: true
        },
      'simulationSetup.total_mode_n': function () {
        this.setIsPlotMode();
        this.setModeNNames();
      },
      units: {
        handler: function () {
          if (!this.isSourceOtherUnits) {
            this.source_units = this.units;
          }
          let u = this.units;
          let u_prev = this.units_prev;
          let layers = this.simulationSetup.layers;
          for (let i=0; i<layers.length; i++) {
            layers[i].R = this.convertUnits(u_prev, u, layers[i].R);
          }
          this.units_prev = this.units;
        }
      },
      source_units: {
        handler: function () {
          // this.$buefy.notification.open({
          //   duration: 4000,
          //   message: 'Source parameters were restored from last simulation for correct plotting.',
          //   type: 'is-info',
          //   position: 'is-top',
          // });
          this.setXaxisTitle();

          let u = this.source_units;
          let u_prev = this.source_units_prev;
          for (let i=0; i<this.plotData.pWLs.length; i++) {
            this.plotData.pWLs[i] = this.convertUnits(u_prev, u,
                    this.plotData.pWLs[i]);
          }
          this.simulationSetup.fromWL = this.convertUnits(u_prev, u,
                  this.simulationSetup.fromWL);
          this.simulationSetup.toWL = this.convertUnits(u_prev, u,
                  this.simulationSetup.toWL);
          if (this.simulationSetup.fromWL > this.simulationSetup.toWL) {
            let tmp = this.simulationSetup.fromWL;
            this.simulationSetup.fromWL = this.simulationSetup.toWL;
            this.simulationSetup.toWL = tmp;
          }
          this.simulationSetup.stepWL = (this.simulationSetup.toWL-this.simulationSetup.fromWL)
                  /this.plotData.pWLs.length;
          this.source_units_prev = this.source_units;

          // Update plotting
          this.plotSelector.isPlotQsca = !this.plotSelector.isPlotQsca
          setTimeout(
                  () => {
                    this.plotSelector.isPlotQsca = !this.plotSelector.isPlotQsca
                  }, 2);

        }
      },
      isSourceOtherUnits: {
        handler: function () {
          // this.setEmptyChart();
          this.setXaxisTitle();
          if (!this.isSourceOtherUnits) {
            this.source_units = this.units;
          } else {
            this.source_units = 'THz'
          }
        }
      },
      window: {
        handler: function () {
            this.chart.layout.width = this.window.width * 0.9;
            this.chart.layout.height = this.window.height * 0.75;
            if (this.window.width < 600) this.chart.layout.width = this.window.width;
            if (this.window.height < 600) this.chart.layout.height = this.window.height;

          },
        deep: true
      }
    },
    methods: {
      handleResize() {
          this.window.width = window.innerWidth;
          this.window.height = window.innerHeight*0.8;
        },
      convertUnits2nm(fromU, val) {
        if (fromU === 'nm') return val;
        if (fromU === 'mkm') return val*1e3;
        if (fromU === 'mm') return val*1e6;
        if (fromU === 'cm') return val*1e7;
        if (fromU === 'm') return val*1e9;
        if (fromU === 'km') return val*1e12;

        let c = 299792458; // m/s
        let hc = 1239841930.092394328e-15; // m*eV
        if (fromU === 'THz') return c/(val*1e12)*1e9;
        if (fromU === 'GHz') return c/(val*1e9)*1e9;
        if (fromU === 'MHz') return c/(val*1e6)*1e9;
        if (fromU === 'kHz') return c/(val*1e3)*1e9;
        if (fromU === 'Hz') return c/(val*1e0)*1e9;

        if (fromU === 'eV') return hc/(val*1e0)*1e9;
        if (fromU === 'meV') return hc/(val/1e3)*1e9;

        if (fromU === 'fs') return (val/1e12)/c*1e9;
        if (fromU === 'ps') return (val/1e15)/c*1e9;

        return undefined;
      },
      convertUnitsFrom_nm(toU, val) {
        if (toU === 'nm') return val;
        if (toU === 'mkm') return val/1e3;
        if (toU === 'mm') return val/1e6;
        if (toU === 'cm') return val/1e7;
        if (toU === 'm') return val/1e9;
        if (toU === 'km') return val/1e12;

        let c = 299792458; // m/s
        let hc = 1239841930.092394328e-15; // m*eV
        if (toU === 'THz') return c/(val/1e9)/1e12;
        if (toU === 'GHz') return c/(val/1e9)/1e9;
        if (toU === 'MHz') return c/(val/1e9)/1e6;
        if (toU === 'kHz') return c/(val/1e9)/1e3;
        if (toU === 'Hz') return c/(val/1e9)/1e0;

        if (toU === 'eV') return hc/(val/1e9);
        if (toU === 'meV') return hc/(val/1e9)*1e3;

        if (toU === 'fs') return (val/1e9)/c*1e12;
        if (toU === 'ps') return (val/1e9)/c*1e15;

        return undefined;
      },
      convertUnits(fromU,toU, val) {
        if (fromU === toU) return val;
        return this.convertUnitsFrom_nm(toU, this.convertUnits2nm(fromU, val));
      },
      runSimulation: function() {
          this.$buefy.notification.open({
            duration: 200,
            message: 'Simulation was started!',
            type: 'is-danger',
            position: 'is-top',
            });
          setTimeout(
                  () => {
                    this.runMie();
                    this.plotResults();
                    this.$buefy.notification.open({
                      duration: 3000,
                      message: 'Finished! '+"It took " + this.ttime + " s.",
                      type: 'is-success',
                      position: 'is-top',
                      })
                    ;
                  }, 20);
        },
      runMie: function () {
        this.simulationRuntime.r_units = this.units;
        this.simulationRuntime.r_source_units = this.source_units;
        this.simulationRuntime.fromWL = this.simulationSetup.fromWL;
        this.simulationRuntime.toWL = this.simulationSetup.toWL;
        this.simulationRuntime.stepWL = this.simulationSetup.stepWL;
        this.simulationRuntime.layers = JSON.parse(JSON.stringify(this.simulationSetup.layers));

        let t0 = performance.now();
        let fromWL = parseFloat(this.simulationSetup.fromWL);
        let toWL = parseFloat(this.simulationSetup.toWL);
        let stepWL = parseFloat(this.simulationSetup.stepWL);
        let host = parseFloat(this.simulationSetup.hostIndex);


        let Qsca = [], Qabs = [];
        let Qsca_n = [[], []], Qabs_n = [[], []];

        let WLs = range(fromWL, toWL, stepWL);
        this.simulationRuntime.WLs = WLs;
        let total_mode_n = this.simulationSetup.total_mode_n;
        let mode_n = [];
        mode_n = rangeInt(Number(total_mode_n), 1);
        this.simulationRuntime.total_mode_n_evaluated = total_mode_n;
        let mode_types = range(0, 2);
        this.simulationSetup.mode_n = mode_n;
        this.simulationSetup.mode_types = mode_types;
        mode_types.forEach(function (mode_type) {
          mode_n.forEach(function () {
            Qsca_n[mode_type].push([]);
            Qabs_n[mode_type].push([]);
            });
        });

        // Provide all sizes (germetry and wavelengths) to nmie in same units [nm]
        const nmie = this.nmie;
        // nmie.ClearTarget();
        // nmie.AddTargetLayerReIm( this.convertUnits2nm(this.units, R)*host,
        //         reN/host, imN/host);
        let WLs_nm = range(fromWL, toWL, stepWL);
        let WL_points = WLs_nm.length;
        for (let i = 0; i < WL_points; i++) {
          let WL = this.convertUnits2nm(this.source_units, WLs[i]);
          nmie.ClearTarget();
          for (let num_layer = 0;
               num_layer < this.simulationRuntime.layers.length; 
               num_layer++) {
            let R = parseFloat(this.simulationRuntime.layers[num_layer].R);
            let reN = parseFloat(this.simulationRuntime.layers[num_layer].reN);
            let imN = parseFloat(this.simulationRuntime.layers[num_layer].imN);
            nmie.AddTargetLayerReIm( this.convertUnits2nm(this.units, R)*host,
                    reN/host, imN/host);
            
          }
          nmie.SetModeNmaxAndType(-1, -1);
          nmie.SetWavelength(WL);
          nmie.RunMieCalculation();
          Qsca.push(nmie.GetQsca());
          Qabs.push(nmie.GetQabs());
          mode_types.forEach(function (mode_type) {
            mode_n.forEach(function (n) {
              nmie.SetModeNmaxAndType(n, mode_type);
              nmie.RunMieCalculation();
              Qsca_n[mode_type][n - 1].push(nmie.GetQsca());
              Qabs_n[mode_type][n - 1].push(nmie.GetQabs());
            });
          });
        }
        this.simulationRuntime.Qsca = Qsca;
        this.simulationRuntime.Qabs = Qabs;
        this.simulationRuntime.Qsca_n = Qsca_n;
        this.simulationRuntime.Qabs_n = Qabs_n;
        this.filterBug();  // TODO: fix the algorithm instead of filtering the final output

        let t1 = performance.now();
        this.ttime = ((t1 - t0) / 1000).toFixed(2);
        // console.log("It took " + this.ttime + " s.");
        
        this.plotData.pWLs = WLs.slice();
        this.plotData.pQsca = Qsca;
        this.plotData.pQabs = Qabs;
        this.plotData.pQsca_n = Qsca_n;
        this.plotData.pQabs_n = Qabs_n;

        this.changes++;

      },
      filterBug: function () {
        let Qsca = this.simulationRuntime.Qsca;
        let Qabs = this.simulationRuntime.Qabs;
        let Qsca_n = this.simulationRuntime.Qsca_n;
        let Qabs_n = this.simulationRuntime.Qabs_n;
        let total_mode_n = this.simulationSetup.total_mode_n;
        for (let mode_type = 0; mode_type < 2; mode_type++) {
          for (let n = 0; n < total_mode_n; n++){
            this.filterMedian(Qsca_n[mode_type][n]);
            this.filterMedian(Qabs_n[mode_type][n]);
          }
        }
        this.filterMedian(Qsca);
        this.filterMedian(Qabs);
        this.simulationRuntime.Qsca = Qsca;
        this.simulationRuntime.Qabs = Qabs;
        this.simulationRuntime.Qsca_n = Qsca_n;
        this.simulationRuntime.Qabs_n = Qabs_n;

      },
      filterMedian: function(spectra) {
        let l = spectra.length;
        let i;
        for (i = 1; i < l-1; i++) {
          let prev = spectra[i-1];
          let curr = spectra[i];
          let next = spectra[i+1];
          let diff1 = Math.abs((prev-curr)/curr);
          let diff2 = Math.abs((next-curr)/curr);
          if (diff1 > 0.3 && diff2 > 0.3) {
            spectra[i] = (prev+next)/2.0;
          }
           // console.log(spectra(i));
        }
        // return spectra;
      },
      setIsPlotMode: function () {
          let total_mode_n = this.simulationSetup.total_mode_n;
          total_mode_n++;
          let np1 = Number(total_mode_n);
          let mode_n = rangeInt(Number(np1), 0);
          let modeE = [];
          let modeH = [];
          mode_n.forEach(function (n) {
            if (n !== 0) {
              modeE.push(false);
              modeH.push(false);
            }
          });

          let saved_modes = this.simulationSetup.total_mode_n > this.simulationRuntime.total_mode_n_evaluated ?
                  this.simulationRuntime.total_mode_n_evaluated : this.simulationSetup.total_mode_n;
          this.simulationRuntime.total_mode_n_evaluated = saved_modes;
          let oldModeE = this.plotSelector.isPlotModeE;
          let oldModeH = this.plotSelector.isPlotModeH;
          let ii;
          for (ii = 0; ii < saved_modes; ++ii) {
            if (ii < oldModeE.length) {
              modeE[ii] = oldModeE[ii];
              modeH[ii] = oldModeH[ii];
            }
          }
          this.plotSelector.isPlotModeE = modeE;
          this.plotSelector.isPlotModeH = modeH;
        },
      setXaxisTitle: function () {
        if (this.source_units.endsWith('Hz')) {
          this.chart.layout.xaxis.title = "Frequency, " + this.source_units;
        } else if (this.source_units.endsWith('eV')) {
          this.chart.layout.xaxis.title = "Energy, " + this.source_units;
        } else if (this.source_units.endsWith('s')) {
          this.chart.layout.xaxis.title = "Period, " + this.source_units;
        } else {
            this.chart.layout.xaxis.title = "Wavelength, " + this.source_units;
        }
      },
      setModeNNames: function () {
        let total_mode_n = this.simulationSetup.total_mode_n;
        total_mode_n++;
        let np1 = Number(total_mode_n);
        let mode_n = rangeInt(Number(np1), 0);
        let mode_n_names = [];
        mode_n.forEach(function (n) {
          if (n === 0) mode_n_names.push({id: 0, name: 'type'});
          else if (n === 1) mode_n_names.push({id: 1, name: 'dipole'});
          else if (n === 2) mode_n_names.push({id: 2, name: 'quadrupole'});
          else if (n === 3) mode_n_names.push({id: 3, name: 'octupole'});
          else mode_n_names.push({id: n, name: (Math.pow(2, n)).toString()});
        });
        this.simulationRuntime.mode_n_names = mode_n_names;
      },
      setQtotalChart: function () {
        let traceQsca, traceQabs;
        traceQsca = {
          x: this.plotData.pWLs,
          y: this.plotData.pQsca,
          type: 'scatter',
          name: 'Qsca'
        };
        traceQabs = {
          x: this.plotData.pWLs,
          y: this.plotData.pQabs,
          type: 'scatter',
          name: 'Qabs'
        };
        this.chart.traces = [];
        if (this.plotSelector.isPlotQsca === true) this.chart.traces.push(traceQsca);
        if (this.plotSelector.isPlotQabs === true) this.chart.traces.push(traceQabs);
      },
      plotResults: function () {
        this.setQtotalChart();
        let mode_type, mode_n;
        let mode_n_names = this.simulationRuntime.mode_n_names;
        let mode_names = this.simulationRuntime.mode_names;

        for (mode_type = 0; mode_type < 2; ++mode_type) {
          for (mode_n = 0; mode_n < this.simulationRuntime.total_mode_n_evaluated; ++mode_n) {
            let is_mode_plot = mode_type === 0 ? this.plotSelector.isPlotModeE : this.plotSelector.isPlotModeH;
            if (is_mode_plot[mode_n] === false) continue;
            if (this.plotSelector.isPlotQsca === true) {
              let trace_sca = {
                x: this.plotData.pWLs,
                y: this.plotData.pQsca_n[mode_type][mode_n],
                  type: 'scatter',
                name: 'Qsca ' + mode_names[mode_type] + ' ' + mode_n_names[mode_n + 1].name
              };
              this.chart.traces.push(trace_sca);
            }
            if (this.plotSelector.isPlotQabs === true) {
              let trace_abs = {
                x: this.plotData.pWLs,
                y: this.plotData.pQabs_n[mode_type][mode_n],
                type: 'scatter',
                name: 'Qabs ' + mode_names[mode_type] + ' ' + mode_n_names[mode_n + 1].name
                };
              this.chart.traces.push(trace_abs);
            }
          }
        }

      }
    },
  };

</script>

<style>
#app {
  font-family: 'Avenir', Helvetica, Arial, sans-serif;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  /*text-align: center;*/
  /*color: #2c3e50;*/
  /*margin-top: 60px;*/
  margin: 0;
}
</style>
