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
      <div class="field is-horizontal">
        <div class="field-label is-normal">
          <label class="label">
            <div v-if="isSourceOtherUnits"> Frequency  </div>
            <div v-else>               Wavelength              </div>
          </label>
        </div>
        <div class="field-body">
          <div class="field is-grouped is-grouped-multiline">
            <input-with-units title="from" v-bind:units="source_units"
                              v-bind:value="simulationSetup.fromWL"
                              @newdata="simulationSetup.fromWL=$event"/>
            <input-with-units title="to" v-bind:units="source_units"
                              v-bind:value="simulationSetup.toWL"
                              @newdata="simulationSetup.toWL=$event"/>
            <input-with-units title="step" v-bind:units="source_units"
                              v-bind:value="simulationSetup.stepWL"
                              @newdata="simulationSetup.stepWL=$event"/>
          </div>
        </div>
      </div>


      <div class="field is-horizontal">
        <div class="field-label is-normal">
          <label class="label">Spherical particle</label>
        </div>
        <div class="field-body">
          <div class="field is-grouped is-grouped-multiline">
            <input-with-units title="R" v-bind:units="units"
                              v-bind:value="simulationSetup.R"
                              @newdata="simulationSetup.R=$event"/>
            <input-with-units title="Re(n)" units=""
                              v-bind:value="simulationSetup.reN"
                              @newdata="simulationSetup.reN=$event"/>
            <input-with-units title="Im(n)" units=""
                              v-bind:value="simulationSetup.imN"
                              @newdata="simulationSetup.imN=$event"/>
          </div>
        </div>
      </div>

      <div  class="field is-horizontal">
        <div class="field-label is-normal">
          <label class="label">Modes to plot</label>
        </div>
        <div class="field-body">
          <b-input v-model="simulationSetup.total_mode_n" type='number' min=1 style="width:7rem"/>
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

  import InputWithUnits from "./components/InputWithUnits.vue";
  import ReactiveChart from "./components/ReactiveChart.vue";
  import ShowInfo from "./components/ShowInfo.vue";
  import GetUnits from "./components/GetUnits.vue";
  export default {
    name: 'app',
    components: {
      GetUnits,
      InputWithUnits,
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
          source_units: 'nm',
          isSourceOtherUnits: false,
          units: 'nm',
          simulationRuntime: {
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
          simulationSetup: {
            stepWL: 0.5,
            fromWL: 300.0,
            toWL: 1000.0,
            R: 100.0,
            reN: 4.0,
            imN: 0.01,
            total_mode_n: 4
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
        }
      },
      source_units: {
        handler: function () {
          this.setXaxisTitle();
        }
      },
      isSourceOtherUnits: {
        handler: function () {
          this.setEmptyChart();
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
          let t0 = performance.now();
          let fromWL = parseFloat(this.simulationSetup.fromWL);
          let toWL = parseFloat(this.simulationSetup.toWL);
          let stepWL = parseFloat(this.simulationSetup.stepWL);
          let R = parseFloat(this.simulationSetup.R);
          let reN = parseFloat(this.simulationSetup.reN);
          let imN = parseFloat(this.simulationSetup.imN);

          let Qsca = [], Qabs = [];
          let Qsca_n = [[], []], Qabs_n = [[], []];
          const nmie = this.nmie;
          nmie.ClearTarget();
          nmie.AddTargetLayerReIm(R, reN, imN);
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
          WLs.forEach(function (WL) {
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
          });
          this.simulationRuntime.Qsca = Qsca;
          this.simulationRuntime.Qabs = Qabs;
          this.simulationRuntime.Qsca_n = Qsca_n;
          this.simulationRuntime.Qabs_n = Qabs_n;

          let t1 = performance.now();
          this.ttime = ((t1 - t0) / 1000).toFixed(2);
          // console.log("It took " + this.ttime + " s.");

          this.changes++;

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
        if (this.isSourceOtherUnits) {
          this.chart.layout.xaxis.title = "Frequency, " + this.source_units;
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
          x: this.simulationRuntime.WLs,
          y: this.simulationRuntime.Qsca,
          type: 'scatter',
          name: 'Qsca'
        };
        traceQabs = {
          x: this.simulationRuntime.WLs,
          y: this.simulationRuntime.Qabs,
          type: 'scatter',
          name: 'Qabs'
        };
        this.chart.traces = [];
        if (this.plotSelector.isPlotQsca === true) this.chart.traces.push(traceQsca);
        if (this.plotSelector.isPlotQabs === true) this.chart.traces.push(traceQabs);
      },
      setEmptyChart: function () {
        this.chart.traces = [
          {
            y: [],
            line: {
              color: "#5e9e7e",
              width: 4,
              shape: "line"
            }
          }];
        this.setXaxisTitle();
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
                x: this.simulationRuntime.WLs,
                y: this.simulationRuntime.Qsca_n[mode_type][mode_n],
                  type: 'scatter',
                name: 'Qsca ' + mode_names[mode_type] + ' ' + mode_n_names[mode_n + 1].name
              };
              this.chart.traces.push(trace_sca);
            }
            if (this.plotSelector.isPlotQabs === true) {
              let trace_abs = {
                x: this.simulationRuntime.WLs,
                y: this.simulationRuntime.Qabs_n[mode_type][mode_n],
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
