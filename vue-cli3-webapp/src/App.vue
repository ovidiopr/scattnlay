<template>
  <div id="app">
    <b-loading :active.sync="isLoading" :can-cancel="isLoading"></b-loading>

    <section style="padding: 1rem">
      <button class="button is-primary is-medium"
              @click="isShowInfo = true">
        Show info
      </button>
      <b-modal :active.sync="isShowInfo">
        <div class="message" style="padding: 2rem">
          <span style="text-decoration: underline;text-emphasis: before;font-weight: bold;">Usage:</span>
          <br><br>
          Feel free to use provided software, however, use it at your own risk. We made our best effort to verify
          it is correct, however, we do not provide any warranty.
          <br><br>
          If it was usefull for your project, please, cite in your related paper the following reference:
          <br><br>
          <article style="margin-left: 1rem">
            <div style="font-style: italic; padding-bottom: 0.5rem">"Mie calculation of electromagnetic near-field
              for a multilayered sphere"</div>
            Konstantin Ladutenko, Umapada Pal, Antonio Rivera, Ovidio Peña-Rodríguez<br>
            <span style="font-weight: bold">Comp. Phys. Comm., vol. 214, pp. 225–230, 2017</span>
          </article>
        </div>
      </b-modal>
    </section>

    <section>
      <div class="field is-horizontal">
        <div class="field-label is-normal">
          <label class="label">Units</label>
        </div>
        <div class="field-body">
          <b-select v-model="units">
            <option units="nm">nm</option>
            <option units="mkm">mkm</option>
            <option units="mm">mm</option>
            <option units="cm">cm</option>
            <option units="m">m</option>
            <option units="km">km</option>
          </b-select>
        </div>
      </div>

      <div class="field is-horizontal">
        <div class="field-label is-normal">
          <label class="label">Wavelength</label>
        </div>
        <div class="field-body">
          <div class="field is-grouped is-grouped-multiline">
            <input-with-units title="from" v-bind:units="units"
                              v-bind:value="simulationSetup.fromWL"
                              @newdata="simulationSetup.fromWL=$event"></input-with-units>
            <input-with-units title="to" v-bind:units="units"
                              v-bind:value="simulationSetup.toWL"
                              @newdata="simulationSetup.toWL=$event"></input-with-units>
            <input-with-units title="step" v-bind:units="units"
                              v-bind:value="simulationSetup.stepWL"
                              @newdata="simulationSetup.stepWL=$event"></input-with-units>
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
                              @newdata="simulationSetup.R=$event"></input-with-units>
            <input-with-units title="Re(n)" units=""
                              v-bind:value="simulationSetup.reN"
                              @newdata="simulationSetup.reN=$event"></input-with-units>
            <input-with-units title="Im(n)" units=""
                              v-bind:value="simulationSetup.imN"
                              @newdata="simulationSetup.imN=$event"></input-with-units>
          </div>
        </div>
      </div>

      <div  class="field is-horizontal">
        <div class="field-label is-normal">
          <label class="label">Modes to plot</label>
        </div>
        <div class="field-body">
          <b-input v-model="simulationSetup.total_mode_n" type='number' min=1 style="width:7rem"></b-input>
        </div>
      </div>
      <div>
        <b-button class="is-primary is-medium" @click="runMie();plotResults();">
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
          <th v-for="mode_name in simulationRuntime.mode_n_names" v-bind:key="mode_name">
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
      <reactive-chart :chart="chart"></reactive-chart>
    </div>

  </div>
</template>

<script>
  import InputWithUnits from "./components/InputWithUnits.vue";
  import ReactiveChart from "./components/ReactiveChart.vue";

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

  // Test the size of wasm file
  fetch('nmiejs.wasm'
  ).then(response =>
          response.arrayBuffer()
  ).then(bytes =>
          console.log(bytes)
  );

  import nmiejs from './nmiejs.js';
  const module = nmiejs({
    locateFile(path) {
      console.log(path);
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

  const range = (start, stop, step = 1) =>
          Array(Math.ceil((stop - start) / step)).fill(start).map((x, y) => x + y * step);

  function rangeInt(size, startAt = 0) {
    return [...Array(size).keys()].map(i => i + startAt);
  }
  export default {
    name: 'app',
    components: {
      InputWithUnits,
      ReactiveChart
      // HelloWorld
    },
    data() {
        return {
          nmie: undefined,
          ttime: "",
          window: {
            width: 0,
            height: 0
          },
          // on change of initial value for __ units __
          // remember to update this.chart.layout.xaxis.title
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
          plotSelectorNames: undefined,
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
              title: 'reactive charts',
              xaxis: {
                // initial value should be same as in units
                title: 'Wavelength, nm'
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
        window.addEventListener('resize', this.handleResize)
        this.handleResize();
      },
    destroyed() {
        window.removeEventListener('resize', this.handleResize)
      },
    mounted() {
      module.onRuntimeInitialized = () => {
        this.nmie = new module.nmie();
        this.runMie();
        this.setIsPlotMode();
        let mode_n_names = this.setModeNNames();
        this.plotResults();
        this.isLoading = false;
      }
      // this.setPlotSelectorDataAndColumns(mode_n_names);
    },
    watch: {
        plotSelector: {
          handler: function () {
            this.plotResults();
          },
          deep: true
        },
        'simulationSetup.total_mode_n': function (newVal, oldVal) {
          this.setIsPlotMode();
          let mode_n_names = this.setModeNNames();
        },
        units: {
          handler: function () {
            this.chart.layout.xaxis.title = "Wavelength, " + this.units;
          }
        },
        window: {
          handler: function () {
            this.chart.layout.width = this.window.width * 0.9;
            this.chart.layout.height = this.window.height * 0.75;
            if (this.window.width < 600) this.chart.layout.width = this.window.width * 1.0;
            if (this.window.height < 600) this.chart.layout.height = this.window.height * 1.0;

          },
          deep: true
        }
      },
    methods: {
        handleResize() {
          this.window.width = window.innerWidth;
          this.window.height = window.innerHeight;
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
          // this.simulationRuntime.mode_n_names = this.setModeNNames(total_mode_n);
          // console.log(this.simulationRuntime.mode_n_names)
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
          console.log("It took " + this.ttime + " s.");

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
          console.log(this.plotSelectorNames)
          console.log("sp,e 'as;ld")
          return mode_n_names;
        },
        plotResults: function () {
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
