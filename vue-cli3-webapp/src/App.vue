<template>
  <div id="app">
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
      <div>From {{simulationSetup.fromWL}} to {{simulationSetup.toWL}} {{ units }}
        <br>
        Changes:{{ changes }}
      </div>

<!--      <div class="chart-container">-->
<!--        <reactive-chart :chart="chart"></reactive-chart>-->
<!--      </div>-->
    </section>

  </div>
</template>

<script>
  import InputWithUnits from "./components/InputWithUnits.vue";

  // You should put *.wasm to public and *.js to src folder

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
  // fetch('fibonacci.wasm'
  // ).then(response =>
  //         response.arrayBuffer()
  // ).then(bytes =>
  //         console.log(bytes)
  // );

  // // TODO: uncomment
  import nmiejs from './nmiejs.js';
  const module = nmiejs({
    locateFile(path) {
      console.log(path);
      return path;
    }
  });
  // import nmiejs from './nmiejs.js';
  // import nmiejsModule from './nmiejs.wasm';
  // const module = nmiejs({
  //   locateFile(path) {
  //     if(path.endsWith('.wasm')) {
  //       return nmiejsModule;
  //     }
  //     return path;
  //   }
  // });
  //
  // nmiejs().then(function(module) {
  //   // this is reached when everything is ready, and you can call methods on Module
  // });
  // let getMethods = (obj) => Object.getOwnPropertyNames(obj);
  // // TODO: end uncomment block



  //
  module.onRuntimeInitialized = () => {
    const nmie = new module.nmie();
    nmie.ClearTarget();
    let R = 100.0;
    let reN = 4.0;
    let imN = 0.01;
    nmie.AddTargetLayerReIm(R, reN, imN)
    nmie.SetModeNmaxAndType(-1, -1);
    let WL = 800;
    nmie.SetWavelength(WL);
    nmie.RunMieCalculation();
    // Qsca.push();
    // Qabs.push(nmie.GetQabs());

    console.log(nmie.GetQsca());
  }
  // console.log(Object.keys(nmie));
  //import ReactiveChart from "./components/ReactiveChart.vue";
// import VueWasm from 'vue-wasm';
// import nmie from './assets/nmie.wasm';
// VueWasm(Vue, { modules: { arithmetic: arithmeticModule } });

const range = (start, stop, step = 1) =>
        Array(Math.ceil((stop - start) / step)).fill(start).map((x, y) => x + y * step);

export default {
  name: 'app',
  components: {
    InputWithUnits,
    // ReactiveChart
    // HelloWorld
  },
  data() {
    return {
      window: {
        width: 0,
        height: 0
      },
      // on change of initial value for __ units __
      // remember to update this.chart.layout.xaxis.title
      units: 'nm',
      simulationRuntime: {
        mode_n : [],
        mode_types: range(0,2),
        mode_names: ['E','H'],
        WLs: [],
        Qsca: [],
        Qabs: [],
        Qsca_n: [[],[]],
        Qabs_n: [[],[]],
        layout: {},
        trace1: {},
        trace2: {}
      },
      simulationSetup: {
        stepWL: 0.5,
        fromWL: 300.0,
        toWL: 1000.0,
        R: 100.0,
        reN: 4.0,
        imN: 0.01,
        total_mode_n:4
      },
      changes: 0,
      isShowInfo: false,
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
          title:'reactive charts',
          xaxis: {
            // initial value should be same as in units
            title: 'Wavelength, nm'
          },
          yaxis: {
            title: 'yaxis title'
          }
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
  mounted(){ this.runMie();},
  watch: {
    simulationSetup: {
      handler: function () {
        this.runMie();
      },
      deep: true
    },
    units:{
      handler: function () {
        this.chart.layout.xaxis.title = "Wavelength, "+this.units;
      }
    }
  },
  methods: {
    handleResize() {
      this.window.width = window.innerWidth;
      this.window.height = window.innerHeight;
    },
    runMie(){

      // let t0 = performance.now();
      // let fromWL = parseFloat(this.simulationSetup.fromWL);
      // let toWL = parseFloat(this.simulationSetup.toWL);
      // let stepWL = parseFloat(this.simulationSetup.stepWL);
      // let R = parseFloat(this.simulationSetup.R);
      // let reN = parseFloat(this.simulationSetup.reN);
      // let imN = parseFloat(this.simulationSetup.imN);
      // this.simulationRuntime.Qsca = [], this.simulationRuntime.Qabs = [];
      // this.simulationRuntime.Qsca_n = [[],[]], this.simulationRuntime.Qabs_n = [[],[]];
      // nmie.ClearTarget();
      // nmie.AddTargetLayerReIm(R, reN, imN)


      this.changes++;
    }
  },

}
</script>

<style>
#app {
  font-family: 'Avenir', Helvetica, Arial, sans-serif;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
  /*text-align: center;*/
  /*color: #2c3e50;*/
  /*margin-top: 60px;*/
}
</style>
