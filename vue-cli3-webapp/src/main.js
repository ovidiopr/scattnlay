import Vue from 'vue'
import App from './App.vue'
import Buefy from 'buefy'
import 'buefy/dist/buefy.css'

// import fibonacci from './fibonacci.js';
// import fibonacciModule from './fibonacci.wasm';
// const module = fibonacci({
//   locateFile(path) {
//     if(path.endsWith('*.wasm')) {
//       return fibonacciModule;
//     }
//     return path;
//   }
// });
// module.onRuntimeInitialized = () => {
//   alert(module._fib(12));
// };

// const response = await fetch("fibonacci.wasm");
// const buffer = await response.arrayBuffer();
// const obj = await WebAssembly.instantiate(buffer);
// obj.instance.exports.exported_func();





// fetch('fibonacci.wasm').then(response =>
//     response.arrayBuffer()
// ).then(bytes =>
// {
//   console.log(bytes);
// }
// );

Vue.use(Buefy)

Vue.config.productionTip = false

new Vue({
  render: h => h(App),
}).$mount('#app')
