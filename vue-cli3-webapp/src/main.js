import Vue from 'vue'
import store from './store'
import App from './App.vue'
import Buefy from 'buefy'
// import 'buefy/dist/buefy.css'
// import Plotly from 'plotly.js'
Vue.use(Buefy)

import { library } from '@fortawesome/fontawesome-svg-core'
import { FontAwesomeIcon } from '@fortawesome/vue-fontawesome'
import { faCoffee, faCheckCircle, faPlusCircle, faBan, faTrash } from '@fortawesome/free-solid-svg-icons'
library.add(faCoffee, faCheckCircle, faPlusCircle, faBan, faTrash);

Vue.component('font-awesome-icon', FontAwesomeIcon);

new Vue({
  render: h => h(App),
  store
}).$mount('#app');
