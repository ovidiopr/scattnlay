import Vue from 'vue'
import App from './App.vue'
import Buefy from 'buefy'
import 'buefy/dist/buefy.css'
// import Plotly from 'plotly.js'
Vue.use(Buefy)

import { library } from '@fortawesome/fontawesome-svg-core'
import { FontAwesomeIcon } from '@fortawesome/vue-fontawesome'
import { faCoffee, faCheckCircle, faBan, faTrash } from '@fortawesome/free-solid-svg-icons'
library.add(faCoffee, faCheckCircle, faBan, faTrash)

Vue.component('font-awesome-icon', FontAwesomeIcon)


Vue.config.productionTip = false

new Vue({
  render: h => h(App),
}).$mount('#app')
