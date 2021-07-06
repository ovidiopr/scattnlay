import Vue from 'vue'
import Vuex from 'vuex'

Vue.use(Vuex)

export default new Vuex.Store({
    state: {
        simulationSetup: {
            hostIndex: 1,
            stepWL: 2,
            fromWL: "100*(1+2)",
            fromWL_hasConflict: false,
            toWL: 1000,
            toWL_hasConflict:false,
            layers: [
                {
                    R: 100.0,
                    material: 'nk',
                    isMaterialLoaded:true,
                    isMaterial_hasConflict:false,
                    reN: "sqrt(16)",
                    imN: 0.01,
                    index: 0
                }
            ],
            total_mode_n: 4
        },
    },
    mutations: {
        increment (state) {
            state.count++
        }
    }
})
