<template>
    <div class="field">
        <!--            Select particle type-->
        <div class="field is-horizontal">
            <div class="field-label is-normal">
                <label class="label">Spherical particle</label>
            </div>
            <div class="field-body">
                <div class="field is-grouped is-grouped-multiline">
                    <b-radio v-model="particle" native-value="bulk"> bulk </b-radio>
                    <b-radio v-model="particle" native-value="core-shell"> core-shell </b-radio>
                    <b-radio v-model="particle" native-value="multilayer"> multilayer </b-radio>
                    <div v-if="particle==='multilayer'">
                        <b-input v-model="layersNum" type='number' min=1 max=10 class="binput"/>
                    </div>
                    <div v-else>
                        <b-input v-model="layersNum" type='number' min=1 max=10 class="binput" disabled/>
                    </div>
                </div>
            </div>
        </div>
        <div v-for="layer in layers" v-bind:key="layer.index">
            <GetLayerParameters v-bind:layer="layer"
                            v-bind:units="units"
                            v-bind:particle="particle"
                            v-bind:index="layer.index"
                            v-bind:materials="materials"
        />
        </div>
    </div>
</template>

<script>
    import GetLayerParameters from "./GetLayerParameters.vue";
    import { mapState } from 'vuex'
    export default {
        name: "GetParticleParameters",
        components: {
            GetLayerParameters
        },
        data () {
            return {
                // TODO: Is it OK to modify Local later?
                particle: 'bulk',
                layersNum: 1,
                index: 1
            }
        },
      computed:
          mapState([
            'simulationSetup'
          ]),
      watch: {
            // emit updated values
            particle: {
                handler: function () {
                    if (this.particle === 'bulk') {this.layersNum = 1;}
                    if (this.particle === 'core-shell') {this.layersNum = 2;}
                    if (this.particle === 'multilayer') {this.layersNum = 3;}
                }
            },
            // layers: {
            //     handler: function () {
            //         if ( this.layers[0].spline_n === undefined ) return;
            //         console.log('GPP1', Object.getOwnPropertyNames(this.layers[0].spline_n));
            //         console.log('GPP2', this.layers[0].spline_n.at(500));
            //     },
            //     deep:true
            // },
            layersNum: {
                handler: function () {
                    while (this.layersNum < this.simulationSetup.layers.length) {
                        this.simulationSetup.layers.pop();
                    }
                    let r_prev = this.simulationSetup.layers[0].R;
                    while (this.layersNum > this.simulationSetup.layers.length) {
                        // r_prev = r_prev*1.1;
                        this.simulationSetup.layers.push(
                            {
                                R: r_prev*0.1,
                                material: 'nk',
                                isMaterialLoaded:true,
                              isMaterial_hasConflict:false,

                              reN: 4.0,
                                imN: 0.01,
                                index: 1,
                                spline_n: undefined,
                                spline_k: undefined,
                            }
                        );
                    }
                    for (let i = 0; i < this.simulationSetup.layers.length; i++) {
                        this.simulationSetup.layers[i].index = i;
                    }
                }
            },
            deep: true
        },
        props: ['layers', 'units','materials']
    }
</script>

<style scoped>
    .binput {
        display:flex;
        align-items:center;
        margin-left: 1rem;
        width:5rem;
    }
    .b-radio.radio {
        /*margin-left: .5em;*/
        margin-top: 0px;
    }
</style>
