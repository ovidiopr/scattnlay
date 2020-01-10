<template>
    <div class="field">
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
        watch: {
            // emit updated values
            particle: {
                handler: function () {
                    if (this.particle === 'bulk') {this.layersNum = 1;}
                    if (this.particle === 'core-shell') {this.layersNum = 2;}
                    if (this.particle === 'multilayer') {this.layersNum = 3;}
                }
            },
            layersNum: {
                handler: function () {
                    while (this.layersNum < this.layers.length) {
                        this.layers.pop();
                    }
                    let r_prev = this.layers[0].R;
                    while (this.layersNum > this.layers.length) {
                        // r_prev = r_prev*1.1;
                        this.layers.push(
                            {
                                R: r_prev*0.1,
                                material: 'nk',
                                reN: 4.0,
                                imN: 0.01,
                                index: 1
                            }
                        );
                    }
                    for (let i = 0; i < this.layers.length; i++) {
                        this.layers[i].index = i;
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
</style>
