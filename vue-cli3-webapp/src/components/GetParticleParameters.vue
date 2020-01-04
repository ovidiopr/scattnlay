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
        <GetLayerParameters v-bind:layer="layers[0]"
                            v-bind:units="units"
                            v-bind:particle="particle"
                            v-bind:index="index"
        />
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
            deep: true
        },
        props: ['layers', 'units']
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
