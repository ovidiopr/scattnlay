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
                                      v-bind:value="layers[0].R"
                                      @newdata="layers[0].R=$event"/>
                    <input-with-units title="Re(n)" units=""
                                      v-bind:value="layers[0].reN"
                                      @newdata="layers[0].reN=$event"/>
                    <input-with-units title="Im(n)" units=""
                                      v-bind:value="layers[0].imN"
                                      @newdata="layers[0].imN=$event"/>
                </div>
            </div>
        </div>
    </div>
</template>

<script>
    import InputWithUnits from "./InputWithUnits.vue";
    export default {
        name: "GetParticleParameters",
        components: {
            InputWithUnits
        },
        data () {
            return {
                // TODO: Is it OK to modify Local later?
                particle: 'bulk',
                layersNum: 1
            }
        },
        props: ['layers', 'units']
    }
</script>

<style scoped>
</style>
