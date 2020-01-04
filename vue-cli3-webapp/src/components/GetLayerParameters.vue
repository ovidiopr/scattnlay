<template>
    <div class="field is-horizontal layer">
        <div class="field-label is-normal">
            <label class="label">
                <div v-if="particle==='bulk'"> bulk </div>
                <div v-else-if="particle==='core-shell'">
                    <div v-if="index == 0"> core </div>
                    <div v-else> shell </div>
                </div>
                <div v-else> layer {{index+1}} </div>
            </label>
        </div>
        <div class="field-body">
            <div class="field is-grouped is-grouped-multiline">
                <input-with-units title="R" v-bind:units="units"
                                  v-bind:value="layer.R"
                                  @newdata="layer.R=$event"/>
                <input-with-units title="Re(n)" units=""
                                  v-bind:value="layer.reN"
                                  @newdata="layer.reN=$event"/>
                <input-with-units title="Im(n)" units=""
                                  v-bind:value="layer.imN"
                                  @newdata="layer.imN=$event"/> 
            </div>
        </div>
    </div>

</template>

<script>
    import InputWithUnits from "./InputWithUnits.vue";
    export default {
        name: "GetLayerParameters",
        components: {
            InputWithUnits
        },
        props: ['layer', 'particle', 'index', 'units']
    }
</script>

<style scoped>
    .binput {
        display:flex;
        align-items:center;
        margin-left: 1rem;
        width:5rem;
    }
    .layer {
        margin: 0.5rem;
    }
</style>
