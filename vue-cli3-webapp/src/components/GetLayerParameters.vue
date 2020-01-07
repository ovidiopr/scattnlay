<template>
    <div class="field">
        <div class="field is-horizontal layer">
            <div class="field-label is-normal">
                <label class="label lnorm">
                    <div v-if="particle==='bulk'"> bulk </div>
                    <div v-else-if="particle==='core-shell'">
                        <div v-if="index == 0"> core </div>
                        <div v-else> shell </div>
                    </div>
                    <div v-else>
                        layer {{index+1}}
                    </div>
                </label>
            </div>
            <div class="field-body">
                <div class="field is-grouped is-grouped-multiline">
                    <div v-if="index == 0" class="rh-input">
                        <input-with-units title="R" v-bind:units="units"
                                          v-bind:value="layer.R"
                                          @newdata="layer.R=$event"/>
                        <span class="tooltiptext">Core radius</span>
                    </div>
                    <div v-else  class="rh-input">
                        <input-with-units title="h" v-bind:units="units"
                                          v-bind:value="layer.R"
                                          @newdata="layer.R=$event"/>
                        <span class="tooltiptext">Layer thickness</span>

                    </div>
                    <b-select v-model="layer.material">
                        <option value="nk">nk-constant</option>
                        <option value="popular">popular</option>
                        <option value="web">refractiveindex.info</option>
                    </b-select>

                </div>
            </div>
        </div>
        <div class="field is-horizontal layer">
            <div class="field-label is-normal">
                <label class="label lnorm">
                    &nbsp;
                </label>
            </div>
            <div class="field-body nk-input">
                <div class="field is-grouped is-grouped-multiline">
                    <input-with-units title="Re(n)" units=""
                                      v-bind:value="layer.reN"
                                      @newdata="layer.reN=$event"
                                      v-bind:isDisabled="isDisabled"
                    />
                    <input-with-units title="Im(n)" units=""
                                      v-bind:value="layer.imN"
                                      @newdata="layer.imN=$event"/>

                </div>
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
        data () {
            return {
                isDisabled: false
            }
        },
        watch: {
            'layer.material': {
                handler: function () {
                    console.log('update');
                    if (this.layer.material == 'nk') {
                        this.isDisabled = false;
                    } else {
                        this.isDisabled = true;
                    }
                }
            }
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
    .lnorm {
        font-weight: normal;
    }
    .rh-input {
        margin-right: 0px;
        position: relative;
        display: inline-block;
    }
    /* Tooltip text */
    .rh-input .tooltiptext {
        visibility: hidden;
        width: 120px;
        background-color: #7957d5;
        color: white;
        text-align: center;
        padding: 5px 0;
        border-radius: 6px;

        /* Position the tooltip text - see examples below! */
        position: absolute;
        z-index: 1;

        bottom: 120%;
        left: 0%;
    }

    /* Show the tooltip text when you mouse over the tooltip container */
    .rh-input:hover .tooltiptext {
        visibility: visible;
    }
    .nk-input {
        margin-bottom: 12px;
    }
</style>
