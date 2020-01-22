<template>
    <div class="field is-horizontal source-params">
        <div class="field-label is-normal">
            <label class="label">
                <div v-if="source_units.endsWith('Hz')"> Frequency  </div>
                <div v-else-if="source_units.endsWith('eV')"> Energy  </div>
                <div v-else-if="source_units.endsWith('s')"> Period  </div>
                <div v-else>               Wavelength              </div>
            </label>
        </div>
        <div class="field-body">
            <div class="field is-grouped is-grouped-multiline">
                <input-with-units title="from" v-bind:units="source_units"
                                  v-bind:value="fromWLLocal"
                                  @newdata="fromWLLocal=$event"/>
                <input-with-units title="to" v-bind:units="source_units"
                                  v-bind:value="toWLLocal"
                                  @newdata="toWLLocal=$event"/>
                <div>
                <div v-if="isStep">
                <input-with-units title="step" v-bind:units="source_units"
                                  v-bind:value="stepWLLocal"
                                  @newdata="stepWLLocal=$event"/>
                </div>
                <div v-else>
                    <input-with-units title="points" units=""
                                      v-bind:value="pointsNum"
                                      @newdata="pointsNum=$event"/>
                </div>
                <b-switch v-model="isStep" class="switch-style"/>
                </div>
            </div>
        </div>
    </div>
</template>

<script>
    import InputWithUnits from "./InputWithUnits.vue";
    export default {
        name: "GetSourceParameters",
        components: {
            InputWithUnits
        },
        data () {
            return {
                // TODO: Is it OK to modify Local later?
                fromWLLocal: this.fromWL,
                toWLLocal: this.toWL,
                stepWLLocal: this.stepWL,
                pointsNum: 100,
                isStep: true
            }
        },
        methods: {

        },
        watch: {
            isStep: {
                handler: function () {
                    if (this.isStep) {
                        this.stepWLLocal = (this.toWLLocal-this.fromWLLocal)/this.pointsNum;
                    } else {
                        this.pointsNum = Math.round((this.toWLLocal-this.fromWLLocal)/this.stepWLLocal);
                    }
                }
            },
            // emit updated values
            fromWLLocal: {
                handler: function () {
                    this.fromWLLocal = parseFloat(this.fromWLLocal);
                    // if (this.fromWLLocal > this.toWLLocal) {
                    //     this.fromWLLocal = this.toWLLocal;
                    // }
                    this.$emit('fromWLData',this.fromWLLocal);
                }
            },
            toWLLocal: {
                handler: function () {
                    this.toWLLocal = parseFloat(this.toWLLocal);
                    // if (this.toWLLocal < this.fromWLLocal) {
                    //     this.toWLLocal = this.fromWLLocal;
                    // }
                    this.$emit('toWLData',this.toWLLocal);
                }
            },
            stepWLLocal: {
                handler: function () {
                    this.$emit('stepWLData',this.stepWLLocal);
                }
            },
            // update local values
            fromWL: {
                handler: function () {
                    this.fromWLLocal = this.fromWL;
                }
            },
            toWL: {
                handler: function () {
                    this.toWLLocal = this.toWL;
                }
            },
            stepWL: {
                handler: function () {
                    this.stepWLLocal = this.stepWL;
                }
            },
            deep: true
        },
        props: ['fromWL', 'toWL', 'stepWL', 'source_units']
    }
</script>

<style scoped>
    .source-params {
        padding-bottom: 5px;
    }
    .switch-style {
        padding: 5px;
    }
</style>
