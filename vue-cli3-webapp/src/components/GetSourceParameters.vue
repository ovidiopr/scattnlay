<template>
    <div class="field is-horizontal">
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
                <input-with-units title="step" v-bind:units="source_units"
                                  v-bind:value="stepWLLocal"
                                  @newdata="stepWLLocal=$event"/>
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
                stepWLLocal: this.stepWL
            }
        },
        watch: {
            // emit updated values
            fromWLLocal: {
                handler: function () {
                    this.$emit('fromWLData',this.fromWLLocal);
                }
            },
            toWLLocal: {
                handler: function () {
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
</style>
