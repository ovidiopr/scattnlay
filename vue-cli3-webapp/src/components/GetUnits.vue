<template>
    <div class="field is-horizontal">
        <div class="field-label is-normal">
            <label class="label">Units</label>
        </div>
        <div class="field-body">
            <div class="columns">
                <div class="column">
                    <b-select v-model="unitsLocal">
                        <option value="nm">nm</option>
                        <option value="mkm">mkm</option>
                        <option value="mm">mm</option>
                        <option value="cm">cm</option>
                        <option value="m">m</option>
                        <option value="km">km</option>
                    </b-select>
                </div>
                <div class="column">
                    <b-checkbox v-model="isSourceOtherUnitsLocal"> use other units for source</b-checkbox>
                </div>
                <div class="column">
                    <div v-if="isSourceOtherUnitsLocal">
                        <b-select v-model="source_unitsLocal">
                            <optgroup label="Frequency">
                                <option value="Hz">Hz</option>
                                <option value="kHz">kHz</option>
                                <option value="MHz">MHz</option>
                                <option value="GHz">GHz</option>
                                <option value="THz">THz</option>
                            </optgroup>
                                <optgroup label="Energy">
                                    <option value="meV">meV</option>
                                    <option value="eV">eV</option>
                                </optgroup>
                                <optgroup label="Period duration">
                                    <option value="ps">ps</option>
                                    <option value="fs">fs</option>
                                </optgroup>
                        </b-select>
                    </div>
                    <div v-else>
                        <b-select v-model="source_unitsLocal" disabled>
                            <option value="nm">nm</option>
                            <option value="mkm">mkm</option>
                            <option value="mm">mm</option>
                            <option value="cm">cm</option>
                            <option value="m">m</option>
                            <option value="km">km</option>
                        </b-select>
                    </div>
                </div>
            </div>
        </div>
    </div>
</template>

<script>
    export default {
        name: "GetUnits",
        data () {
            return {
                unitsLocal: this.units,
                source_unitsLocal: this.source_units,
                isSourceOtherUnitsLocal: this.isSourceOtherUnits
            }
        },
        watch: {
            // emit updated values
            unitsLocal: {
                handler: function () {
                    this.$emit('unitsData',this.unitsLocal);
                }
            },
            source_unitsLocal: {
                handler: function () {
                    this.$emit('source_unitsData',this.source_unitsLocal);
                }
            },
            isSourceOtherUnitsLocal: {
                handler: function () {
                    this.$emit('isSourceOtherUnitsData',this.isSourceOtherUnitsLocal);
                }
            },
            // update local values
            units: {
                handler: function () {
                    this.unitsLocal = this.units;
                }
            },
            source_units: {
                handler: function () {
                    this.source_unitsLocal = this.source_units;
                }
            },
            isSourceOtherUnits: {
                handler: function () {
                    this.isSourceOtherUnitsLocal = this.isSourceOtherUnits;
                }
            },
            deep: true
        },
        props: ['units', 'source_units', 'isSourceOtherUnits']
    }
</script>

<style scoped>
    .checkbox {
        margin-top: 0px;
    }
</style>
