<template>
    <div class="field is-horizontal">
        <div class="field-label is-normal">
            <label class="label">Host media</label>
        </div>
        <div class="field-body">
            <div class="field is-grouped is-grouped-multiline">
                <input-with-units title="Re(n)" units=""
                                  v-bind:value="hostIndexLocal"
                                  @newdata="hostIndexLocal=$event"/>
            </div>
        </div>
    </div>
</template>

<script>
    import InputWithUnits from "./InputWithUnits.vue";
    export default {
        name: "GetHostIndex",
        components: {
            InputWithUnits
        },
        data () {
            return {
                // TODO: Is it OK to modify Local later?
                hostIndexLocal: this.hostIndex,
            }
        },
        watch: {
            // emit updated values
            hostIndexLocal: {
                handler: function () {
                    this.$emit('hostIndexData',this.hostIndexLocal);
                }
            },
            // update local values
            hostIndex: {
                handler: function () {
                    this.hostIndexLocal = this.hostIndex;
                }
            },
            deep: true
        },
        props: ['hostIndex']
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
