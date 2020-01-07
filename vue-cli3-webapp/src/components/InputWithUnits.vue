<template>
    <div class="field has-addons">
        <p class="control">
            <a class="button is-static input-with-units-title">
                {{ title }}
            </a>
        </p>
        <div v-if="isDisabledLocal">
            <p class="control">
            <b-input v-model="valueLocal" type="number" step="any"
                 class="input-with-units-value" disabled/>

            </p>
        </div>
        <div v-else>
            <p class="control">
                <b-input v-model="valueLocal" type="number" step="any"
                         class="input-with-units-value"/>

            </p>
        </div>
        <p class="control">
            <a class="button is-static input-with-units-units">
                 {{ units }}
            </a>
        </p>
    </div>
</template>

<script>
    export default {
        name: "InputWithUnits",
        watch: {
            valueLocal: {
                handler: function () {
                    this.$emit('newdata',this.valueLocal);
                }
            },
            value: {
                handler: function () {
                    this.valueLocal = this.value;
                }
            },
            isDisabled: {
                handler: function () {
                    this.isDisabledLocal = this.isDisabled;
                }
            },
            deep: true
        },
        data () {
            return {
                // TODO: Is it OK to modify valueLocal later in <b-input>?
                valueLocal: this.value,
                isDisabledLocal: false
            }
        },

        props: ['title', 'units', 'value', 'isDisabled']
    }
</script>

<style scoped>
.input-with-units-title {
    width:4rem;
    z-index: 1;
}

.input-with-units-value {
    width:6rem;
}
.input-with-units-units {
    width:3rem;
}
</style>
