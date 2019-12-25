<template>
    <div class="field has-addons">
        <p class="control">
            <a class="button is-static input-with-units-title">
                {{ title }}
            </a>
        </p>
        <p class="control">
            <b-input v-model="valueLocal" type="number" step="any"
                     class="input-with-units-value"></b-input>
        </p>
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
            deep: true
        },
        data () {
            return {
                // TODO: Is it OK to modify valueLocal later in <b-input>?
                valueLocal: this.value
            }
        },

        props: ['title', 'units', 'value']
    }
</script>

<style scoped>
.input-with-units-title {
    width:4rem;
}

.input-with-units-value {
    width:6rem;
}

.input-with-units-units {
    width:3rem;
}
</style>
