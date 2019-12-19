<template>
    <div class="field has-addons">
        <p class="control">
            <a class="button is-static input-with-units-title">
                {{ title }}
            </a>
        </p>
        <p class="control">
            <b-input v-model="value" type="number" step="any"
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
            value: {
                handler: function () {
                    this.$emit('newdata',this.value);
                }
            },
            deep: true
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