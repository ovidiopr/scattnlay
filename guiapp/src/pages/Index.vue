<template>
  <q-page class="row items-center justify-evenly">

    <input-with-units
        v-model:name="name"
        title="Re(n)"
        units="nm"
        tooltip_text="help text"
        active
    ></input-with-units>
    Input result: {{name}}
<!--    tooltip_text="help text"-->
  </q-page>
</template>

<script lang="ts">
// import { Todo, Meta } from 'components/models';
import { defineComponent, ref } from 'vue';
import InputWithUnits from 'components/CompositionComponent.vue';

export default defineComponent({
  name: 'PageIndex',
  components: {InputWithUnits },
  setup() {
    let name = ref('')
    return {name};
  }
});
</script>
