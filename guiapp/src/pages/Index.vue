<template>
  <q-page class="row items-center justify-evenly">

    <input-with-units
        v-model:output-value="someValue"
        v-model:eval-expr="someExpr"
        v-model:is-show-help="isShowHelp"
        title="Re(n)"
        units="nm"
        active
    ></input-with-units>
    Input result: {{someValue}}
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
    let someValue = ref(10);
    let someExpr = ref('10');
    let isShowHelp = ref(true)
    return { someValue, someExpr, isShowHelp };
  }
});
</script>
