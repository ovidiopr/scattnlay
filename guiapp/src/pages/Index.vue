<template>
  <q-page class="row items-center justify-evenly">

    <input-with-units
        v-model:output-value="someValue"
        v-model:eval-expr="someExpr"
        v-model:is-showing-help="isShowingHelp"
        title="Re(n)"
        units="nm"
        active
    ></input-with-units>
    <input-with-units
        v-model:output-value="someValue"
        v-model:eval-expr="someExpr"
        v-model:is-showing-help="isShowingHelp"
        title="Im(n)"
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
    let someExpr = ref('10+2');
    // let isShowingHelp = ref(false)
    let isShowingHelp = ref(true)
    return { someValue, someExpr, isShowingHelp };
  }
});
</script>
