<template>
  <q-page class="row items-center justify-evenly">

    <input-with-units
        v-model:input-result="someValue"
        v-model:initial-expression="someExpr"
        v-model:is-showing-help="isShowingHelpInputWithUnits"
        title="Re(n)"
        units="nm"
        active
    ></input-with-units>
    <input-with-units
        v-model:input-result="someValue"
        v-model:initial-expression="someExpr"
        v-model:is-showing-help="isShowingHelpInputWithUnits"
        title=""
        units=""
        active
    ></input-with-units>
    Input result: {{someValue}}
<!--    tooltip_text="help text"-->
  </q-page>
</template>

<script lang="ts">
import { defineComponent, ref } from 'vue';
import InputWithUnits from 'components/InputWithUnits.vue';

export default defineComponent({
  name: 'PageIndex',
  components: {InputWithUnits },
  setup() {
    let someValue = ref(10);
    let someExpr = ref('10+3');
    // InputWithUnits component will disable showing help after first input
    let isShowingHelpInputWithUnits = ref(true)
    return { someValue, someExpr, isShowingHelpInputWithUnits };
  }
});
</script>
