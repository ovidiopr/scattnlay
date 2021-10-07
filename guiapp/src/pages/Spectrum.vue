<template>
  <q-page class="row items-center justify-evenly">

    <input-with-units
        v-model:input-result="someValue"
        v-model:is-showing-help="isShowingHelpInputWithUnits"
        :initial-expression="someExpr"
        title="Re(n)"
        units="nm"
    ></input-with-units>
    <input-with-units
        v-model:input-result="someValue"
        v-model:is-showing-help="isShowingHelpInputWithUnits"
        :initial-expression="someExpr"
        title=""
        units=""
        active
    ></input-with-units>
    Input result: {{someValue}}
<!--    tooltip_text="help text"-->
  </q-page>
</template>

<script lang="ts">
import { defineComponent, ref,
  // computed
} from 'vue';
import { useStore } from 'vuex'
import InputWithUnits from 'components/InputWithUnits.vue';
// import {store} from "quasar/wrappers";
import {storeKey} from "src/store";


export default defineComponent({
  name: 'PageIndex',
  components: {InputWithUnits },
  setup() {
    const $store = useStore(storeKey)
    console.log( 'store value',$store.state.example.prop)
    let someValue = ref(10);
    let someExpr = ref('10');
    // InputWithUnits component will disable showing help after first input
    let isShowingHelpInputWithUnits = ref(true)
    return { someValue, someExpr, isShowingHelpInputWithUnits };
  }
});
</script>
