<template>
  <q-page class="row items-center justify-evenly">
<!--    Your code-->
    <h4>Far-field</h4>
    <input-with-units
        v-model:input-result="someValue"
        v-model:is-showing-help="isShowingHelpForInputWithUnits"
        :initial-expression="someExpr"
        title="Re(n)"
        units="nm"
    ></input-with-units>
    <input-with-units
        v-model:input-result="someValue"
        v-model:is-showing-help="isShowingHelpForInputWithUnits"
        :initial-expression="someExpr"
        title=""
        units=""
        active
    ></input-with-units>
    Input result: {{someValue}}
<!--    tooltip_text="help text"-->
  </q-page>
</template>

<script lang='ts'>
import { defineComponent, ref,
  computed
} from 'vue';
import { useStore } from 'src/store';
import InputWithUnits from 'components/InputWithUnits.vue';


export default defineComponent({
  name: 'PageIndex',
  components: {InputWithUnits },
  setup() {
    const $store = useStore()
    let someValue = ref(10);
    let someExpr = ref('10');
    // InputWithUnits component will disable showing help after first input
    const isShowingHelpForInputWithUnits = computed({
      get: () => $store.state.guiRuntime.isShowingHelpForInputWithUnits,
      set: val => {
        $store.commit('guiRuntime/setIsShowingHelpForInputWithUnits', val)
      }
    })
    // let isShowingHelpForInputWithUnits = ref(true)
    return { someValue, someExpr, isShowingHelpForInputWithUnits};
  }
});
</script>
