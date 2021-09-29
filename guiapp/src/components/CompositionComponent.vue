<template>
  <div>
    <q-card
        bordered
        flat
    >
      <q-card-section
          class="items-center bg-grey-2"
          horizontal>
        <div v-if='title' class="side_note text-grey-9 q-px-xs"
             style="width: 4rem"
        >
          {{ title }}
        </div>
        <div
            @mousemove="setTooltip"
        >

        <q-tooltip
            v-model = "isShowingTooltip"
                   anchor="top middle"
                   self="center middle"
        >
          = {{round_number(localTooltipText,5)}}
        </q-tooltip>

        <q-select
            :model-value="localQSelectValue"
            :options="options"
            bg-color="white"
            class="q-py-none"
            dense
            fill-input
            hide-selected
            input-debounce="0"
            options-dense
            style="width: 10rem"
            use-input
            behavior="menu"
            @focus="setTooltip"
            @blur="handleQSelectBlur"
            @keydown.enter="handleQSelectBlur"
            @keydown.tab="handleQSelectBlur"
            @input-value="localQSelectValue=$event"
        >
          <template
                    v-if="localTooltipText&&!isShowingTooltip"
                    #append
          >
            <div
                style="font-size: 12px"
            >
              ={{round_number(localTooltipText,1)}}
            </div>
          </template>
<!--          <template #no-option>-->
<!--            <q-item>-->
<!--              <q-item-section class="text-grey">-->
<!--                No results-->
<!--              </q-item-section>-->
<!--            </q-item>-->
<!--          </template>-->
        </q-select>
        </div>

        <div
            class="side_note text-grey-9 q-px-xs"
            style="width: 3rem"
        >
          {{ units }}
        </div>
      </q-card-section>
    </q-card>
<!--    local options: {{ inputValue }}-->

    local: {{localQSelectValue}}
  </div>
</template>

<script lang="ts">
// import { useModelWrapper } from 'components/modelWrapper'
import {evaluate, isNumeric} from 'mathjs';
import {defineComponent,
  ref,
  watch,
  onMounted,
  } from 'vue';

export default defineComponent({
  name: 'InputWithUnits',
  props: {
    outputValue: {
      type: Number,
      default: 0
    },
    evalExpr: {
      type: String,
      default: ''
    },
    title: {
      type: String,
      required: true,
      default: ''
    },
    units: {
      type: String,
      default: ''
    },
    tooltipText: {
      type: String,
      default: ''
    }
  },
  emits: [
    'update:output-value'
    // <!--            @input-value="$emit('update:input-expr', $event)"-->
    // <!--        >-->
  ],
  setup(props, {emit}) {
    let localQSelectValue = ref('')
    let localTooltipText = ref('')
    let isShowingTooltip = ref(false)


    onMounted(()=>{
      localQSelectValue.value = props.evalExpr
      localTooltipText.value = props.evalExpr
      if (props.tooltipText) localTooltipText.value = props.tooltipText
      console.log(localTooltipText.value)
    });


    let options = ref(['1+2', 'sqrt(6)*20'])
    let evaluatedValue = ref(0)


    function setTooltip(){
      if (evaluatedValue.value != Number(localQSelectValue.value)) {
        localTooltipText.value = evaluatedValue.value.toString()
        isShowingTooltip.value = true
      } else {
        localTooltipText.value = ''
        isShowingTooltip.value = false
      }
    }


    function handleQSelectBlur(){
      isShowingTooltip.value=false
      const expr = localQSelectValue.value
      if (!options.value.includes(expr)) options.value.unshift(expr);
      if (options.value.length > 5) options.value.pop();
      // localQSelectValue.value = evaluatedValue.value.toString()
    }


    function round_number (value:string, digits:number):number {
      // TODO manage special cases of very big and very
      // small numbers assuming digits value is small
      return Number(Math.round(parseFloat(value + 'e' + digits.toString())) + 'e-' + digits.toString())
    }


    watch(localQSelectValue, () => {
      let tryEvaluate:any
      try {
        tryEvaluate = evaluate(localQSelectValue.value)
        if (isNumeric(tryEvaluate)) evaluatedValue.value = Number(tryEvaluate)
      } catch { }
    })




    watch(evaluatedValue, () => {
      emit('update:output-value', evaluatedValue.value)
      setTooltip()
    })


    return {
      options,  localQSelectValue, localTooltipText, isShowingTooltip,
      setTooltip, handleQSelectBlur, round_number
    };
  },
});
</script>
