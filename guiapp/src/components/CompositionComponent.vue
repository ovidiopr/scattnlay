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
        <div>

        <q-tooltip
            v-model = "isShowingTooltip"
                   anchor="top middle"
                   self="center middle"
        >
          {{format_number(localTooltipText,2)}}
        </q-tooltip>
          <q-tooltip v-if="isShowingHelp"
              v-model = "isShowingHelpLocal"
              anchor="bottom middle"
              self="top middle"
          >
              Input example: <b>(1+2)*sqrt(2)</b><br>
          </q-tooltip>

        <q-select
            :model-value="localQSelectValue"
            :options="options"
            bg-color="white"
            dense
            fill-input
            hide-selected
            input-debounce="0"
            options-dense
            style="width: 10rem"
            use-input
            behavior="menu"
            @filter="filterFn"
            @blur="handleQSelectBlur"
            @keydown.enter="handleQSelectBlur"
            @input-value="localQSelectValue=$event"
        >
          <template
                    v-if="isShowingTooltipAppend&&!isShowingTooltip"
                    #append
          >
            <div
                style="font-size: 12px"
                class="q-py-xl"
            >
              {{format_number(localTooltipText,1)}}
            </div>
          </template>
          <template #option="scope">
            <q-item v-bind="scope.itemProps">
              <q-item-section>
                {{scope.opt}}
              </q-item-section>
              <q-item-section side>
                {{format_number(evalString(scope.opt),1)}}
              </q-item-section>
            </q-item>
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
import {evaluate} from 'mathjs';
import {defineComponent,
  ref,
  watch,
  onMounted
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
    isShowingHelp: {
      type: Boolean,
      default: false
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
    let isShowingTooltipAppend = ref(false)
    let isShowingHelpLocal = ref(false)

    let options = ref(['a','b'])
    let optionsHistory = ref(['a','b'])
    let evaluatedValue = ref(0)

    function filterFn (val:string,
                       update:(data: ()=>void) => void) {
      update(() => {
        // To remove the selection from previously
        // selection option - we recreate the options list
        options.value = optionsHistory.value
      })
    }

    function runEvaluate() {
      try {
        const tryEvaluate = Number(evaluate(localQSelectValue.value))
        if (!isNaN(tryEvaluate)) evaluatedValue.value = tryEvaluate
      } catch { }
    }

    function evalString(val:string):string {
      try {
        const tryEvaluate = Number(evaluate(val))
        if (!isNaN(tryEvaluate) && tryEvaluate != Number(val)) {
          return tryEvaluate.toString()
        }
      } catch { }
      return ''
    }

    onMounted(()=>{
      isShowingTooltip.value = false
    });


    function setTooltipVisibility(){
      if (evaluatedValue.value != Number(localQSelectValue.value)) {
        isShowingTooltip.value = true
        isShowingTooltipAppend.value = true
      } else {
        isShowingTooltip.value = false
        isShowingTooltipAppend.value = false
      }
    }

    watch(isShowingTooltip, ()=>{
      // For a trivial case we would like to switch off showing tooltip
      if (isShowingTooltip.value) setTooltipVisibility()
    })

    watch(isShowingHelpLocal, ()=>{
      // isShowingHelpLocal.value is set to be
      // true on hover. Disable it if needed.
      if (isShowingHelpLocal.value) {
        if (options.value.length>0) isShowingHelpLocal.value = false
      }
    })


    watch(options, ()=>{
      if (options.value.length>0) isShowingHelpLocal.value = false
    })


    function setTooltip(){
      localTooltipText.value = evaluatedValue.value.toString()
      setTooltipVisibility()
    }

    function handleQSelectBlur(){
      isShowingTooltip.value = false
      isShowingHelpLocal.value = false
      const expr = localQSelectValue.value
      if (!optionsHistory.value.includes(expr)) optionsHistory.value.unshift(expr)
      if (optionsHistory.value.length > 5) optionsHistory.value.pop()
    }



    function format_number (value:string, digits:number):string {
      if (value==='') return ''
      // const help_string = 'text\n text'
      const num = parseFloat(value)
      if ( num < Math.pow(10, -digits) ||
           num > 5*Math.pow(10,  digits+2)
         ) return '='+num.toExponential(digits)
      return '='+Number(Math.round(
          parseFloat(value + 'e' + digits.toString())).toString()
          + 'e-' + digits.toString()).toString()
    }


    watch(localQSelectValue, () => {
      runEvaluate()
    })




    watch(evaluatedValue, () => {
      emit('update:output-value', evaluatedValue.value)
      setTooltip()
    })


    // eslint-disable-next-line vue/no-setup-props-destructure
    localQSelectValue.value = props.evalExpr // TODO do we need reactivity from props.evalExpr?
    runEvaluate()
    localTooltipText.value = localQSelectValue.value
    setTooltip()
    options.value.pop()
    options.value.pop()
    optionsHistory.value.pop()
    optionsHistory.value.pop()
    optionsHistory.value.unshift(localQSelectValue.value)

    return {
      options,  localQSelectValue, localTooltipText, isShowingTooltip,
      isShowingHelpLocal, isShowingTooltipAppend,
      setTooltip, handleQSelectBlur, format_number,
      filterFn, evalString
    };
  },
});
</script>
