<template>
  <div>
    <q-tooltip
        v-model = "isShowingTooltip"
        anchor="top middle"
        self="center middle"
    >
      {{formatNumber(localTooltipText,digits+1)}}
    </q-tooltip>
    <q-tooltip v-if="isShowingHelp"
               v-model = "isShowingHelpLocal"
               anchor="bottom middle"
               self="top middle"
    >
      Input example: <b>(1+2)*sqrt(2)</b><br>
    </q-tooltip>
    <q-card
        bordered
        flat
    >
      <q-card-section
          class="items-center bg-grey-2"
          horizontal>
        <div class="side_note text-grey-9 q-px-xs"
             style="width: 4rem"
        >
          {{title}}
        </div>
        <div>
          <q-select
              :model-value="localQSelectValue"
              :options="qSelectOptions"
              bg-color="white"
              dense
              fill-input
              hide-selected
              input-debounce="0"
              options-dense
              style="width: 10rem"
              use-input
              behavior="menu"
              @filter="filterQSelectOptions"
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
                  class="q-py-sm"
              >
                {{formatNumber(localTooltipText,digits)}}
              </div>
            </template>
            <template #option="scope">
              <q-item v-bind="scope.itemProps">
                <q-item-section>
                  {{scope.opt}}
                </q-item-section>
                <q-item-section side>
                  {{formatNumber(evalString(scope.opt),digits)}}
                </q-item-section>
              </q-item>
            </template>
          </q-select>
        </div>
        <div
            class="side_note text-grey-9 q-px-xs"
            style="width: 3rem"
        >
          {{units}}
        </div>
      </q-card-section>
    </q-card>
  </div>
</template>

<script lang="ts">
import {evaluate} from 'mathjs';
import {
  defineComponent,
  ref,
  watch,
  } from 'vue';

export default defineComponent({
  name: 'InputWithUnits',
  props: {
    inputResult: {
      type: Number,
      required: true,
      default: 0
    },
    initialExpression: {
      type: String,
      required: true,
      default: ''
    },
    title: {
      type: String,
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
    'update:input-result',
    'update:is-showing-help'
  ],
  setup(props, {emit}) {

    let localQSelectValue = ref('')
    let localTooltipText = ref('')
    let isShowingTooltip = ref(false)
    let isShowingTooltipAppend = ref(false)
    let isShowingHelpLocal = ref(false)

    let evaluatedValue = ref(0)
    let count_updates = 0
    const digits = 1

    // Set some random values to get correct typing with
    // TypeScript and remove them after initialization.
    let qSelectOptions = ref(['a','b'])
    let qSelectOptionsHistory = ref(['a','b'])
    qSelectOptions.value.pop()
    qSelectOptions.value.pop()
    qSelectOptionsHistory.value.pop()
    qSelectOptionsHistory.value.pop()

    function filterQSelectOptions (val:string,
                       update:(data: ()=>void) => void) {
      update(() => {
        // To remove the selection from previously
        // selected option - we refill the options list
        qSelectOptions.value = qSelectOptionsHistory.value
      })
    }

    function formatNumber (value:string, digits:number):string {
      if (value==='') return ''
      const num = parseFloat(value)
      if ( num < Math.pow(10, -digits) ||
          num > 5*Math.pow(10,  digits+2)
      ) return '='+num.toExponential(digits)
      return '='+Number(Math.round(
              parseFloat(value + 'e' + digits.toString())).toString()
          + 'e-' + digits.toString()).toString()
    }

    // evaluate option items, returns empty string for trivial evaluations and errors
    function evalString(val:string):string {
      // Using try{} block to drop silently invalid input
      try {
        const tryEvaluate = Number(evaluate(val))
        if (!isNaN(tryEvaluate) && tryEvaluate != Number(val)) {
          return tryEvaluate.toString()
        }
      } catch { }
      return ''
    }

    // evaluate current input, keeps the previous evaluateValue for invalid input
    function runEvaluate() {
      // Using try{} block to drop silently invalid input
      try {
        const tryEvaluate = Number(evaluate(localQSelectValue.value))
        if (!isNaN(tryEvaluate)) evaluatedValue.value = tryEvaluate
      } catch { }
    }

    watch(localQSelectValue, () => {
      runEvaluate()
    })

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

    function setTooltip(){
      localTooltipText.value = evaluatedValue.value.toString()
      setTooltipVisibility()
    }

    watch(evaluatedValue, () => {
      emit('update:input-result', evaluatedValue.value)
      setTooltip()
      // Switch off showing help as soon as we have some input from user
      const threshold = 1
      if (count_updates < threshold+1) { // emit only once
        count_updates += 1
        if (props.isShowingHelp && count_updates > threshold) emit('update:is-showing-help', false)
      }
    })

    watch(isShowingHelpLocal, ()=>{
      // isShowingHelpLocal.value is set to be
      // true on hover. Disable it if needed.
      if (isShowingHelpLocal.value) {
        if (qSelectOptions.value.length>0) isShowingHelpLocal.value = false
      }
    })

    watch(qSelectOptions, ()=>{
      if (qSelectOptions.value.length>0) isShowingHelpLocal.value = false
    })

    function handleQSelectBlur(){
      isShowingTooltip.value = false
      const expr = localQSelectValue.value
      if (!qSelectOptionsHistory.value.includes(expr)) qSelectOptionsHistory.value.unshift(expr)
      if (qSelectOptionsHistory.value.length > 5) qSelectOptionsHistory.value.pop()
    }

    // eslint-disable-next-line vue/no-setup-props-destructure
    localQSelectValue.value = props.initialExpression // TODO do we need reactivity from props.initialExpression?
    runEvaluate()
    localTooltipText.value = localQSelectValue.value
    setTooltip()
    isShowingTooltip.value = false

    return {
      localTooltipText, isShowingTooltip,
      isShowingTooltipAppend, isShowingHelpLocal,
      qSelectOptions,  localQSelectValue,
      handleQSelectBlur, filterQSelectOptions,
      digits, formatNumber, evalString
    };
  },
});
</script>
