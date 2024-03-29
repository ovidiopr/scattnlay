<template>
  <div class="q-pa-xs">
    <q-tooltip
      v-model="isShowingTooltip"
      anchor="top middle"
      self="center middle"
    >
      {{ formatNumber(localTooltipText, inputWithUnitsTooltipDigits) }}
    </q-tooltip>
    <q-tooltip
      v-if="isShowingHelp && !isInfoMode"
      v-model="isShowingHelpLocal"
      anchor="bottom middle"
      self="center middle"
    >
      Input example: <b>{{ helpExpr }}</b
      ><br />
    </q-tooltip>
    <q-tooltip v-if="isInfoMode" anchor="top middle" self="center middle">
      current settings<br />
    </q-tooltip>
    <q-card bordered flat>
      <q-card-section class="items-center bg-grey-2" horizontal>
        <div
          class="side_note text-grey-9 q-px-xs text-center"
          :style="'width: ' + inputWithUnitsTitleWidthStyle"
        >
          {{ title }}
        </div>
        <div>
          <q-select
            :model-value="localQSelectModel"
            :options="qSelectOptions"
            :disable="isInfoMode"
            bg-color="white"
            dense
            fill-input
            hide-selected
            input-debounce="0"
            options-dense
            :style="'width: ' + inputWithUnitsBodyWidthStyle"
            use-input
            input-class="q-pl-xs"
            behavior="menu"
            @filter="filterQSelectOptions"
            @blur="handleQSelectBlur"
            @keydown.enter="handleQSelectBlur"
            @input-value="localQSelectModel = $event"
          >
            <template v-if="isError" #prepend>
              <q-tooltip> Input conflict </q-tooltip>
              <q-icon name="error" class="text-warning q-pl-xs" />
            </template>
            <template
              v-if="isShowingTooltipAppend && !isShowingTooltip"
              #append
            >
              <div style="font-size: 12px" class="q-py-sm">
                {{ formatNumber(localTooltipText, inputWithUnitsInlineDigits) }}
              </div>
            </template>
            <template #option="scope">
              <q-item v-bind="scope.itemProps" class="q-px-sm">
                <q-item-section>
                  {{ scope.opt }}
                </q-item-section>
                <q-item-section side>
                  {{
                    formatNumber(
                      evalString(scope.opt),
                      inputWithUnitsInlineDigits
                    )
                  }}
                </q-item-section>
              </q-item>
            </template>
          </q-select>
        </div>
        <div
          class="side_note text-grey-9 q-px-xs text-center"
          :style="'width: ' + inputWithUnitsUnitsWidthStyle"
        >
          {{ units }}
        </div>
      </q-card-section>
    </q-card>
  </div>
</template>

<script lang="ts">
import { evaluate } from 'mathjs';
import { defineComponent, ref, watch } from 'vue';
import {
  inputWithUnitsTitleWidthStyle,
  inputWithUnitsBodyWidthStyle,
  inputWithUnitsUnitsWidthStyle,
  inputWithUnitsHistoryLength,
  inputWithUnitsInlineDigits,
  inputWithUnitsTooltipDigits,
} from 'components/config';

export default defineComponent({
  name: 'InputWithUnits',
  props: {
    inputResult: {
      type: Number,
      required: true,
      default: 0,
    },
    initialExpression: {
      type: String,
      required: true,
      default: '',
    },
    title: {
      type: String,
      default: '',
    },
    units: {
      type: String,
      default: '',
    },
    isShowingHelp: {
      type: Boolean,
      default: false,
    },
    isError: {
      type: Boolean,
      default: false,
    },
    isInfoMode: {
      type: Boolean,
      default: false,
    },
  },
  emits: ['update:input-result', 'update:is-showing-help'],
  setup(props, { emit }) {
    let localQSelectModel = ref('');
    let localTooltipText = ref('');
    let isShowingTooltip = ref(false);
    let isShowingTooltipAppend = ref(false);
    let isShowingHelpLocal = ref(false);
    let helpExpr = ref('(1+2)*sqrt(2)');

    let evaluated = ref(0);
    let count_updates = 0;

    // Set some random values to get correct typing with
    // TypeScript and remove them after initialization.
    let qSelectOptions = ref(['a', 'b']);
    let qSelectOptionsHistory = ref(['a', 'b']);
    qSelectOptions.value.pop();
    qSelectOptions.value.pop();
    qSelectOptionsHistory.value.pop();
    qSelectOptionsHistory.value.pop();

    // evaluate current input, keeps the previous evaluateValue for invalid input
    function runEvaluate() {
      // Using try{} block to drop silently invalid input
      try {
        const tryEvaluate = Number(evaluate(localQSelectModel.value));
        if (!isNaN(tryEvaluate)) evaluated.value = tryEvaluate;
      } catch {}
    }

    watch(localQSelectModel, () => {
      runEvaluate();
    });

    function setTooltipVisibility() {
      if (evaluated.value != Number(localQSelectModel.value)) {
        isShowingTooltip.value = true;
        isShowingTooltipAppend.value = true;
      } else {
        isShowingTooltip.value = false;
        isShowingTooltipAppend.value = false;
      }
    }

    watch(isShowingTooltip, () => {
      // For a trivial case we would like to switch off showing tooltip
      if (isShowingTooltip.value) setTooltipVisibility();
    });

    function setTooltip() {
      localTooltipText.value = evaluated.value.toString();
      setTooltipVisibility();
    }

    watch(evaluated, () => {
      emit('update:input-result', evaluated.value);
      setTooltip();
      // Switch off showing help as soon as we have some input from user
      const threshold = 1;
      if (count_updates < threshold + 1) {
        // limit the unbound grow of count_updates
        count_updates += 1;
        if (props.isShowingHelp && count_updates > threshold) {
          qSelectOptionsHistory.value.unshift(helpExpr.value);
          emit('update:is-showing-help', false);
        }
      }
    });

    watch(isShowingHelpLocal, () => {
      // isShowingHelpLocal.value is set to be
      // true on hover. Disable it if needed.
      if (isShowingHelpLocal.value) {
        if (qSelectOptions.value.length > 0) isShowingHelpLocal.value = false;
      }
    });

    watch(qSelectOptions, () => {
      if (qSelectOptions.value.length > 0) isShowingHelpLocal.value = false;
    });

    watch(props, () => {
      // Using try{} block to drop silently invalid input
      try {
        // If props.inputResults changed and is not equal to local
        // expression localQSelectModel.value, then update local expression
        const tryEvaluate = Number(evaluate(localQSelectModel.value));
        if (!isNaN(tryEvaluate) && props.inputResult != tryEvaluate) {
          localQSelectModel.value = props.inputResult.toString();
        }
      } catch {}
    });

    localQSelectModel.value = props.initialExpression.toString();
    runEvaluate();
    localTooltipText.value = localQSelectModel.value;
    setTooltip();
    isShowingTooltip.value = false;

    return {
      inputWithUnitsTitleWidthStyle,
      inputWithUnitsBodyWidthStyle,
      inputWithUnitsUnitsWidthStyle,
      inputWithUnitsHistoryLength,
      inputWithUnitsInlineDigits,
      inputWithUnitsTooltipDigits,
      localTooltipText,
      isShowingTooltip,
      isShowingTooltipAppend,
      isShowingHelpLocal,
      helpExpr,
      qSelectOptions,
      localQSelectModel,

      handleQSelectBlur() {
        isShowingTooltip.value = false;
        const expr = localQSelectModel.value;
        if (!qSelectOptionsHistory.value.includes(expr))
          qSelectOptionsHistory.value.unshift(expr);
        if (qSelectOptionsHistory.value.length > inputWithUnitsHistoryLength)
          qSelectOptionsHistory.value.pop();
      },

      filterQSelectOptions(val: string, update: (data: () => void) => void) {
        update(() => {
          // To remove the selection from previously
          // selected option - we refill the options list
          qSelectOptions.value = qSelectOptionsHistory.value;
        });
      },

      formatNumber(value: string, digits: number, prepend: string): string {
        if (!prepend) prepend = '=';
        if (value === '') return '';
        const num = parseFloat(value);
        if (num < Math.pow(10, -digits) || num > 5 * Math.pow(10, digits + 2))
          return prepend + num.toExponential(digits);
        return (
          prepend +
          Number(
            Math.round(parseFloat(value + 'e' + digits.toString())).toString() +
              'e-' +
              digits.toString()
          ).toString()
        );
      },

      // evaluate option items, returns empty string for trivial evaluations and errors
      evalString(val: string): string {
        // Using try{} block to drop silently invalid input
        try {
          const tryEvaluate = Number(evaluate(val));
          if (!isNaN(tryEvaluate) && tryEvaluate != Number(val)) {
            return tryEvaluate.toString();
          }
        } catch {}
        return '';
      },
    };
  },
});
</script>
