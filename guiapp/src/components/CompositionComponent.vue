<template>
  <div>
    <q-card
            flat
            bordered
    >
      <q-tooltip v-if="tooltipText">
        {{tooltip_text}}
      </q-tooltip>
      <q-card-section
          horizontal
          class="items-center bg-grey-2">
        <div class="side_note text-grey-9 q-px-xs"
             style="width: 4rem"
        >
          {{title}}
        </div>

        <q-select
            :model-value="name"
            @input-value="$emit('update:name', $event)"
            use-input
            hide-selected
            fill-input
            input-debounce="0"
            :options="options"
            style="width: 10rem"
            options-dense
            dense
            bg-color="white"
            class="q-py-none"
        >
          <template v-slot:no-option>
            <q-item>
              <q-item-section class="text-grey">
                No results
              </q-item-section>
            </q-item>
          </template>
        </q-select>

        <div
            class="side_note text-grey-9 q-px-xs"
            style="width: 3rem"
        >
          {{ units}}
        </div>
      </q-card-section>
    </q-card>
    local options: {{name}}
  </div>
</template>

<script lang="ts">
// import { useModelWrapper } from 'components/modelWrapper'
import {
  defineComponent,
  //   watch,
  // PropType,
  // computed,
  ref,
  // toRef,
  // Ref,
} from 'vue';
// import { Todo, Meta } from './models';

// function useModelWrapper(props, emit, name = 'modelValue') {
//   return computed({
//     get: () => props[name],
//     set: (value) => emit(`update:${name}`, value)
//   })
// }
const stringOptions = [
  'Google', 'Facebook', 'Twitter', 'Apple', 'Oracle'
]

export default defineComponent({
  name: 'InputWithUnits',
  props: {
    name: {
      type: String,
      default: ''
    },
    title: {
      type: String,
      required: true,
      default: ''
    },
    units: {
      type:String,
      default: ''
    },
    tooltipText: {
      type: String,
      default: ''
    }
  },
  emits : [
    'update:name'
  ],
  setup(props, {emit}) {
  //   // import { ref } from 'vue'
  //   // const return_value = ref(0)
  //   // var options = ref([3, 2])
    const options = ref(stringOptions)
  //   // var return_value = ref(null)
  //   var model = ref(null)
    return {
  //     model,
      options
    }
  },
});
</script>
