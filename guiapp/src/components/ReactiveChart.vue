<template>
    <div :ref="uuidLocal"></div>
</template>
<script lang="ts">
import {Data, Layout, newPlot, react} from 'plotly.js'
import {
  defineComponent,
    PropType,
  ref,
    onMounted
  // watch,
} from 'vue'

export default defineComponent({
  name: 'ReactiveChart',
  props: {
    uuid: {
      type: String,
      required: true,
    },
    traces: {
      type: Object as PropType<Data[]>,
      required: true
    },
    layout: {
      type: Object as PropType<Layout>,
      required: true
    },
  },
  // emits: [
  //   'update:input-result',
  //   'update:is-showing-help'
  // ],
  setup(props,/* {emit} */) {
    const uuidLocal = ref(props.uuid)
    onMounted(()=>
        newPlot(uuidLocal.value, props.traces,
            props.layout,
            {responsive: true, showSendToCloud: true, displaylogo: false}
        )
    )

      // watch: {
      //   chart: {
      //     handler: function () {
      //       react(
      //           this.$refs[uuid],
      //           traces,
      //           layout
      //       );
      //     },
      //     deep: true
      //   }
      // }
    return {
      uuidLocal
    }
    }
})
</script>
