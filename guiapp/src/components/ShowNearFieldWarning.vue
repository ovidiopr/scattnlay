<template>
  <div>
    <q-card v-if="showWarning">
      <div class="q-pa-sm">
      <div class="text-h6 q-ml-sm text-warning"><q-icon name="warning"/> Warning</div>
      <div class="q-ma-sm">Near-field evaluation is an experimental feature. In general, it provides a correct result. However, it is not as tested as spectrum computations, especially for the case of large, multilayer, and absorbing spheres. Please, verify the result before using it.</div>
      <div class="q-ma-xs"/>
      <div class="q-ml-sm text-h6"><q-icon size='sm' name="o_info"/> Usage</div>
      <div class="q-ma-sm">Near-field simulation uses settings from spectrum simulation.</div>
              <q-separator inset/>
        <div class="q-ma-sm"/>
        <q-btn icon="settings" flat text-color="primary"
               label="Change spectrum settings"
               @click="void router.push({path: '/spectrum'})"/>
      <q-btn icon="close" flat text-color="primary"
               label="Close"
               @click="showWarning=false"/>
      </div>
    </q-card>
    <div class="q-ma-md"/>
  </div>
</template>

<script lang="ts">
import { useRouter } from 'vue-router'
import {
  defineComponent,
  ref,
} from 'vue'

export default defineComponent({
  name: 'ShowNearFieldWarning',
  setup () {
    const router = useRouter()
    let showWarning = ref(true)

    // import { useQuasar } from 'quasar'
    // const $q = useQuasar()
    // onActivated(()=>{
    //   if (count.value < 2) {
    //     $q.notify({
    //       message: 'Near-field simulation uses settings from spectrum simulation.',
    //       multiLine: true,
    //       position: 'top-right',
    //       icon: 'o_info',
    //       progress: true,
    //       color: 'white',
    //       textColor: 'black',
    //       actions: [
    //         { icon: 'settings', label: 'Change spectrum settings', color: 'primary',
    //           handler: () => { void router.push({path: '/spectrum'}) } },
    //         { icon: 'close', label: 'Close', color: 'primary'},
    //
    //       ],
    //       timeout: 10000
    //     })
    //     count.value += 1
    //   }
    // })
    //

    return { showWarning, router }
  }
})
</script>
