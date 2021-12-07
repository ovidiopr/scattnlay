<template>
  <q-layout view="lHh lpR fFf">
<!--  <q-layout view="hHh Lpr lFr">-->
    <q-header elevated>
      <q-toolbar>
        <q-btn
          flat
          dense
          round
          icon="menu"
          aria-label="Menu"
          @click="toggleLeftDrawer"
        />
        <q-tabs align="right">
          <q-route-tab to="/spectrum" label="Spectrum" name="spectrum"/>
          <q-route-tab to="/nearfield" label="Near-field" name="nearfield"/>
<!--          <q-route-tab to="/farfield" label="Far-field" name="farfield"/>-->
        </q-tabs>
      </q-toolbar>
    </q-header>

    <q-drawer
      v-model="leftDrawerOpen"
      bordered
    >
      <q-list>
        <q-item clickable @click="toggleLeftDrawer">
          <q-item-section></q-item-section>
          <q-item-section avatar side> <q-icon name="close" /> </q-item-section>
        </q-item>

        <q-item
            clickable
            to="/materials"
            @click="toggleLeftDrawer"
        >
          <q-item-section
              avatar
          >
            <q-icon name="o_tune" />
          </q-item-section>
          <q-item-section>
            <q-item-label>            Materials
            </q-item-label>
          </q-item-section>
        </q-item>

        <q-item
            clickable
            to="/info"
            @click="toggleLeftDrawer"
        >
          <q-item-section
              avatar
          >
            <q-icon name="o_info" />
          </q-item-section>
          <q-item-section>
            <q-item-label>            Info
            </q-item-label>
          </q-item-section>
        </q-item>

        <q-separator inset spaced/>
        <q-item> <q-item-section> <q-item-label> External links</q-item-label> </q-item-section> </q-item>
        <q-item clickable tag="a" target="_blank" href="https://github.com/ovidiopr/scattnlay" >
          <q-item-section avatar > <q-icon name="code" /> </q-item-section>
          <q-item-section> <q-item-label>Project at GitHub</q-item-label>
            <q-item-label caption> Open-source software </q-item-label>
          </q-item-section>
        </q-item>
        <q-item clickable tag="a" target="_blank" href="https://github.com/ovidiopr/scattnlay/issues/new?title=[webapp]" >
          <q-item-section avatar > <q-icon name="support" /> </q-item-section>
          <q-item-section> <q-item-label>Support</q-item-label>
            <q-item-label caption> Requires a GitHub account</q-item-label></q-item-section>
        </q-item>
        <q-separator inset spaced/>

        <q-item class="q-mt-auto">
          Last simulation took {{$store.state.simulationSetup.nmieTotalRunTime.toFixed(2)}} s.
        </q-item>
      </q-list>

    </q-drawer>

    <q-page-container>
      <router-view />
    </q-page-container>
  </q-layout>
</template>

<script lang="ts">

import { defineComponent, ref } from 'vue'

export default defineComponent({
  name: 'MainLayout',

  // components: {
  //   EssentialLink
  // },

  setup () {
    const leftDrawerOpen = ref(false)
    return {
      leftDrawerOpen,
      toggleLeftDrawer () {
        leftDrawerOpen.value = !leftDrawerOpen.value
      }
    }
  }
})
</script>
