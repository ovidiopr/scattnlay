<template>
  <q-layout view="hHh lpr lFf">
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
        <q-tabs align="left">
          <q-route-tab to="/" label="Spectrum" />
          <q-route-tab to="/page2" label="Near-field" />
          <q-route-tab to="/page3" label="Far-field" />
        </q-tabs>
      </q-toolbar>
    </q-header>

    <q-drawer
      v-model="leftDrawerOpen"
      show-if-above
      bordered
    >
      <q-list>

        <q-item
            clickable
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
        <q-item clickable tag="a" target="_blank" :href="`https://github.com/ovidiopr/scattnlay`" >
          <q-item-section avatar > <q-icon name="code" /> </q-item-section>
          <q-item-section> <q-item-label>Project at GitHub</q-item-label>
            <q-item-label caption> Open-source software </q-item-label>
          </q-item-section>
        </q-item>
        <q-item clickable tag="a" target="_blank" :href="`https://github.com/ovidiopr/scattnlay/issues/new?title=[webapp]`" >
          <q-item-section avatar > <q-icon name="support" /> </q-item-section>
          <q-item-section> <q-item-label>Support</q-item-label>
            <q-item-label caption> Requires a GitHub account</q-item-label></q-item-section>
        </q-item>
<!--        <q-item-label-->
<!--          header-->
<!--        >-->
<!--          Essential Links-->
<!--        </q-item-label>-->

<!--        <EssentialLink-->
<!--          v-for="link in essentialLinks"-->
<!--          :key="link.title"-->
<!--          v-bind="link"-->
<!--        />-->
      </q-list>
    </q-drawer>

    <q-page-container>
      <router-view />
    </q-page-container>
  </q-layout>
</template>

<script lang="ts">
// import EssentialLink from 'components/EssentialLink.vue'

// const linksList = [
//   {
//     title: 'Docs',
//     caption: 'quasar.dev',
//     icon: 'school',
//     link: 'https://quasar.dev'
//   },
//   {
//     title: 'Github',
//     caption: 'github.com/quasarframework',
//     icon: 'code',
//     link: 'https://github.com/quasarframework'
//   },
//   {
//     title: 'Discord Chat Channel',
//     caption: 'chat.quasar.dev',
//     icon: 'chat',
//     link: 'https://chat.quasar.dev'
//   },
//   {
//     title: 'Forum',
//     caption: 'forum.quasar.dev',
//     icon: 'record_voice_over',
//     link: 'https://forum.quasar.dev'
//   },
//   {
//     title: 'Twitter',
//     caption: '@quasarframework',
//     icon: 'rss_feed',
//     link: 'https://twitter.quasar.dev'
//   },
//   {
//     title: 'Facebook',
//     caption: '@QuasarFramework',
//     icon: 'public',
//     link: 'https://facebook.quasar.dev'
//   },
//   {
//     title: 'Quasar Awesome',
//     caption: 'Community Quasar projects',
//     icon: 'favorite',
//     link: 'https://awesome.quasar.dev'
//   }
// ];

import { defineComponent, ref } from 'vue'

export default defineComponent({
  name: 'MainLayout',

  // components: {
  //   EssentialLink
  // },

  setup () {
    const leftDrawerOpen = ref(false)

    return {
      // essentialLinks: linksList,
      leftDrawerOpen,
      toggleLeftDrawer () {
        leftDrawerOpen.value = !leftDrawerOpen.value
      }
    }
  }
})
</script>
