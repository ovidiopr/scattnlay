import { RouteRecordRaw } from 'vue-router';

const routes: RouteRecordRaw[] = [
  {
    path: '/',
    component: () => import('layouts/MainLayout.vue'),
    children: [
      { path: '', redirect: 'spectrum' },
      { path: 'spectrum', component: () => import('pages/Spectrum.vue') },
      { path: 'nearfield', component: () => import('pages/Near-field.vue') },
      { path: 'farfield', component: () => import('pages/Far-field.vue') },
      { path: 'info', component: () => import('pages/Info.vue') },
      { path: 'materials', component: () => import('pages/Materials.vue') },
    ],
  },

  // Always leave this as last one,
  // but you can also remove it
  {
    path: '/:catchAll(.*)*',
    component: () => import('pages/Error404.vue'),
  },
];

export default routes;
