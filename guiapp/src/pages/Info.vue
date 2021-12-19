<template>
  <q-page padding>
    <!-- content -->
    <q-card class="q-mt-sm q-mx-sm">
      <q-card-section>
        <div class="text-h6">Usage:</div>
<!--        <div class="text-subtitle2">by John Doe</div>-->
      </q-card-section>
      <q-card-section>
      Feel free to use provided software, however, use it at your own risk. We made our best effort to verify
      it is correct, however, we do not provide any warranty.
      <br><br> Report bugs, issues, and feature requests at project's
      <a href="https://github.com/ovidiopr/scattnlay">GitHub<q-icon name="launch" style="padding-left: 0.1rem; padding-bottom: 1rem" /></a>
      (preferred) or mail to <a href="mailto:k.ladutenko@metalab.ifmo.ru?subject=Mie calculator webapp (Scattnlay)">k.ladutenko@metalab.ifmo.ru</a>.
      <br><br>
      If it was useful for your project, please, cite the following reference:
      <br><br>
      <article style="margin-left: 1rem">
        <div style="font-style: italic; padding-bottom: 0.5rem">"Mie calculation of electromagnetic near-field
          for a multilayered sphere"</div>
        Konstantin Ladutenko, Umapada Pal, Antonio Rivera, Ovidio Peña-Rodríguez<br>
        <span style="font-weight: bold">Comp. Phys. Comm., vol. 214, pp. 225–230, 2017</span>
      </article>
      <br>
      and the related URL:

      <article style="margin: 1rem">
        <div style="font-style: italic">https://physics.itmo.ru/mie/</div>
      </article>
      </q-card-section>
    </q-card>
    <q-card class="q-mt-md  q-mx-sm">
      <q-card-section>
        <div class="text-h6">Alternatives:</div>
      </q-card-section>
      <q-card-section>
        <ul>
          <li><a href="https://nanocomposix.com/pages/mie-theory-calculator">Bulk or core-shell<q-icon name="launch" style="padding-left: 0.1rem; padding-bottom: 1rem" />,</a>   only dipole and quadrupole contributions by Nanocomposix. </li>
          <li><a href="https://de.ifmo.ru/miecalculator/">Bulk particles<q-icon name="launch" style="padding-left: 0.1rem; padding-bottom: 1rem" /></a>   by Ivan Toftul.</li>
          <li><a href="https://saviot.cnrs.fr/mie/index.en.html">Bulk<q-icon name="launch" style="padding-left: 0.1rem; padding-bottom: 1rem" /></a> and <a href="https://saviot.cnrs.fr/miecoat/index.en.html">core-shell<q-icon name="launch" style="padding-left: 0.1rem; padding-bottom: 1rem" /></a> calculators by Lucien Saviot.</li>
          <li><a href="https://omlc.org/calc/mie_calc.html">Angle distribution<q-icon name="launch" style="padding-left: 0.1rem; padding-bottom: 1rem" /></a>  by Scott Prahl.</li>
        </ul>
      </q-card-section>
    </q-card>
    <!-- place QPageScroller at end of page -->
    <q-page-scroller position="bottom-right" :scroll-offset="150" :offset="[18, 18]">
      <q-btn size="xs" fab icon="keyboard_arrow_up" color="primary" />
    </q-page-scroller>
  </q-page>
</template>

<script>
export default {
  name: 'InfoPage',
}
</script>
