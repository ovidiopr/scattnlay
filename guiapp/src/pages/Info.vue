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
      <a href="https://github.com/ovidiopr/scattnlay">GitHub</a>
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
          <li><a href="https://nanocomposix.com/pages/mie-theory-calculator">Bulk or core-shell,</a>   only dipole and quadrupole contributions by Nanocomposix. </li>
          <li><a href="https://de.ifmo.ru/miecalculator/">Bulk particles</a>   by Ivan Toftul.</li>
          <li><a href="https://saviot.cnrs.fr/mie/index.en.html">Bulk</a> and <a href="https://saviot.cnrs.fr/miecoat/index.en.html">core-shell</a> calculators by Lucien Saviot.</li>
          <li><a href="https://omlc.org/calc/mie_calc.html">Angle distribution</a>  by Scott Prahl.</li>
        </ul>
      </q-card-section>
    </q-card>
  </q-page>
</template>

<script>
export default {
  // name: 'PageName',
}
</script>
<style lang="scss">
 a {
   margin: 0.2rem;

 }
</style>
