<template>
    <section class="show-info-section">
        <button class="button is-primary is-medium"
                @click="isShowInfo = true">
            Show info
        </button>
        <b-modal :active.sync="isShowInfo">
            <div class="message" style="padding: 2rem">
                <span style="text-decoration: underline;text-emphasis: before;font-weight: bold;">Usage:</span>
                <br><br>
                Feel free to use provided software, however, use it at your own risk. We made our best effort to verify
                it is correct, however, we do not provide any warranty.
                <br><br>
                Report bugs, issues, and feature requests at project's <a href="https://github.com/ovidiopr/scattnlay">GitHub</a>

                <br><br>
                If it was useful for your project, please, cite in your related paper the following reference:
                <br><br>
                <article style="margin-left: 1rem">
                    <div style="font-style: italic; padding-bottom: 0.5rem">"Mie calculation of electromagnetic near-field
                        for a multilayered sphere"</div>
                    Konstantin Ladutenko, Umapada Pal, Antonio Rivera, Ovidio Peña-Rodríguez<br>
                    <span style="font-weight: bold">Comp. Phys. Comm., vol. 214, pp. 225–230, 2017</span>
                </article>
            </div>
            <div class="message" style="padding: 2rem">
                <span style="text-decoration: underline;text-emphasis: before;font-weight: bold;">Alternatives:</span>
                <div class="content">
                    <ul>
                        <li><a href="https://nanocomposix.com/pages/mie-theory-calculator">Nanocomposix.</a>  Bulk or core-shell, only dipole and quadrupole contributions.</li>
                        <li><a href="https://de.ifmo.ru/miecalculator/">Ivan Toftul.</a>  Only bulk particles.</li>
                        <li>Lucien Saviot. <a href="https://saviot.cnrs.fr/mie/index.en.html">Bulk</a> and <a href="https://saviot.cnrs.fr/miecoat/index.en.html">core-shell</a> calculators.</li>
                        <li><a href="https://omlc.org/calc/mie_calc.html">Scott Prahl.</a> Angle distribution. </li>
                    </ul>
                </div>
            </div>
        </b-modal>
    </section>
</template>

<script>
    export default {
        name: "ShowInfo",
        data () {
            return {
                isShowInfo: false
            }
        }
    }
</script>

<style scoped>
    .show-info-section {
        padding: 1rem;
    }

</style>
