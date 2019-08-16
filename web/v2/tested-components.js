Vue.component("reactive-chart", {
    props: ["chart"],
    template: '<div :ref="chart.uuid"></div>',
    mounted() {
        Plotly.newPlot(this.$refs[this.chart.uuid], this.chart.traces,
            this.chart.layout,
            {responsive: true, showSendToCloud: true, displaylogo: false}
        );
    },
    watch: {
        chart: {
            handler: function() {
                Plotly.react(
                    this.$refs[this.chart.uuid],
                    this.chart.traces,
                    this.chart.layout
                );
            },
            deep: true
        }
    }
});


Vue.component('input-with-units',{
    // data: function() {return {value: 303.0}},
    watch: {
        value: {
            handler: function () {
                this.$emit('newdata',this.value);
            }
        },
        deep: true
    },
    props: ['title', 'units', 'value'],
    template: `
                    <div class="field has-addons">
                        <p class="control">
                            <a class="button is-static" style="width:4rem">
                                {{ title }}
                            </a>
                        </p>
                        <p class="control">
                            <b-input v-model="value" type="number" step="any" style="width:6rem"></b-input>
                        </p>
                        <p class="control">
                            <a class="button is-static" style="width:3rem">
                                {{ units }}
                            </a>
                        </p>
                    </div>

    `
})
