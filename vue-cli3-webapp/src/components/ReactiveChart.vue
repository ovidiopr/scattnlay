<template>
    <div :ref="chart.uuid"></div>
</template>
<script>
    export default {
        name: "ReactiveCart",
        props: ["chart"],
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
    }
</script>