const example = {
    data() {
        return {
            window: {
                width: 0,
                height: 0
            },
            textStyle:"underline",
            styleObject: {
                color: 'red',
                fontSize: '13px'
            },
            units: 'nm',
            stepWL: 0.5,
            fromWL: 300.0,
            toWL:1000.0,
            R: 100.0,
            reN: 4.0,
            imN: 0.01,
            isShowInfo: false,
            total_mode_n:4,
            chart: {
                uuid: "123",
                traces: [
                    {
                        y: [],
                        line: {
                            color: "#5e9e7e",
                            width: 4,
                            shape: "line"
                        }
                    }
                ],
                layout: {
                    title:'reactive charts',
                    xaxis: {
                        title: 'xaxis title'
                    },
                    yaxis: {
                        title: 'yaxis title'
                    }
                }
            }
        };

    },
    created() {
        window.addEventListener('resize', this.handleResize)
        this.handleResize();
    },
    destroyed() {
        window.removeEventListener('resize', this.handleResize)
    },
    methods: {
        handleResize() {
            this.window.width = window.innerWidth;
            this.window.height = window.innerHeight;
        }
    },
    watch: {
        window: function () {
            this.styleObject.color ='blue';
        }
    }

};


Vue.component("reactive-chart", {
    props: ["chart"],
    template: '<div :ref="chart.uuid"></div>',
    mounted() {
        Plotly.newPlot(this.$refs[this.chart.uuid], this.chart.traces, this.chart.layout, {responsive: true});
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

const app = new Vue(example);

app.$mount('#app');