'use strict';
const example = {
    data() {
        return {
            window: {
                width: 0,
                height: 0
            },
            // on change of initial value for __ units __
            // remember to update this.chart.layout.xaxis.title
            units: 'nm',
            simulationSetup: {
                stepWL: 0.5,
                fromWL: 300.0,
                toWL: 1000.0,
                R: 100.0,
                reN: 4.0,
                imN: 0.01,
                total_mode_n:4
            },
            changes: 0,
            isShowInfo: false,
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
                        // initial value should be same as in units
                        title: 'Wavelength, nm'
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
    mounted(){ this.runMie();},
    watch: {
        simulationSetup: {
            handler: function () {
                this.runMie();
            },
            deep: true
        },
        units:{
            handler: function () {
                this.chart.layout.xaxis.title = "Wavelength, "+this.units;
            }
        }
    },
    methods: {
        handleResize() {
            this.window.width = window.innerWidth;
            this.window.height = window.innerHeight;
        },
        runMie(){
            this.changes++;
        }
    },

};


const app = new Vue(example);
app.$mount('#app');