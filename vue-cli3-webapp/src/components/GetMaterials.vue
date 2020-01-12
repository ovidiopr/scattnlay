<template>
    <div class="field">
        <div class="field is-horizontal">
            <div class="field-label is-normal">
                <label class="label">Materials</label>
            </div>
            <div class="field-body">
                <b-switch v-model="isVisible">
                    <div v-if="isVisible==true">Hide</div>
                    <div v-else>Show</div>
                    </b-switch>
            </div>
        </div>
        <transition name="slide">
            <div class="list" v-if="isVisible">
                <p style="margin-left: 1rem">Data files were taken from <a href="https://refractiveindex.info" class="is-family-code">refractiveindex.info</a></p>
                <p style="margin-left: 1rem">Request at <a href="https://github.com/ovidiopr/scattnlay/issues/20" class="is-family-code">GitHub</a> to add new materails to the default list.</p>

                <table class="table is-striped is-narrow is-hoverable">
                    <thead>
                    <tr>
                        <th>Use</th>
                        <th>Name</th>
                        <th>Plot</th>
                        <th>File or URL</th>
                    </tr>
                    </thead>
                    <tr v-for="material in materials" v-bind:key="material.name">
                        <td><b-switch v-model="material.isUsed"/>
<!--                                <font-awesome-icon icon="check-circle" class="rh-input has-text-success" v-if="material.isLoaded">-->
<!--                                    <span class="tooltiptext">test</span>-->
<!--                                </font-awesome-icon>-->
                            <span class="rh-input"  v-if="material.isLoaded">
                                <font-awesome-icon icon="check-circle" class="has-text-success"/>
                                <span class="tooltiptext tip-success">Loaded.</span>
                            </span>
                            <span class="rh-input"  v-else>
                                <font-awesome-icon icon="ban" class="has-text-danger"/>
                                <span class="tooltiptext tip-danger">Failed to load.</span>
                            </span>
                        </td>
                        <td>{{material.name}}</td>
                        <td><b-checkbox v-model="material.isPlot"/></td>
                        <td>{{material.fname}}
                            <span class="rh-input">
                                <b-button @click="deleteMaterial(material.fname);" class="is-danger is-small is-outlined">
                                    <font-awesome-icon icon="trash"/>
                                </b-button>
                                <span class="tooltiptext tip-danger">Delete.</span>
                            </span>
                        </td>
                    </tr>
                    <tr>
                        <td><b-button @click="clickAdd()">
                            <font-awesome-icon icon="plus-circle" class="has-text-success"/>
                            Add</b-button>
                        </td>
                        <td><b-input v-model="new_name" placeholder="New name"/></td>
                        <td></td>
                        <td><b-field message="full database record from refractiveindex.info">
                            <b-input v-model="new_fname" placeholder="New URL of *.yml file"/>

                        </b-field>
                        </td>
                    </tr>
                </table>
                <div class="field is-horizontal">
                    <div class="field-label is-normal">
                        <label class="label"> Plot options </label>
                    </div>
                    <div class="field-body">
                        <div class="columns">
                            <div class="column">
                                <b-switch v-model="isPlot_n"> Re(n)</b-switch>
                            </div>
                            <div class="column">
                                <b-switch v-model="isPlot_k"> Im(n)</b-switch>
                            </div>
                            <div class="column">
                                <b-switch v-model="isPlot_spline"> Interpolation</b-switch>
                            </div>
                        </div>
                    </div>
                </div>
                <div class="chart-container">
                    <reactive-chart :chart="chart"/>
                </div>
                <div class="field is-horizontal add-margin">
                    <div class="field-label is-normal">
                        <label class="label"> Materials </label>
                    </div>
                    <div class="field-body">
                        <b-switch v-model="isVisible">
                            <div v-if="isVisible==true">Hide</div>
                            <div v-else>Show</div>
                        </b-switch>
                    </div>
                </div>
            </div>
        </transition>
    </div>
</template>

<!--https://refractiveindex.info/database/data/main/Au/Johnson.yml-->
<script>
    import ReactiveChart from "./ReactiveChart.vue";
    export default {
        name: "GetMaterials",
        components: {
            ReactiveChart,
        },
        data () {
            return {
                isVisible: false,
                chart: {
                    uuid: "materials",
                    traces: [],
                    layout: {
                        margin: {
                            l:10,
                            r:40,
                            b:50,
                            t:10
                        },
                        // title: 'reactive charts',
                        xaxis: {
                            title: 'Wavelength, nm'
                        },
                        yaxis: {
                            title: 'refractive index'
                        },
                        showlegend: true,
                        legend: {
                            orientation:"h",
                            x: -.1,
                            y: 1.2
                        },
                        width: this.plot_width,
                        height: this.plot_height}
                },
                isPlot_n: true,
                isPlot_k: true,
                isPlot_spline: true,
                new_name: '',
                new_fname:'',
            }
        },
        mounted() {
            let files = ['Ag-Johnson-1972.yml',
                'Al-McPeak-2015.yml',
                'Au-McPeak-2015.yml',
                'Cu-McPeak-2015.yml',
                'Si-Green-2008.yml',
                'Ag-McPeak-2015.yml',
                'Au-Johnson-1972.yml',
                'Cu-Johnson-1972.yml',
                'Si-Aspnes-1983.yml'];
            let names = ['Ag (Silver) Johnson',
                'Al (Aluminium) McPeak',
                'Au (Gold) McPeak',
                'Cu (Copper) McPeak',
                'Si (Silicon) Green',
                'Ag (Silver) McPeak',
                'Au (Gold) Johnson',
                'Cu (Copper) Johnson',
                'Si (Silicon) Aspnes'];
            this.addMaterial(files,names);
            for (const mat of this.materials) {
                if (mat.name.includes('Johnson')) mat.isUsed=false;
                if (mat.name.includes('Aspnes')) mat.isUsed=false;
            }
        },
        watch: {
            materials: {
                handler: function () {
                    this.updateChart();
                },
                deep:true
            },
            isPlot_n: {
                handler: function () {
                    this.updateChart();
                }
            },
            isPlot_k: {
                handler: function () {
                    this.updateChart();
                }
            },
            isPlot_spline: {
                handler: function () {
                    this.updateChart();
                }
            },
          plot_width: {
              handler: function () {
                  this.chart.layout.width = this.plot_width;
              }
          },
          plot_height: {
              handler: function () {
                  this.chart.layout.height = this.plot_height;
              }
          },
        },
        methods: {
            transpose(array) {
                return array[0].map((col, i) => array.map(row => row[i]));
            },
            updateChart() {
                this.chart.traces = [];
                for (const mat of this.materials) {
                    if (!mat.isPlot) continue;
                    if (this.isPlot_n) {
                        this.chart.traces.push({
                            x: mat.data_nk[0],
                            y: mat.data_nk[1],
                            mode: 'markers',
                            type: 'scatter',
                            name: mat.name + ' data Re(n)'
                        });}
                    if (this.isPlot_k) {
                        this.chart.traces.push({
                            x: mat.data_nk[0],
                            y: mat.data_nk[2],
                            mode: 'markers',
                            type: 'scatter',
                            name: mat.name + ' data Im(n)'
                        });}
                    // blocks sline interpolation before sline objects are initialized
                    if (!mat.isLoaded) continue;
                    if (!this.isPlot_spline) continue;

                    let x_nk = mat.spline_n.xs;
                    let from_x = parseFloat(x_nk[0]);
                    let to_x = parseFloat(x_nk[x_nk.length-1]);
                    let steps = 1000;
                    let step_x = Math.abs(to_x-from_x)/parseFloat(steps);
                    let spline_x = [];
                    let spline_n =[];
                    let spline_k = [];
                    for (let i = 0; i<steps+1; i++) {
                        let new_x = i * step_x+from_x;
                        if (new_x > to_x) new_x = to_x;
                        spline_x.push(new_x);
                        spline_n.push(mat.spline_n.at(new_x));
                        spline_k.push(mat.spline_k.at(new_x));
                    }
                    if (this.isPlot_n) {
                        this.chart.traces.push({
                            x: spline_x,
                            y: spline_n,
                            type: 'scatter',
                            name: mat.name + ' spline Re(n)'
                        });}
                    if (this.isPlot_k) {
                        this.chart.traces.push({
                            x: spline_x,
                            y: spline_k,
                            mode: 'line',
                            type: 'scatter',
                            name: mat.name + ' spline Im(n)'
                        });}

                }
            },
            sortMaterials() {
                this.materials.sort((a,b) => (a.name > b.name) ? 1 : ((b.name > a.name) ? -1 : 0));
            },
            clickAdd() {
                this.addMaterial([this.new_fname],[this.new_name])
            },
            addMaterial(files, names){
                let old_names=[];
                for (const mat of this.materials) old_names.push(mat.name);
                for (let i = 0; i < files.length; i++) {
                    if (old_names.includes(names[i])) continue; //Filter out reloads during development
                    this.materials.push({
                        fname: files[i],
                        name: names[i],
                        isUsed: true,
                        isPlot: false,
                        isLoaded: false,
                        data_nk: []
                    });
                }
                this.sortMaterials();
                this.materials[0].isPlot= true;
                for (const material of this.materials) {
                    this.loadMaterial(material);
                }
            },
            deleteMaterial(fname) {
                for (let i = 0; i < this.materials.length; i++) {
                    if (this.materials[i].fname == fname) {
                        this.materials.splice(i,1);
                    }
                }
            },
            async loadMaterial(material) {
                const data_nk = await this.loadMaterialData(material.fname);
                material.data_nk = data_nk;
                const Spline = require('cubic-spline');
                const xs = data_nk[0];
                const ys1 = data_nk[1];
                const ys2 = data_nk[2];
                const spline_n = new Spline(xs, ys1);
                const spline_k = new Spline(xs, ys2);
                material.spline_n = spline_n;
                material.spline_k = spline_k;
                material.isLoaded = true;
                // const spline = new Spline(xs, ys);
                // // get Y at arbitrary X
                // console.log(spline.at(1.4));
                // // interpolate a line at a higher resolution
                // for (let i = 0; i < 50; i++) {
                //     console.log(spline.at(i * 0.1));
                // }
            },
            async loadMaterialData(URL){
                const yaml = require('js-yaml');
                // Get document, or throw exception on error
                let Ag_data;

                try {
                    let proxyUrl = 'https://cors-anywhere.herokuapp.com/';
                    if (URL.includes('https://refractiveindex.info')) URL = proxyUrl+URL;
                    let response = await fetch(URL);
                    let Ag_data = await response.text();

                    const doc = await yaml.safeLoad(Ag_data);
                    if (doc.DATA[0].type == "tabulated nk") {
                        let csv = doc.DATA[0].data;
                        let rows = csv.split("\n");

                        let data =  rows.map(function (row) {
                            return row.split(" ");
                        });
                        data.pop();
                        let data_num = data.map(function(elem) {
                            return elem.map(function(elem2) {
                                return parseFloat(elem2);
                            });
                        })
                        let data_columns = this.transpose(data_num);
                        // Convert from default refractiveindex.info mkm to nm

                        for (let i=0; i<data_columns[0].length; i++)
                            data_columns[0][i] *= 1000;
                        return data_columns;
                    }
                } catch (e) {
                    console.log(e);
                }
            }
        },
        props: ['materials', 'plot_width', 'plot_height']
    }
</script>

<style scoped>
    .list{
        transform-origin: top;
        transition: transform .4s ease-in-out;
    }
    .slide-enter, .slide-leave-to{
        transform: scaleY(0);
    }
    .rh-input {
        margin-right: 0px;
        position: relative;
        display: inline-block;
    }
    /* Tooltip text */
    .tip-danger {
        background-color:         #ff3860;
    }
    .tip-success {
        background-color:         #23d160;

    }
    .rh-input .tooltiptext {
        visibility: hidden;
        width: 120px;
        /*background-color: #7957d5;*/
        /*background-color:         #23d160;*/

        color: white;
        text-align: center;
        padding: 5px 0;
        border-radius: 6px;

        /* Position the tooltip text - see examples below! */
        position: absolute;
        z-index: 1;

        bottom: 110%;
        left: 0%;
    }

    /* Show the tooltip text when you mouse over the tooltip container */
    .rh-input:hover .tooltiptext {
        visibility: visible;
    }

    .add-margin {
        padding-bottom: 1rem;
    }
</style>
