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
                <table class="table is-striped is-narrow is-hoverable">
                    <thead>
                    <tr>
                        <th>Use</th>
                        <th>Name</th>
                        <th>Plot</th>
                        <th>File or URL</th>
                    </tr>
                    </thead>
                    <tfoot>
                    <tr>
                        <th>Use</th>
                        <th>Name</th>
                        <th>Plot</th>
                        <th>File or URL</th>
                    </tr>
                    </tfoot>
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
                </table>
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
                isVisible: true,
                chart: {
                    uuid: "materials",
                    traces: [],
                    layout: {
                        // title: 'reactive charts',
                        xaxis: {
                            title: 'Wavelength, mkm'
                        },
                        yaxis: {
                            title: 'refractive index'
                        },
                        width: this.plot_width,
                        height: this.plot_height}
                },
            }
        },
        mounted() {
            let files = ['Au-Johnson-1972.yml','Ag-Johnson-1972.yml'];
            let names = ['Au (Gold) Johnson','Ag (Silver) Johnson'];
            let old_names=[];
            for (const mat of this.materials) old_names.push(mat.name);
            for (let i = 0; i < files.length; i++) {
                if (old_names.includes(names[i])) continue; //Filter out reloads during development
                this.materials.push({
                    fname: files[i],
                    name: names[i],
                    isUsed: true,
                    isPlot: false,
                    isLoaded: false
                });
            }
            this.sortMaterials();
            for (const material of this.materials) {
                this.loadMaterial(material);
            }
        },
        watch: {
          materials: {
              handler: function () {
                  this.chart.traces = [];
                  for (const mat of this.materials) {
                      if (mat.isPlot) {
                          this.chart.traces.push({
                              x: mat.data_nk[0],
                              y: mat.data_nk[1],
                              type: 'scatter',
                              name: mat.name + ' data Re(n)'
                          });
                          this.chart.traces.push({
                              x: mat.data_nk[0],
                              y: mat.data_nk[2],
                              type: 'scatter',
                              name: mat.name + ' data Im(n)'
                          });
                      }

                  }

                  // console.log('update material');
              },
              deep:true
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
            sortMaterials() {
                this.materials.sort((a,b) => (a.name > b.name) ? 1 : ((b.name > a.name) ? -1 : 0));
            },
            deleteMaterial(fname) {
                for (let i = 0; i < this.materials.length; i++) {
                    if (this.materials[i].fname == fname) {
                        this.materials.splice(i);
                    }
                }
            },
            async loadMaterial(material) {
                const data_nk = await this.loadMaterialData(material.fname);
                material.data_nk = data_nk;
            },
            async loadMaterialData(URL){
                const yaml = require('js-yaml');
                // Get document, or throw exception on error
                let Ag_data;

                try {
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
                        return this.transpose(data_num);
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
