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
        <transition name="slide" class="list">
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
                            <td>{{material.fname}}</td>
                        </tr>
                    </table>
                </div>
            </transition>
        </transition>
    </div>
</template>

<!--https://refractiveindex.info/database/data/main/Au/Johnson.yml-->
<script>
    export default {
        name: "GetMaterials",
        data () {
            return {
                isVisible: true,
            }
        },
        mounted() {
            let files = ['Ag-Johnson-1972.yml','Au-Johnson-1972.yml'];
            let names = ['Ag (Silver) Johnson', 'Au (Gold) Johnson'];
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
            console.log(this.materials)
        },
        props: ['materials']
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

</style>
