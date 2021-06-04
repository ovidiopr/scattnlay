<template>
    <div class="field has-addons">
        <p class="control">
            <a class="button is-static input-with-units-title">
                {{ title }}
            </a>
        </p>
        <div v-if="isDisabledLocal">
            <p class="control">
            <b-autocomplete v-model="showValue"
                 class="input-with-units-value" disabled/>

            </p>
        </div>
        <div v-else>
            <p class="control">
                <b-autocomplete v-on:focus="handleFocus($event)" v-on:blur="handleBlur"
                         v-on:keyup.native.enter="focusNext($event)" v-model="showValue"
                                :data="expr_list"
                                :open-on-focus="true"
                                class="input-with-units-value"/>

            </p>
        </div>
        <p class="control">
            <a class="button is-static input-with-units-units">
                 {{ units }}
            </a>
        </p>
    </div>
</template>

<script>
    const math = require('mathjs')
    export default {
        name: "InputWithUnits",
      created() {
        document.addEventListener('focusin', this.focusChanged)
      },
      beforeDestroy() {
        document.removeEventListener('focusin', this.focusChanged)
      },
      methods: {
        focusChanged () {
          if (!this.isFocused) {
            // To manage click on autocomplete item we need to wait
            // autocomplete to update showValue.
            setTimeout(() => {
              this.showValue = this.valueLocal;
            }, 500);
            if (!this.expr_list.includes(this.expr)) this.expr_list.unshift(this.expr);
            if (this.expr_list.length > 5) this.expr_list.pop();
          }
        },
        handleFocus(event){
          if (this.expr == "") this.expr = this.valueLocal;
          this.showValue = this.expr;
          event.target.selectionStart = event.target.selectionEnd;
          this.isFocused = true;
        },
        handleBlur(){
          this.isFocused = false;
          this.focusChanged();
        },
        focusNext(event) {
          const inputs = Array.from(document.querySelectorAll('input[type="text"]'));
          const index = inputs.indexOf(event.target);
          if (index < inputs.length) {
            inputs[index + 1].focus();
            document.activeElement.click();
          }
        }
      },
        watch: {
            valueLocal: {
                handler: function () {
                    this.$emit('newdata',math.evaluate(this.valueLocal));
                }
            },
          showValue: {
              handler: function () {
                try {
                  if (math.evaluate(this.showValue) != this.showValue) this.expr = this.showValue;
                  if (math.evaluate(this.showValue) != math.evaluate(this.expr)) this.expr = this.showValue;
                } catch (e) {return}
              }
            },
          expr: {
            handler: function () {
              this.valueLocal = math.evaluate(this.expr);
            }
          },
            isDisabled: {
                handler: function () {
                    this.isDisabledLocal = this.isDisabled;
                }
            },
            deep: true
        },
        data () {
            return {
                // TODO: Is it OK to modify valueLocal later in <b-input>?
                valueLocal: math.evaluate(this.value),
              expr: this.value,
              expr_list: [this.value],
              showValue: math.evaluate(this.value),
              isFocused: false,
                isDisabledLocal: false,
            }
        },
      mounted() {
        this.$emit('newdata',math.evaluate(this.value));

      },
        props: ['title', 'units', 'value', 'isDisabled']
    }
</script>

<style scoped>
.input-with-units-title {
    width:4rem;
    z-index: 1;
}

.input-with-units-value {
    width:10rem;
}
.input-with-units-units {
    width:3rem;
}
</style>
