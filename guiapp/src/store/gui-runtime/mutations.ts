import { MutationTree } from 'vuex';
import { guiRuntimeStateInterface } from './state';

const mutation: MutationTree<guiRuntimeStateInterface> = {
  someMutation (/* state: guiRuntimeStateInterface */) {
    // your code
  }
};

export default mutation;
