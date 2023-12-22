import numpy as np
import time
import copy

from fock_state_class import fock_state
from block_class import block
from hop_mat import hopping_matrix

class hubbard:

    def __init__(self, N : int = 2, t = -1, U = 4, mu = 2, hopping_matrix = np.matrix([[0,1],[1,0]])):

        self._N = N
        self._t = t
        self._U = U
        self._mu = mu
        self._hopping_matrix = hopping_matrix

        self._blocks = self.__calculate_blocks()
        self._gs_blocks = self.__find_gs()

    def __calculate_associations(self, initial_state : fock_state):
        block_states = [initial_state]
        associations = []
        
        for state in block_states:
            for i, j in np.ndindex( (self._N, self._N,) ):
                for spin in ['+','-']:
                    if not i == j and not self._hopping_matrix[i,j] == 0:

                        new_state = state.copy_object()
                        new_state.destroy(i,spin)
                        new_state.create(j,spin)

                        if not new_state._void:
                            associations.append((state,new_state))

                            if not new_state._state in [state._state for state in block_states]:
                                new_block_state = new_state.copy_object()
                                new_block_state._sign = 1
                                block_states.append(new_block_state)
        
        if not associations:
            associations.append((initial_state,initial_state))
            
        return associations

    def __calculate_blocks(self):
        bank = [state for state in range(2**(2*self._N))]
        blocks = [] * (self._N * 2 + 1) 
        while bank:
            initial_state = fock_state(self._N, bank[0])
            associations = self.__calculate_associations(initial_state)
            new_block = block(associations, self._t, self._U, self._mu, self._hopping_matrix) 
            blocks.append(new_block)

            to_remove = [state._state for state in blocks[-1]._states]
            for state in to_remove:
                bank.remove(state)
        
        blocks_Nsorted = []
        for n in range(2 * self._N + 1):
            blocks_Nsorted.append( [block for block in blocks if block._num_electrons == n] )

        key = lambda block : block._total_spin
        for group_index in range(2 * self._N):
            blocks_Nsorted[group_index].sort(key = key)

        return blocks_Nsorted
            
    def __find_gs(self):
        gs_energy = self._blocks[0][0]._energy
        gs_blocks = [self._blocks[0][0]]

        for NIndex, block_group in enumerate(self._blocks):
            for SpinIndex, block in enumerate(block_group):

                if np.abs(block._energy - gs_energy) < 0.001:
                    gs_blocks.append(block)

                elif block._energy < gs_energy:
                    gs_energy = block._energy
                    gs_blocks.clear()
                    gs_blocks.append(block)
        
        if not len(gs_blocks) == 1:
            gs_blocks.sort(key = lambda b: b._total_spin) 

        return gs_blocks

    def get_block(self, num_electrons : int, spin : int):
        try:
            block_spin_index = None
            for index, block in enumerate(self._blocks[num_electrons]):
                if block._total_spin == spin:
                    block_spin_index = index
                    break

            return self._blocks[num_electrons][block_spin_index]

        except (IndexError, TypeError):
            return None
            
    def get_gs_block(self, spin : int):
        try:
            block_spin_index = None
            for index, block in enumerate(self._gs_blocks):
                if block._total_spin == spin:
                    block_spin_index = index
                    break

            return self._gs_blocks[block_spin_index]

        except (IndexError, TypeError):
            return None

    def get_biggest_spin_gs_block(self):
        return self._gs_blocks[-1]

    def get_smallest_spin_gs_block(self):
        return self._gs_blocks[0]

    def __repr__(self):
        string = ''
        i = 0

        for block_group in self._blocks:
            for block in block_group:
                string += f">>>BLOCK #{i}<<<\n"
                string += str(block) + "\n\n"
                i += 1

        return string

if __name__ == "__main__":
    start = time.time()
    result = hubbard(N = 2, t = -1, U = 4, mu = 2, hopping_matrix = hopping_matrix(1,2))
    end = time.time()
    print(f"Time : {end - start} seconds.")

    try:
        import resource
        print(f"Max Memory Usage : {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1000} megabytes.")
    except ModuleNotFoundError:
        pass

