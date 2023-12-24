import numpy as np
import bisect
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

        self._bank = [state for state in range(2**(2*N))]

        # The index for this list corresponds to the number of electrons.
        # For example, the first element corresponds to all blocks that have one electron.
        self._blocks = []
        for n in range(2 * N + 1):
            self._blocks.append([])

        # Calculated only if self.calculate_all_blocks() is ran.
        self._gs_blocks = None
        self._gs_energy = None

    # This function is ran by another function called self.calculate_all_blocks().
    # It proceeds to scan all the calculated blocks, and finds the blocks with the lowest energy.
    def __find_gs(self):
        self._gs_energy = self._blocks[0][0]._energy
        self._gs_blocks = [self._blocks[0][0]]

        # Finding gs_blocks
        for NIndex, block_group in enumerate(self._blocks):
            for SpinIndex, block in enumerate(block_group):

                if np.abs(block._energy - self._gs_energy) < 0.001:
                    self._gs_blocks.append(block)

                elif block._energy < self._gs_energy:
                    self._gs_energy = block._energy
                    self._gs_blocks.clear()
                    self._gs_blocks.append(block)
        
        # Sorting by spin value
        if not len(self._gs_blocks) == 1:
            self._gs_blocks.sort(key = lambda b: b._total_spin) 

    # Calculates all the states of a block, starting from an initial_state, and returns
    # all associations made. These will correspond the indexes of "t" in the hubbard hamiltonian.
    def __calculate_associations(self, initial_state : fock_state):
        block_states = [initial_state]
        associations = []
        
        # Calculation is derived from the fermi-hubbard hamiltonian formula.
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
        
        # To prevent the code from crashing, the block is of size 1x1, then
        # we return an association with itself, for example, (0,0) or (1,1)...
        if not associations:
            associations.append((initial_state,initial_state))
            
        return associations

    # Calculates a block by giving it the initial state it needs to start calculating.
    # Updates the object's attributes accordingly.
    def __calculate_block_from_init_state(self, initial_state : fock_state):
        associations = self.__calculate_associations(initial_state)
        new_block = block(associations, self._t, self._U, self._mu, self._hopping_matrix) 
        bisect.insort(
            self._blocks[new_block._num_electrons],
            new_block,
            key = lambda block : block._total_spin)

        to_remove = [state._state for state in new_block._states]
        for state in to_remove:
            self._bank.remove(state)

    # Public method allowing user to calculate one specific bloc by feeding it the number of electrons
    # and the total spin.
    def calculate_block(self, N_electrons, spin):
        # Checking if paramters make for a possible block
        if not (N_electrons % 2) == (spin % 2):
            print("Impossible block. Skipping instruction.")
            return
            
        # Checks if block is already calculated
        for spin_index, block in enumerate(self._blocks[N_electrons]):
            if block._total_spin == spin:
                return

        # Finds the first state that has these paramters as attributes
        found = False
        for num in self._bank:
            initial_state = fock_state(self._N, num)
            if initial_state._num_electrons == N_electrons and initial_state._total_spin == spin:
                found = True
                break

        # If this state is not found, then the block doesn't exist.
        # If found, then calculates the block with the first state.
        if not found:
            print("Impossible block. Skipping instruction.")
            return

        else: 
            self.__calculate_block_from_init_state(initial_state)

    # Calculates all blocks and finds the ground state by calling __find_gs()
    def calculate_all_blocks(self):
        while self._bank:
            initial_state = fock_state(self._N, self._bank[0])
            self.__calculate_block_from_init_state(initial_state)
        self.__find_gs()
            
    # Getter for one block
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
            
    # Getter for one of the ground block(s)
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

        if not self._gs_energy is None:
            string += f"Ground state energy : {self._gs_energy}.\nThese blocks are part of the GS:\n"
            for block in self._gs_blocks:
                string += f"\tN:{block._num_electrons}/S:{block._total_spin}\n"

        return string
