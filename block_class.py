import numpy as np
from fock_state_class import fock_state

class block:

    def __init__(self, associations : list, t : float, U : float, mu : float, hopping_matrix):
        self._states = self.__find_unique_states(associations)
        self._matrix = self.__calculate_matrix(associations, t, U, mu, hopping_matrix)
        self._energy, self._eigen_vector = self.__diagonalisation()
        self._total_spin = self._states[0]._total_spin
        self._num_electrons = self._states[0]._num_electrons

        for state in self._states:
            state._sign = 1

    def __find_unique_states(self, associations : list):
        states_list = [state for asso in associations for state in asso]
        key = lambda f : f._state 
        states_list = list(set({key(i): i for i in states_list}.values()))
        states_list.sort(key = key)
        #states_list = [copy.deepcopy(state) for state in states_list]
        #for state in states_list:
        #    state._sign = 1
        return states_list

    def __find_permutation_sites(self, association : list):
        if not len(association) == 2:
            raise TypeError("Parameter \"association\" must be of size 2.")

        permutation_indexes = []
        permutation_int = association[0]._state ^ association[1]._state
        for i in range(2 * association[0]._N):

            if (permutation_int & 1) == 1:
                permutation_indexes.append(i)

            permutation_int >>= 1
        
        return (int(permutation_indexes[0] / 2), int(permutation_indexes[1] / 2))

    def __calculate_matrix(self, associations : list, t : float, U : float, mu : float, hopping_matrix):
        N = associations[0][0]._N
        dim_matrix = len(self._states)
        matrix = np.zeros((dim_matrix, dim_matrix))
        states_num = [stateObj._state for stateObj in self._states]

        # t
        for association in associations:
            if not association[0]._state == association[1]._state:
                site1, site2 = self.__find_permutation_sites(association)
                
                if not hopping_matrix[site1,site2] == 0:
                    row = states_num.index(association[0]._state)
                    column = states_num.index(association[1]._state)
                    matrix[row,column] += t * association[1]._sign

        # U
        for index, state in enumerate(states_num):
            new_state = fock_state(N, state)
            for site in range(N):
                for spin in ["+","-"]:
                    new_state.destroy(site,spin)
                    new_state.create(site,spin)

                if not new_state._void:
                    matrix[index,index] += U

                new_state.reset(N, state)
    
        # mu
        for index, state in enumerate(states_num):
            new_state = fock_state(N, state)
            for site in range(N):
                for spin in ["+","-"]:
                    new_state.destroy(site,spin)
                    new_state.create(site,spin)

                    if not new_state._void:
                        matrix[index,index] += -mu

                    new_state.reset(N, state)

        return matrix

    def __diagonalisation(self):
        gs_energies, gs_eigen_vectors = np.linalg.eigh(self._matrix)
        return gs_energies[0], gs_eigen_vectors[:,0]

    def __repr__(self):
        string = f"Lowest Energy : {self._energy}\nEigen Vector : {self._eigen_vector}\n"
        string += f"This block has {self._num_electrons} electron(s) for a total spin of {self._total_spin}.\n"
        string += "States : [" + " ".join(str(state) for state in self._states) + "]\n"
        string += str(self._matrix)
        return string
