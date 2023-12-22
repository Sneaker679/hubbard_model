import numpy as np
import copy

class fock_state: 

    def __init__(self, N : int, initial_state : int = 0):

        self._state = initial_state
        self._num_electrons = self.__ones(initial_state)
        self._N = N
        self._sign = 1
        self._total_spin = self.__calculate_total_spin(initial_state)
        self._void = False
        
    def __calculate_total_spin(self, state):
        spin = 0

        for bit in range(2 * self._N):
            if bit % 2 == 0:
                spin -= (state >> bit) & 1
            else:
                spin += (state >> bit) & 1

        return spin

    def __ones(self, num : int):
        return bin(num).replace("0b", "").count("1")

    def create(self, site : int, spin : str):
        
        if not self._void:

            index = site*2 + (spin == '+')
            create_op = 1 << index
            
            if self._state | create_op == self._state + create_op:
                self._state += create_op
                truncated_state = self._state & (( 2**(2 * self._N) - 1) << index + 1)

                self._sign *= (-1)**self.__ones(truncated_state)
                self._total_spin += -1 if index % 2 == 0 else 1
                self._num_electrons += 1

                # successful creation
            else:
                self.__init__(self._N, 0)
                self._void = True # unsuccesful creation | state would be destroyed - void
        
    def destroy(self, site : int, spin : str):

        if not self._void:

            #check parameter validity

            index = site*2 + (spin == '+')
            destroy_op = 1 << index
            
            if self._state ^ destroy_op == self._state - destroy_op:
                self._state -= destroy_op
                truncated_state = self._state & (( 2**(2 * self._N) - 1) << index + 1)

                self._sign  *= (-1)**self.__ones(truncated_state)
                self._total_spin -= -1 if index % 2 == 0 else 1
                self._num_electrons -= 1

                # successful annihilation
            else:
                self.__init__(self._N, 0)
                self._void = True # unsuccesful annihilation | state would be destroyed - void
    
    def reset(self, N : int, initial_state : int = 0):
        self.__init__( N, initial_state)

    def copy_object(self):
        return copy.deepcopy(self)

    def __repr__(self):
        if not self._void:
            return str(self._sign * self._state)
        else:
            return "void"
