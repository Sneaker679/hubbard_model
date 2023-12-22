import numpy as np

def hopping_matrix(num_rows : int, num_columns : int):
    lattice_fock = []
    for row in range(num_rows):
        for column in range(num_columns):
            lattice_fock.append(np.matrix([row,column]))

    N = num_rows * num_columns
    hopping_matrix = np.zeros((N,N))
    for index1,site1 in enumerate(lattice_fock):
        for index2,site2 in enumerate(lattice_fock):
            
            hop = np.subtract(site1,site2)
            if (np.array_equal(hop,np.matrix([1,0]))
                or np.array_equal(hop,np.matrix([-1,0]))
                or np.array_equal(hop,np.matrix([0,1]))
                or np.array_equal(hop,np.matrix([0,-1]))):

                hopping_matrix[index1,index2] = 1
    return hopping_matrix
