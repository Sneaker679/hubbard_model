from hubbard_class import hubbard
from hop_mat import hopping_matrix
import time
import matplotlib.pyplot as plt

def find_num_row_col(N):
    factors = []
    for i in range(1, N + 1):
        if N % i == 0:
            factors.append(i)

    if len(factors) % 2 == 0:
        num_columns = factors[int(len(factors)/2 - 1)]
        num_rows = factors[int(len(factors)/2)]
    else:
        num_columns = factors[int((len(factors)-1)/2)]
        num_rows = num_columns

    return (num_rows,num_columns)

N_values = list(range(1,9))
t_values = []
for N in N_values:
    row, col = find_num_row_col(N)

    start = time.time()
    hub_object = hubbard(N = N, hopping_matrix = hopping_matrix(row,col))
    hub_object.calculate_all_blocks()
    end = time.time()

    exec_time = end - start
    t_values.append(exec_time)
    print(f"N = {N} is done. Time : {exec_time} seconds.")

plt.plot(N_values,t_values)
plt.show()
