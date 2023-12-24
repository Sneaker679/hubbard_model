from time import time
from hubbard_class import hubbard
from hop_mat import hopping_matrix


start = time()

hub_object = hubbard(N = 2, t = -1, U = 4, mu = 2, hopping_matrix = hopping_matrix(2,1))
#hub_object.calculate_block(2,0)
hub_object.calculate_all_blocks()

end = time()

#print(hub_object)
print(f"Time : {end - start} seconds.")

try:
    import resource
    print(f"Max Memory Usage : {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1000} megabytes.")
except ModuleNotFoundError:
    pass
