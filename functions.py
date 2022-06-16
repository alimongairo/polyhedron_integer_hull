import numpy as np
from itertools import combinations, product
from scipy.optimize import linprog


def get_data(file):
    with open(file, mode='r') as f:
        line = f.readline().split(sep=' ')
        n = int(line[0])
        m = int(line[1])
        f.readline()
        ineq_list = f.readline().split(sep=' ')
        ineq_list[-1] = ineq_list[-1].strip('\n')
        if ineq_list[-1] == '':
            ineq_list.pop(-1)
        f.readline()
        a = []
        for _ in range(m):
            tmp = []
            line = f.readline().split(sep=' ')
            if line[-1] == '\n':
                line.pop(-1)
            for i in line:
                tmp.append(int(i))
            a.append(tmp)
        f.readline()
        b = []
        line = f.readline().split(sep=' ')
        if line[-1] == '\n':
            line.pop(-1)
        for i in line:
            b.append(int(i))
    return n,m,ineq_list,a,b


def system_solve(a, b):     # solver for system Ax=B
    res = np.linalg.solve(a, b)
    return res


def enum_bases(A, B, N, M, IL):     # creation of NxN bases
    row_idx = [i for i in range(M)]  # all row indexes
    row_sets_list = [i for i in combinations(row_idx, N)]  # all subsets I, where |I| = N

    a_sub_list = []  # set of sub-matrices of A for each I
    b_sub_list = []  # set of sub-matrices of B for each A
    il_sub_list = []  # set of sub-matrices of inequality type for each A
    for row_set in row_sets_list:  # for each I
        a_sub = []  # sub-matrix of A
        b_sub = []  # sub-matrix of B
        il_sub = []  # sub-matrix of il
        for s in row_set:
            a_sub.append(A[s])
            b_sub.append(B[s])
            il_sub.append(IL[s])
        if np.linalg.det(a_sub) != 0:  # filter non-singular bases
            a_sub_list.append(a_sub)
            b_sub_list.append(b_sub)
            il_sub_list.append(il_sub)
    return a_sub_list, b_sub_list, il_sub_list   # all subsystems of A and B (dim NxN)


def is_node_natural(x, tolerance=0.001):    # filter for not negative int vertices
    res = []
    for c in x:
        if abs(c-round(c)) > tolerance: # or c < 0:
            return False
        else:
            res.append(round(c))
    return res


def enum_subsystem(a_sub, b_sub, N, il_sub):
    x_list = []
    det_i = abs(np.linalg.det(a_sub))   # det I = |det A_i|
    xi_set = [i for i in product([j for j in range(round(det_i))], repeat=N)]  # set of all xi vectors
    for xi in xi_set:  # solving system Ax=B-xi for each xi
        b_new = [0] * N
        for i in range(N):
            if il_sub[i] == 'u':
                b_new[i] = b_sub[i] - xi[i]
            else:
                b_new[i] = b_sub[i] + xi[i]
        x = system_solve(a_sub, b_new).tolist()
        x_list.append(x)
    return x_list   # list of all x-candidates


def filter_outer(a, x_list, b, n, m, il):   # check if Ax<=B
    result = []
    for x in x_list:
        ax = np.array(a).dot(np.array(x))
        flag = 0
        for i in range(m):
            if ax[i] > b[i] and il[i] == 'u':
                flag = 1
                break
            elif ax[i] < b[i] and il[i] == 'b':
                flag = 1
                break
        if flag == 0:
            result.append(x)
    return result   # desired and inner x-candidates


def filter_inner(v, v_idx, v_full_list):    # check if v_i is not in conv M\{v_i}
    v_list = v_full_list.copy()
    v_list.pop(v_idx)
    v_matrix = np.transpose(np.array(v_list))

    v_matrix = list(v_matrix)
    v_matrix.append(np.array([1.0 for _ in range(len(v_list))]))
    v_matrix = np.array(v_matrix)
    # np.append(v_matrix, [[1.0 for _ in range(len(v_list))]], axis=0)

    xi_matrix = np.transpose(np.array(v))
    xi_matrix = list(xi_matrix)
    xi_matrix.append(1.0)
    xi_matrix = np.array(xi_matrix)
    # np.append(xci_matrix, [1.0], axis=0)

    bounds = [(0,None) for _ in range(len(v_matrix[0]))]
    alpha_list = linprog([1]*len(v_matrix[0]), A_eq=v_matrix, b_eq=xi_matrix, bounds=bounds)
    return alpha_list.success
