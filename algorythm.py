import functions
import time

import warnings
warnings.filterwarnings("ignore")


def algorythm(dir, filename):
    N_0, M_0, IL_0, A_0, B_0 = functions.get_data(f'{dir}\\{filename}.txt')
    print('DATA')
    print(f'N = {N_0}   M = {M_0}')
    print(f'IL = {IL_0}')
    print(f'A = {A_0}')
    print(f'B = {B_0}')

    start_time = time.time()

    A = A_0.copy()
    B = B_0.copy()
    IL = IL_0.copy()
    # for i in range(N_0):
    #     tmp = [0]*N_0
    #     tmp[i] = 1
    #     A.append(tmp)
    #     B.append(0)
    N = N_0
    M = len(B)

    x_candidates = []

    A_sub_list, B_sub_list, IL_sub_list = functions.enum_bases(A, B, N, M, IL)

    for i in range(len(A_sub_list)):
        x_res = functions.enum_subsystem(A_sub_list[i], B_sub_list[i], N, IL_sub_list[i])
        for j in range(len(x_res)):
            if x_res[j] not in x_candidates:
                x_candidates.append(x_res[j])

    x_tmp_list = []
    for x in x_candidates:
        tmp = functions.is_node_natural(x)
        if tmp:
            x_tmp_list.append(tmp)
    x_candidates = x_tmp_list

    x_candidates = functions.filter_outer(A_0, x_candidates, B_0, N_0, M_0, IL_0)

    result = []
    idx = 0
    while idx < len(x_candidates):
        flag = functions.filter_inner(x_candidates[idx], idx, x_candidates)
        if not flag:
            result.append(x_candidates[idx])
        idx += 1

    end_time = time.time()

    print(f'RESULT   {result}')
    print(f'TIME    {end_time - start_time}')
    return end_time - start_time
