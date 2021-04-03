import sympy
import random
import time
import copy
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pylab
pylab.rcParams['figure.figsize'] = (7, 5)
import scipy.optimize as opt

from logger import init_logger
logger = init_logger()

POOLS = 3
Placeholder = (1 << 6)
MAX_POWER = (1 << 6)
Threshold = 1/(1 << 18)

def print3d(x, y, z, label_x, label_y, label_z, com = ''):
    fig = plt.figure()
    title = 'Setting: ' + com
    fig.suptitle(title, fontsize=40)
    ax = fig.add_subplot(111, projection='3d')
    for yy in y:
        ax.scatter(x, yy, z)
    ax.set_xlabel(label_x, fontsize=40)
    ax.set_ylabel(label_y, fontsize=40)
    ax.set_zlabel(label_z, fontsize=40)
    ax.legend()
    timestamp = time.strftime("%Y_%m_%d_%H_%M_%S", time.localtime())
    plt.savefig(f"fig_{timestamp}.png")
    return

def print2d(x, y, label_x, label_y, com = ''):
    color=['red','yellow','green','blue','black']
    plt.figure()
    title = 'Setting: ' + com
    plt.title(title)
    if y.shape[0] == POOLS:
        for i in range(POOLS):
            plt.plot(x, y[i], c = color[i % len(color)], label = f'pool {i}')
        plt.legend()
        timestamp = time.strftime("%Y_%m_%d_%H_%M_%S", time.localtime()) + '_revenue'
    else:
        plt.plot(x, y)
        timestamp = time.strftime("%Y_%m_%d_%H_%M_%S", time.localtime()) + '_ppoa'
    plt.xlabel(label_x)
    plt.ylabel(label_y)
    plt.savefig(f"fig_{timestamp}.png")
    return

def get_x(x, i, j):
    return x[i*POOLS + j] if i != j else 0

def print_r(res):
    logger.info('#'*Placeholder)
    for i in range(POOLS):
        logger.info(f"The average revenue of pool {i} is {sympy.simplify(res[r[i]])}")
    logger.info('#'*Placeholder)
    logger.info("")
    return

def argmax(func, upper_bound):
    fun = sympy.lambdify([y], - func / Threshold)
    cons = ({'type': 'ineq', 'fun': lambda z: upper_bound - sum(z)})
    bnds = tuple([(0, None) for i in range(POOLS-1)])
    stas = tuple([0 for i in range(POOLS-1)])
    res = opt.minimize(fun, stas, method='SLSQP', bounds=bnds, constraints=cons)
    # logger.info('#'*Placeholder)
    # logger.info(res)
    # logger.info('#'*Placeholder)
    # logger.info("")
    return res

def take_actions(res, value_m, A, new_A):
    res_copy = copy.deepcopy(res)
    for i in range(POOLS):
        for j in range(POOLS):
            if i == j:
                continue
            for k in range(POOLS):
                if j == k:
                    continue
                res_copy[r[i]] = res_copy[r[i]].subs(get_x(x, j, k), A[j][k])
    # print_r(res_copy)
    max_revenue = [0] * POOLS
    for i in range(POOLS):
        cnt = 0
        for j in range(POOLS):
            if i == j:
                continue
            res_copy[r[i]] = res_copy[r[i]].subs(get_x(x, i, j), y[cnt])
            # logger.info(f"The average revenue of pool {i} is {sympy.simplify(res_copy[r[i]])}")
            cnt += 1
        value_f = argmax(res_copy[r[i]], value_m[i])
        max_revenue[i] = - value_f.fun * Threshold
        cnt = 0
        for j in range(POOLS):
            if i == j:
                continue
            new_A[i][j] = value_f.x[cnt]
            cnt += 1
    # print_r(res_copy)
    return max_revenue

def calculate_equilibrium(res, value_m, value_t, value_p, show_converge = False):
    res_copy = copy.deepcopy(res)
    for i in range(POOLS):
        res_copy[r[i]] = res_copy[r[i]].subs(t, value_t)
        res_copy[r[i]] = res_copy[r[i]].subs(p, value_p)
        for j in range(POOLS):
            res_copy[r[i]] = res_copy[r[i]].subs(m[j], value_m[j])
        # res_copy[r[i]] = sympy.simplify(res_copy[r[i]])
    # print_r(res_copy)

    A = np.zeros((POOLS, POOLS))
    start_time = time.time()
    step = 1
    revenue_converge = [[0]*POOLS]
    for i in range(POOLS):
        revenue_converge[0][i] = float(res_copy[r[i]].subs([(xx, 0) for xx in x]))
    ppoa_converge = [1]
    while time.time() - start_time < 12:
        new_A = np.zeros((POOLS, POOLS))
        max_revenue = take_actions(res_copy, value_m, A, new_A)
        step += 1
        revenue_converge.append(max_revenue)
        ppoa_converge.append(sum(value_m) / (sum(value_m) - np.sum(new_A)))
        if (np.abs(new_A-A) < Threshold * MAX_POWER).all():
            break
        else:
            A = new_A

    if show_converge:
        com = ""
        for i in range(POOLS):
            com += f"m{i} = {value_m[i]}   "
        com += f"t = {value_t}   p = {value_p}"
        print2d(list(range(step)), np.transpose(revenue_converge), 'step', 'Average Revenue', com = com)
        print2d(list(range(step)), np.array(ppoa_converge), 'step', 'PPoA', com = com)

    logger.info('#'*Placeholder)
    logger.info(f"The mining powers are {value_m}")
    logger.info(f"The left power is {value_t}, the betrayal parameter is {value_p}")
    logger.info('#'*Placeholder)
    logger.info(f'The strategy matrix is:\n{A}')
    logger.info('#'*Placeholder)
    ppoa = sum(value_m) / (sum(value_m) - np.sum(A))
    logger.info(f'PPoA is: {ppoa}')
    logger.info('#'*Placeholder)
    logger.info("")
    return ppoa

def simulate_x_m1_y_m2(res):
    axis_x = []
    axis_y = []
    ppoa_set = []
    # value_t = random.randint(1, MAX_POWER)
    value_t = 0
    # value_p = random.random()
    value_p = 0
    for i in range(1, MAX_POWER):
        for j in range(1, MAX_POWER):
            value_m = [MAX_POWER>>1, i, j]
            ppoa = calculate_equilibrium(res, value_m, value_t, value_p)
            axis_x.append(i)
            axis_y.append(j)
            ppoa_set.append(ppoa)
    logger.info('#'*Placeholder)
    logger.info(f"Fixing m0, t, p, the x axis is m1 and the y axis is m2")
    logger.info(ppoa_set)
    logger.info('#'*Placeholder)
    logger.info("")
    print3d(axis_x, axis_y, ppoa_set, 'm1', 'm2', 'PPoA', com = f'm0 = {value_m[0]} t = {value_t} p = {value_p}')
    return

def simulate_x_m0m1_y_m2(res):
    axis_x = []
    axis_y = []
    ppoa_set = []
    # value_t = random.randint(1, MAX_POWER)
    value_t = MAX_POWER>>1
    # value_p = random.random()
    value_p = 0
    for i in range(1, MAX_POWER):
        for j in range(1, MAX_POWER):
            value_m = [i, i, j]
            ppoa = calculate_equilibrium(res, value_m, value_t, value_p)
            axis_x.append(i)
            axis_y.append(j)
            ppoa_set.append(ppoa)
    logger.info('#'*Placeholder)
    logger.info(f"Fixing t, p, the x axis is m0 = m1 and the y axis is m2")
    logger.info(ppoa_set)
    logger.info('#'*Placeholder)
    logger.info("")
    print3d(axis_x, axis_y, ppoa_set, 'm0 = m1', 'm2', 'PPoA', com = f't = {value_t} p = {value_p}')
    return

def simulate_x_m0m1m2_y_t(res):
    axis_x = []
    axis_y = []
    ppoa_set = []
    # value_p = random.random()
    value_p = 0
    for i in range(1, MAX_POWER):
        for j in range(1, 2 * MAX_POWER):
            value_m = [i, i, i]
            value_t = j
            ppoa = calculate_equilibrium(res, value_m, value_t, value_p)
            axis_x.append(i)
            axis_y.append(j)
            ppoa_set.append(ppoa)
    logger.info('#'*Placeholder)
    logger.info(f"Fixing p, the x axis is m0 = m1 = m2 and the y axis is t")
    logger.info(ppoa_set)
    logger.info('#'*Placeholder)
    logger.info("")
    print3d(axis_x, axis_y, ppoa_set, 'm0 = m1 = m2', 't', 'PPoA', com = f'p = {value_p}')
    return

def simulate_x_m2_y_t(res):
    axis_x = []
    axis_y = []
    ppoa_set = []
    # value_p = random.random()
    value_p = 0
    for i in range(1, MAX_POWER):
        for j in range(1, 2 * MAX_POWER):
            value_m = [MAX_POWER>>1, MAX_POWER>>1, i]
            value_t = j
            ppoa = calculate_equilibrium(res, value_m, value_t, value_p)
            axis_x.append(i)
            axis_y.append(j)
            ppoa_set.append(ppoa)
    logger.info('#'*Placeholder)
    logger.info(f"Fixing m0, m1, p, the x axis is m2 and the y axis is t")
    logger.info(ppoa_set)
    logger.info('#'*Placeholder)
    logger.info("")
    print3d(axis_x, axis_y, ppoa_set, 'm2', 't', 'PPoA', com = f'm0 = {value_m[0]} m1 = {value_m[1]} p = {value_p}')
    return

def simulate_x_p_y_t(res):
    axis_x = []
    axis_y = []
    ppoa_set = []
    value_m = [MAX_POWER>>1, MAX_POWER>>1, MAX_POWER>>1]
    for i in range(0, MAX_POWER):
        for j in range(1, 2 * MAX_POWER):
            value_p = i / float(MAX_POWER)
            value_t = j
            ppoa = calculate_equilibrium(res, value_m, value_t, value_p)
            axis_x.append(i)
            axis_y.append(j)
            ppoa_set.append(ppoa)
    logger.info('#'*Placeholder)
    logger.info(f"Fixing m0, m1, m2, the x axis is p and the y axis is t")
    logger.info(ppoa_set)
    logger.info('#'*Placeholder)
    logger.info("")
    print3d(axis_x, axis_y, ppoa_set, 'p', 't', 'PPoA', com = f'm0 = {value_m[0]} m1 = {value_m[1]} m2 = {value_m[2]}')
    return

def simulate_equilibrium_3(res, choice = 0):
    assert POOLS == 3
    if choice in [0, 1]:
        simulate_x_m1_y_m2(res)
    if choice in [0, 2]:
        simulate_x_m0m1_y_m2(res)
    if choice in [0, 3]:
        simulate_x_m0m1m2_y_t(res)
    if choice in [0, 4]:
        simulate_x_m2_y_t(res)
    if choice in [0, 5]:
        simulate_x_p_y_t(res)
    return

def simulate_equilibrium(res):
    value_m = [random.randint(1, MAX_POWER) for i in range(POOLS)]
    # value_m = [MAX_POWER for i in range(POOLS)]
    value_t = random.randint(1, MAX_POWER)
    # value_t = MAX_POWER
    value_p = random.random()
    # value_p = 0
    ppoa = calculate_equilibrium(res, value_m, value_t, value_p, show_converge = True)
    return ppoa

if __name__ == "__main__":

    x = sympy.symbols(f'x:{POOLS}:{POOLS}')
    y = sympy.symbols(f'y:{POOLS-1}')
    m = sympy.symbols(f'm:{POOLS}')
    t, p = sympy.symbols(f't, p')
    R = sympy.symbols(f'R:{POOLS}')
    r = sympy.symbols(f'r:{POOLS}')

    E = [(m[i] + sum([get_x(x, j, i) for j in range(POOLS)])) * r[i] - (R[i] + sum([get_x(x, i, j) * r[j] for j in range(POOLS)])) for i in range(POOLS)]
    res = sympy.solve(E, r)
    print_r(res)
    TOTAL = sum(m) + t - sum([(1-p) * get_x(x, i, j) for i in range(POOLS) for j in range(POOLS)])
    for i in range(POOLS):
        tmp = (m[i] - sum([get_x(x, i, j) for j in range(POOLS)]) + sum([p * get_x(x, j, i) for j in range(POOLS)])) / TOTAL
        for j in range(POOLS):
            res[r[j]] = res[r[j]].subs(R[i], tmp)
    print_r(res)

    if POOLS == 3:
        if len(sys.argv) > 1:
            simulate_equilibrium_3(res, int(sys.argv[1]))
        else:
            simulate_equilibrium(res)
    else:
        simulate_equilibrium(res)
