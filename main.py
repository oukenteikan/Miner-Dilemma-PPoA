from sympy import *
from logger import init_logger
logger = init_logger()

Placeholder = 64
Simulation = 32
Extra_Simulation = 64

def log_list(msg):
    for m in msg:
        logger.warning(m)
    
def log_bare(msg):
    if isinstance(msg, list) or isinstance(msg, tuple):
        log_list(msg)
    else:
        logger.warning(msg)

def log_unit(msg):
    logger.warning('#'*Placeholder)
    log_bare(msg)

def log_end():
    log_unit("")

def log_wrap(msg):
    log_unit(msg)
    log_end()

def log_pair(s, t):
    log_unit(s)
    log_bare(t)
    log_end()

def log_set(doing, before, str1, after, str2):
    log_unit(doing)
    log_unit(before)
    log_bare(str1)
    log_unit(after)
    log_bare(str2)
    log_end()
    return
    
def first_derivative_of_r(r, x):
    pr = simplify(expand(simplify(diff(r, x) * ((m1*m2+m1*x1+m2*x2)*(m-(1-p)*(x1+x2)))**2 / (1-p))))
    log_set("Differentiating the reward function ...", "The primitive function r is:", r, f"The first derivative to {x} is:", pr)
    return pr

def second_derivative_of_r(r, x, y):
    ppr = simplify(expand(simplify(diff(diff(r, x), y) * ((m + (p - 1)*(x1 + x2))**3*(m1*m2 + m1*x1 + m2*x2)**3/(1-p)))))
    log_set("Differentiating the reward function ...", "The primitive function r is:", r, f"The second derivative to {x} then {y} is:", ppr)
    return ppr

def get_order(equ, var, deg):
    equ = expand(equ)
    coe = factor(simplify(equ.coeff(var, deg)))
    return coe

def print_by_order(equ, var, deg = 5, com = ""):
    log_unit(com)
    equ = expand(equ)
    for i in range(deg, -1, -1):
        logger.info(f"The {i}th order coefficient of {var} is {get_order(equ, var, i)}")
    log_end()
    return None

def make_U(u11, u12, u21, u22):
    U = [u11, u12, u21, u22]
    for j in range(len(U)):
        U[j] = expand(U[j].subs(m, m1+m2+t))
        print_by_order(U[j], t, com = f"The {j+1}th element of U is a function of t:")
    return U

def first_derivative(equ, var, com = ""):
    log_unit(com)
    pequ = diff(equ, var)
    log_pair(f"The derivative to {var} is:", pequ)
    print_by_order(pequ, var, com = f"The derivative is still a function of {var}:")
    return pequ

def study_on_diagonal_element(u11):
    log_wrap(f"We aim to show the diagonal elements are negative. WLOG, consider u11:")
    u11_t_1 = get_order(u11, t, 1)
    print_by_order(u11_t_1, x1, com = f"The linear coefficient of t in u11 is a function of {x1}:")
    u11_t_1_x1_0 = get_order(u11_t_1, x1, 0)
    print_by_order(u11_t_1_x1_0, x2, com = f"The constant term without x1 in u11_t_1 is a function of {x2}:")
    log_wrap(["The constant term without x1 in u11_t_1 is negative by analysis.",
             "Easy to see both u11_t_2 and u11_t_1 are negative.",
             "Since t>=0, to prove u11<=0, we only need to prove u11_t_0<=0."])
    
    u11_t_0 = get_order(u11, t, 0) / 2 / m2
    print_by_order(u11_t_0, m2, com = f"The constant term without t in u11 is a function of {m2}:")
    u11_t_0_pr_m2 = first_derivative(u11_t_0, m2, com = "Differentiating u11_t_0 to m2 ...")
    u11_t_0_pr_m2_m2_x2 = factor(simplify(u11_t_0_pr_m2.subs(m2, x2)))
    log_unit(["We find this quadratic function opens downwards, the axis of symmetry is on the left of x2.",
              "It obtains its maximum value at x2 because m2 is in [x2, infty)."])
    log_pair("Setting m2=x2, we get function:", u11_t_0_pr_m2_m2_x2)
    print_by_order(u11_t_0_pr_m2_m2_x2, m1, com = f"u11_t_0_pr_m2_m2_x2 is a function of m1:")
    
    log_wrap(["We can prove u11_t_0_pr_m2_m2_x2 decrease on [x1, infty), so it obtains its maximum value at x1.",
              "Setting m1=x1, we find u11_t_0_pr_m2_m2_x2_m1_x1<0. Therefore, u11_t_0_pr_m2<0 holds."])
    u11_t_0_m2_x2 = expand(simplify(u11_t_0.subs(m2, x2)))
    str1 = u11_t_0_m2_x2
    log_pair("Setting m2=x2 in u11_t_0, we get:", str1)
    u11_t_0_m2_x2 = expand(simplify((u11_t_0_m2_x2 - m1*x2*(-m1+x1)**3)/p))
    str2 = u11_t_0_m2_x2
    log_pair("Removing a negative term m1*x2*(-m1+x1)**3 and dividing p from u11_t_0_m2_x2, we get:", str2)
    print_by_order(u11_t_0_m2_x2, m1, com = f"u11_t_0_m2_x2 is a function of m1:")
    u11_t_0_m2_x2_pr_m1 = first_derivative(u11_t_0_m2_x2, m1, com = "Differentiating u11_t_0_m2_x2 to m1 ...")

    log_unit(["We find this quadratic function opens downwards, the axis of symmetry is on the left of x1.",
              "It obtains its maximum value at x1 because m1 is in [x1, infty)."])
    u11_t_0_m2_x2_pr_m1_m1_x1 = factor(expand(u11_t_0_m2_x2_pr_m1.subs(m1, x1)))
    log_pair("Setting m1=x1, we get function:", u11_t_0_m2_x2_pr_m1_m1_x1)
    u11_t_0_m2_x2_m1_x1 = u11_t_0_m2_x2.subs(m1, x1)
    log_pair("We can prove u11_t_0_m2_x2_pr_m1_m1_x1<0 holds, so u11_t_0_m2_x2 decreases on [x1, infty). Setting m1=x1, we get function:", u11_t_0_m2_x2_m1_x1)
    log_wrap(["We can prove u11_t_0_m2_x2_m1_x1<0 holds, so u11_t_0<0 holds. Finally, we have u11<0.",
              "Symmetrically, we have u22<0. Thus we know existence of Nash equilibrium can be derived by the concavity."])
    return None

def study_on_offdiagonal_element(u11, u12, u21):
    S = factor(simplify(u12 + u21))
    log_unit("We aim to show U is a diagonally dominant matrix.")
    log_pair("The sum of off-diagonal elements:", S)
    S1 = expand(2 * u11 - S)
    S2 = expand(2 * u11 + S)
    print_by_order(S1, t, com = "2u11 - u12 - u21 is a function of t:")
    print_by_order(S2, t, com = "2u11 + u12 + u21 is a function of t:")
    log_wrap("Unfortunately, S1 and S2 are not always negetive!!!!!!")
    return None

def study_on_det(U):
    det = factor(expand(4*U[0]*U[3] - (U[1]+U[2])**2))
    log_unit("We aim to show det U is positive.")
    log_pair("det U is:", det)
    log_unit("Simulating det U ...")
    for i in range(1, Simulation):
        for j in range(1, Simulation):
            for k in range(Simulation):
                for y in range(i+1):
                    for z in range(j+1):
                        det_v = det
                        det_v = det_v.subs(p, 0)
                        det_v = det_v.subs(m1, i)
                        det_v = det_v.subs(m2, j)
                        det_v = det_v.subs(t,  k)
                        det_v = det_v.subs(x1, y)
                        det_v = det_v.subs(x2, z)
                        if (int(det_v) < 0):
                            logger.fatal('We find a set of parameter to make det<0!!!!!!')
                            logger.fatal(f"det={det_v} m1={i}, m2={j}, t={k}, x1={y}, x2={z}")
    log_end()
    log_wrap("Unfortunately, det U is not always positive!!!!!!")
    return None

def simulate_rosen(pr1, pr2):
    log_wrap("We aim to show rosen criterion holds.")
    pr1 = pr1.subs(p, 0)
    pr2 = pr2.subs(p, 0)
    pr1 = pr1.subs(m, m1+m2+t)
    pr2 = pr2.subs(m, m1+m2+t)
    log_unit("Simulating Rosen ...")
    for i in range(1, Simulation):
        for j in range(1, Simulation):
            for k in range(Simulation):
                for p1 in range(i+1):
                    for q1 in range(j+1):
                        for p2 in range(i+1):
                            for q2 in range(j+1):
                                pr1_v1 = pr1
                                pr1_v1 = pr1_v1.subs(m1, i)
                                pr1_v1 = pr1_v1.subs(m2, j)
                                pr1_v1 = pr1_v1.subs(t,  k)
                                pr1_v1 = pr1_v1.subs(x1, p1)
                                pr1_v1 = pr1_v1.subs(x2, q1)

                                pr1_v2 = pr1
                                pr1_v2 = pr1_v2.subs(m1, i)
                                pr1_v2 = pr1_v2.subs(m2, j)
                                pr1_v2 = pr1_v2.subs(t,  k)
                                pr1_v2 = pr1_v2.subs(x1, p2)
                                pr1_v2 = pr1_v2.subs(x2, q2)

                                pr2_v1 = pr2
                                pr2_v1 = pr2_v1.subs(m1, i)
                                pr2_v1 = pr2_v1.subs(m2, j)
                                pr2_v1 = pr2_v1.subs(t,  k)
                                pr2_v1 = pr2_v1.subs(x1, p1)
                                pr2_v1 = pr2_v1.subs(x2, q1)

                                pr2_v2 = pr2
                                pr2_v2 = pr2_v2.subs(m1, i)
                                pr2_v2 = pr2_v2.subs(m2, j)
                                pr2_v2 = pr2_v2.subs(t,  k)
                                pr2_v2 = pr2_v2.subs(x1, p2)
                                pr2_v2 = pr2_v2.subs(x2, q2)
                
                                res = (p1 - p2)*(pr1_v1 - pr1_v2) + (q1 - q2)*(pr2_v1 - pr2_v2)
                                if (int(res) > 0):
                                    log_pair('We find a set of parameter to make a counter-example!!!!!!', f"res={res} m1={i}, m2={j}, t={k}, x1*={p1}, x2*{q1}, x1-={p2}, x2-{q2}")
    log_end()
    log_wrap("Unfortunately, Rosen criterion does not always hold!!!!!!")
    return None

def simulate_uniqueness(pr1, pr2):
    log_wrap("We aim to show pure NE is unique.")
    pr1 = pr1.subs(p, 0)
    pr2 = pr2.subs(p, 0)
    pr1 = pr1.subs(m, m1+m2+t)
    pr2 = pr2.subs(m, m1+m2+t)
    log_unit("Simulating uniqueness ...")
    for i in range(1, Simulation):
        for j in range(1, Simulation):
            for k in range(Extra_Simulation):
                pr1_v = pr1
                pr1_v = pr1_v.subs(m1, i)
                pr1_v = pr1_v.subs(m2, j)
                pr1_v = pr1_v.subs(t,  k)

                pr2_v = pr2
                pr2_v = pr2_v.subs(m1, i)
                pr2_v = pr2_v.subs(m2, j)
                pr2_v = pr2_v.subs(t,  k)

                logger.fatal('#'*Placeholder)
                logger.info(f'pr1={pr1_v}')
                logger.info(f'pr2={pr2_v}')
                raw_res = solve([pr1_v, pr2_v], [x1, x2], domain=S.Reals)
                res = []
                logger.info(f'm1={i}, m2={j}, t={k}......')
                for r in raw_res:
                    if r[0] > 0 and r[0] < i and r[1] > 0 and r[1] < j:
                        res.append(r)
                        logger.info(f'x1={float(r[0])}, x2={float(r[1])}')
                if len(res) > 1:
                    logger.fatal(f'Too many solutions????????????????????????????????????????????????')
                if len(res) > 0:
                    pr2_v = pr2_v.subs(x1, 0)
                    x2_raw_res = solve(pr2_v, x2, domain=S.Reals)
                    x2_res = []
                    for x2_r in x2_raw_res:
                        if x2_r >= 0 and x2_r <= j:
                            x2_res.append(x2_r)
                    if len(x2_res) != 1:
                        logger.fatal(f'Something must be wrong????????????????????????????????????????????????')
                    pr1_v = pr1_v.subs(x1, 0)
                    pr1_v = pr1_v.subs(x2, x2_res[0])
                    if pr1_v <= 0:
                        logger.info(f'x1={0}, x2={x2_res[0]}')
                        logger.fatal(f'Another result is found????????????????????????????????????????????????')
                        exit(1)
    log_end()
    log_wrap("Thankfully, it is likely to be not.")
    return None

def study_on_extreme(pr1, pr2):
    log_wrap("Studying on extreme cases...")
    E1 = simplify(pr1.subs(x2, 0))
    log_pair("Consider extreme case when x2=0:", E1)
    E2 = solve(E1, x1)
    log_pair("Solve pr1=0 to get x1, we get:", E2)
    E3 = E2[1]
    log_pair("Getting root for x1", E3)
    E4 = pr2.subs(x1, E3)
    E4 = E4.subs(x2, 0)
    E4 = factor(simplify(E4*(m - m1 + m2*p - m2)**3/(m1*m2**2)))
    E4 = simplify(E4.subs(m, m1+m2+t))
    log_pair("Getting pr2", E4)

    # Fail to remember what the rest of function is doing.
    E5 = (m1*m2*p + m1*t + m2**2 + 2*m2*t + t**2) * (m1*m2**2*p**3 - 2*m1*m2**2*p**2 + 3*m1*m2**2*p + m1*m2*p*t + 3*m1*m2*t - m1*p*t**2 + 3*m1*t**2 - 2*m2**3*p + 4*m2**3 - 4*m2**2*p*t + 10*m2**2*t - 2*m2*p*t**2 + 8*m2*t**2 + 2*t**3) ** 2 - (m1**2*m2**2*p**3 - m1**2*m2**2*p**2 + 2*m1**2*m2*p**2*t - 2*m1**2*m2*p*t + m1**2*p*t**2 - m1**2*t**2 + 3*m1*m2**3*p**2 - 5*m1*m2**3*p - m1*m2**2*p**3*t + 6*m1*m2**2*p**2*t - 6*m1*m2**2*p*t - 5*m1*m2**2*t + 3*m1*m2*p*t**2 - 9*m1*m2*t**2 + m1*p*t**3 - 3*m1*t**3 + 2*m2**4*p - 4*m2**4 + 6*m2**3*p*t - 14*m2**3*t + 6*m2**2*p*t**2 - 18*m2**2*t**2 + 2*m2*p*t**3 - 10*m2*t**3 - 2*t**4) ** 2
    E5 = factor(simplify(E5/(m1*(m1 - m2*p + 2*m2 + t)*(m2*p + t)**3)))
    E6 = diff(E5, t)
    print_by_order(E5, t, com = 'negative required')
    E8 = (1-p)*f**2 - (2*m2+(1-p)*m1+2*t)*f + 2*m1*m2
    E9 = solve(E8, f)
    E9 = E9[1]
    E10 = f**4*p**2*t - 2*f**4*p*t + f**4*t + f**3*m1*m2*p**3 + f**3*m1*m2*p**2 - f**3*m1*m2*p - f**3*m1*m2 + 4*f**3*m1*p*t - 4*f**3*m1*t + 4*f**3*m2*p*t - 4*f**3*m2*t + 4*f**3*p*t**2 - 4*f**3*t**2 + f**2*m1**2*m2*p**2 + 2*f**2*m1**2*m2*p + f**2*m1**2*m2 + 4*f**2*m1**2*t + f**2*m1*m2**2*p**2 + 2*f**2*m1*m2**2*p + f**2*m1*m2**2 + 3*f**2*m1*m2*p**2*t - 2*f**2*m1*m2*p*t + 11*f**2*m1*m2*t + 8*f**2*m1*t**2 + 4*f**2*m2**2*t + 8*f**2*m2*t**2 + 4*f**2*t**3 - f*m1**2*m2**2*p**3 - f*m1**2*m2**2*p**2 + f*m1**2*m2**2*p + f*m1**2*m2**2 + 4*f*m1**2*m2*p*t - 4*f*m1**2*m2*t + 4*f*m1*m2**2*p*t - 4*f*m1*m2**2*t + 4*f*m1*m2*p*t**2 - 4*f*m1*m2*t**2 - m1**3*m2**2*p**2 - 2*m1**3*m2**2*p - m1**3*m2**2 - m1**2*m2**3*p**2 - 2*m1**2*m2**3*p - m1**2*m2**3 - 4*m1**2*m2**2*p*t
    E11 = E10.subs(f, E9)
    E11 = factor(simplify(E11*(2*(p - 1)**2)/(m1*(m1 - m2)*(p + 1)**2)))
    E12 = (m1**2*m2*p**3 - 2*m1**2*m2*p**2 + m1**2*m2*p + m1**2*p**2*t - 2*m1**2*p*t + m1**2*t - 2*m1*m2*p**2*t + 2*m1*m2*t - 4*m1*p*t**2 + 4*m1*t**2 + 4*m2**3 + 12*m2**2*t + 12*m2*t**2 + 4*t**3) ** 2 - (m1**2*p**2 - 2*m1**2*p + m1**2 + 4*m1*m2*p - 4*m1*m2 - 4*m1*p*t + 4*m1*t + 4*m2**2 + 8*m2*t + 4*t**2) * (m1*m2*p**2 - m1*m2*p + m1*p*t - m1*t - 2*m2**2 - 4*m2*t - 2*t**2) ** 2
    E12 = factor(simplify(E12/(4*m1*m2**2*(p - 1)**2)))
    print_by_order(E12, t, com = 'positive required')
    return None

def study_on_symmetric(pr1, pr2):
    log_wrap('Consider symmetric case when m1=m2=m/k.')
    S3 = expand(pr1).subs(m1, m/k)
    S3 = S3.subs(m2, m/k)
    S3 = simplify(S3 * k**4 / m)
    S4 = expand(pr2).subs(m2, m/k)
    S4 = S4.subs(m1, m/k)
    S4 = simplify(S4* k**4 / m)
    log_pair(S3, S4)
    S6 = factor(S3 - S4)
    log_pair('Subtract pr2=0 from pr1=0.', S6)
    S5 = factor((expand(S3) - expand(S4))/(k*(x1 - x2)))
    log_pair('The result can be factorized but only x1=x2 is possible. The discriminant of following equation is negative.', S5)
    S7 = S3.subs(x2, x1)
    S7 = solve(S7, x1)
    log_pair('Setting x2=x1, we solve pr1=0 to get x1. x=-m/(2*k) is impossible. Then the positive root is also impossible.', S7)
    return None

def study_on_binary(G7, G8):
    log_wrap('Consider binary case when t=0.')
    N8 = simplify(G8.subs(t, 0)/(m1*m2*(p + 1)**2))
    log_pair('Setting t=0, we get:', N8)
    print_by_order(N8, f, com = 'g(f) is a still function of f:')
    N9 = solve(N8, f)
    log_pair('This cubic function can be solved to get f.', N9)
    N10 = G7.subs(t, 0)
    N11 = simplify(N10.subs(f, sqrt(m1*m2)))
    log_pair('Only f=sqrt(m1*m2) is possible, so setting f=m1*m2 we get:', N11)
    log_wrap('Only the negative root is possible. We can get x2 symmetrically.')
    return None

def study_on_general(pr1, pr2):
    log_wrap('Consider general case when m1!=m2. We will show the binary case by the way.')
    G5 = simplify((pr1 + pr2) / (m1*m2 + m1*x1 + m2*x2))
    G5 = factor(G5.subs(m, m1+m2+t))
    log_pair('The sum of first derivatives of r.', G5)
    log_unit('Replace (1-p)*(x1+x2)*2 in pr1 and pr2 to lower the degree using pr1+pr2=0')
    G_pr1 = expand(m1*(m2+2*x1)*pr1)
    G_pr2 = expand(m2*(m1+2*x2)*pr2)
    G6 = simplify((G_pr1 - G_pr2)/m1/m2/(m1*m2 + m1*x1 + m2*x2))
    G6 = factor(G6.subs(m, m1+m2+t))
    log_pair('The difference of first derivatives of r.', G6)
    G5 = factor(G5.subs(x2, f-x1))
    G6 = factor(G6.subs(x2, f-x1))
    log_unit('The sum and difference of first derivatives of r using f.')
    log_pair(G5, G6)
    print_by_order(G5, f)
    print_by_order(G6, f)
    G7 = solve(G5, x1)
    log_pair('Solve x1 in pr1+pr2=0 using f.', G7)
    G7 = G7[0]
    G8 = G6.subs(x1, G7)
    G8 = simplify(expand(G8*((p**2 + 2*p + 1)*m1*m2*(m2-m1))/2))
    log_pair('Substituting the above equation we get a quartic function g(f):', G8)
    print_by_order(G8, f, com = 'g(f) is a function of f:')

    study_on_binary(G7, G8)

    log_unit('If t!=0, since g(f) cannot be factorized, we differentiate it twice.')
    G11 = collect(diff(G8, f), f)
    G12 = collect(diff(G11, f), f)
    log_pair(G11, G12)
    print_by_order(G11, f, com = 'The first derivative of g(f) is still a function of f')
    print_by_order(G12, f, com = 'The second derivative of g(f) is still a function of f')
    G13 = solve(G12, f)
    log_pair("Solving g''(f)=0, we get:", G13)
    G14 = 9*m1**2*m2**2*p**4 + 36*m1**2*m2**2*p**3 + 54*m1**2*m2**2*p**2 + 36*m1**2*m2**2*p + 9*m1**2*m2**2 + 48*m1**2*m2*p**2*t + 96*m1**2*m2*p*t + 48*m1**2*m2*t + 48*m1**2*t**2 + 48*m1*m2**2*p**2*t + 96*m1*m2**2*p*t + 48*m1*m2**2*t + 192*m1*m2*p*t**2 + 96*m1*m2*t**2 + 96*m1*t**3 + 48*m2**2*t**2 + 96*m2*t**3 + 48*t**4
    print_by_order(G14, t, com = "The discriminant is function of t:")
    G15 = G11.subs(f, (m1+m2)/2)
    G15 = collect(simplify(expand(G15)), t)
    log_pair("g'(f) increases, then decreases. The middle point is positive.", G15)
    print_by_order(G15, t, com = 'The middle point is a function of t:')
    G16 = G8.subs(f, (m1+m2)/2)
    G16 = collect(expand(G16), t)
    log_pair("g(f) increases only in a continuous interval containing (m1+m2)/2. The middle point is positve.", G16)
    print_by_order(G16, t, com = 'The middle point is a function of t:')
    log_wrap("g(0)<0 and g(m1+m2)>0, so there exists only one root in (0, m1+m2) and it is strictly less than (m1+m2)/2.")
    return None

if __name__ == "__main__":

    # Symbols
    x1, x2, m1, m2, t, m, k, p, f, d = symbols('x1 x2 m1 m2 t m k p f d')

    # Average Rewards
    r1 = ((m1-x1+p*x2)*(m2+x1) + x1*(m2-x2+p*x1)) / ((m1*m2+m1*x1+m2*x2)*(m-(1-p)*(x1+x2)))
    r2 = ((m2-x2+p*x1)*(m1+x2) + x2*(m1-x1+p*x2)) / ((m2*m1+m2*x2+m1*x1)*(m-(1-p)*(x2+x1)))

    # First Partial
    pr1 = first_derivative_of_r(r1, x1)
    pr2 = first_derivative_of_r(r2, x2)

    # Scond Partial
    u11 = second_derivative_of_r(r1, x1, x1)
    u12 = second_derivative_of_r(r1, x1, x2)
    u21 = second_derivative_of_r(r2, x2, x1)
    u22 = second_derivative_of_r(r2, x2, x2)
    U = make_U(u11, u12, u21, u22)

    # Show Concavity
    study_on_diagonal_element(U[0])
    study_on_offdiagonal_element(U[0], U[1], U[2])
    # study_on_det(U)
    # simulate_rosen(pr1, pr2)
    # simulate_uniqueness(pr1, pr2)

    # Main Study
    study_on_extreme(pr1, pr2)
    study_on_symmetric(pr1, pr2)
    study_on_general(pr1, pr2)




