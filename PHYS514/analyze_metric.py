"""Analyze a given metric. 
"""

from metric_handler import * 
import sympy as sp 
import math 

analysis = MetricAnalysis() 

problem_2 = False  
problem_3 = True 
problem_4 = False

if problem_2: 
    g_schwarzchild = analysis.init_schwarzchild()
    schwarzchild_analysis = analysis.analyze_metric(g=g_schwarzchild, index_dict=analysis.SCHWARZCHILD)
    analysis.print_results(g_schwarzchild, analysis.SCHWARZCHILD)

    g_deSitter = analysis.init_deSitter()
    deSitter_analysis = analysis.analyze_metric(g=g_deSitter, index_dict=analysis.DESITTER)
    analysis.print_results(g_deSitter, analysis.DESITTER)

    einstein_deSitter = deSitter_analysis['einstein'] 
    product = 0.0 
    for i in range(len(analysis.DESITTER)):
        for j in range(len(analysis.DESITTER)):
            m = analysis.DESITTER[i] 
            n = analysis.DESITTER[j] 
            mn = m+n
            try: 
                product += (g_deSitter[mn])**(-1)*einstein_deSitter[mn] 
            except ZeroDivisionError: 
                print('g[{mn}] is zero, cannot evaluate.'.format(mn=m+n))
    print('For de Sitter space, the product (g^mn)*(G_mn) is {product}.'.format(product=product))

if problem_3:
    t = sp.symbols('t')
    r = sp.symbols('r')
    theta = sp.symbols('theta')
    phi = sp.symbols('phi')
    a = sp.symbols('a')
    J = sp.symbols('J')
    M = sp.symbols('M')
    G = sp.symbols('G')

    g_kerr = analysis.init_kerr()
    print('Sanity check on the Kerr metric...') 
    for key in g_kerr:
        print('\tg[{key}] = {expr}'.format(key=key, expr=sp.latex(g_kerr[key])))
    print('Sanity check on inverting the metric...') 
    g_kerr_inv = analysis.invert_metric(g_kerr, analysis.KERR)
    for key in g_kerr_inv:
        print('\tg_inv[{key}] = {expr}'.format(key=key, expr=sp.latex(g_kerr_inv[key])))
    simplify = {'Gamma':False, 'riemann':False, 'ricci_tensor':False, 'ricci_scalar':False, 'einstein':False}
    kerr_analysis = analysis.analyze_metric(g=g_kerr, index_dict=analysis.KERR, simplify=simplify)
    #analysis.print_results(kerr_analysis)
    einstein = kerr_analysis['einstein'] 
    point = {t:100.0, r:7e8, theta:math.pi/3, phi:math.pi/4, a:100, M:1e42, G:6.67e-11} 
    print('Evaluating the Einstein tensor at this point:')
    for var in point:
        print('\t{var} = {point}'.format(var=var, point=point[var]))
    einstein_eval = analysis.evaluate(name='einstein', symbol=einstein, point=point)
    print('The result is...')
    for key in einstein_eval:
        print('\tG[{key}] = {expr}'.format(key=key, expr=einstein_eval[key]))


