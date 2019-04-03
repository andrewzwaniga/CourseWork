"""
Analysis of the FRW metric for problem set 10, problem 1b. 

Compute the Einstein tensor of the FRW metric, given by

\[
ds^2 = -dt^2 + r^2d\sigma_\kappa^2 
\]

where 

\[
d\sigma_\kappa^2 = \frac{dr^2}{1-\kappa r^2} + r^2 d\Omega^2
\]

and $d\Omega^2 = d\theta^2 + \sin^2\theta d\phi^2$.  
"""

from metric import Metric
from christoffel import Christoffel
from riemann import Riemann
from ricci import Ricci
from scalar import Scalar
from einstein import Einstein
from metric_analysis import MetricAnalysis

import sympy as sp
import math 

analysis = MetricAnalysis() 

simplify = {'Gamma':True, 'R':True, 'P':True, 'Q':True, 'G':True}

g_FRW = analysis.init_FRW() 
g_FRW_inv = analysis.invert_metric(g=g_FRW)
FRW_analysis = analysis.analyze_metric(g=g_FRW, g_inv=g_FRW_inv, simplify=simplify)
analysis.print_results(FRW_analysis) 


