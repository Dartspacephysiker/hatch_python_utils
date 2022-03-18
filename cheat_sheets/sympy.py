

# Derivative
from sympy import diff, sin, exp 
from sympy.abc import x,y 
expr=x*sin(x*x)+1
diff(expr,x)
