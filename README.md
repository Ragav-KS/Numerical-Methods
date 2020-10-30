# Numerical-Methods
Python Algorithms for numerically solving problems which are difficult to solve analytically

## How to Use
Just place the modules in your code directory import them in your code

```python
import Roots
import LinEq 
import Differentiation
```
 and use the functions as it is:
 
 For example,
 ```python
f = lambda x: (np.exp(x) - 4*x)
root, table = Roots.Newton_Raphson_Method(f, fd, x=0)
```

