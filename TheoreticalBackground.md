This project attempts to find a solution to laplace equation $\nabla^2 \Psi(r) = 0$ for a given electrode geometry. In particular, we are interested in quadrapole geometry used in electron microscopy or linear paul traps. 
Project uses established numerical methods to find such a solution. In particular, we use Finite Elements Methods (FEM) or Boundary Elements Methods (BEM) to find the physical solution, and then we interpolate the result using solid harmonics basis. 

The following text describes the technical details of each step. 



## 2. Boundary Elements Methods



## 3. Interpolation 

The BEM returns us function of potential, stored as a set of ordered pairs of coordinate:potential. Mathematically, let position be $\mathbf{r} = (x,y,z)$, with $(x,y,z)$ being the coordinates in the cartesian coordinate basis and $u$ be the potential. If an ordered set of coordinates be $X = \{ \mathbf{r}_i \mid i \in \{0, \dots, N\} \}$, and ordered set of potential be $Y = \{ u_i \mid i \in \{0, \dots, N\} \}$, then the function of potential $f(X) = Y$ maps linearly all points from potential Y to X. 


$$ Fitting the solid harmonics basis and the linear algebra equation problem Ac = y


### Code implementation
In code, the following relation is expressed by 

(C++ code snippet)







there is a vector function 

