# Basis Analysis

1. Solve the baseline LP and record your optimal **basis** $B$.
2. Form the **basic matrix** $A_B$ (columns of $A$ in $B$) and compute

   $$x^*_B = A_B^{-1}D^0,
   \quad
   M = A_B^{-1}.$$
3. By definition,

   $$\mathcal K_B
   =\{\;\delta D: M\,\delta D \;\ge\;-\,x^*_B\}.$$

   That’s just a system of linear inequalities (an H-representation of your cone).
4. **Practical handling**:

   * **Half-space form**: store the rows of $M$ and the rhs vector $-x^*_B$.  Any candidate shock $\delta D$ is in $\mathcal K_B$ iff every row-inequality holds.
   * **Independent bounds**: for each coordinate $j$, you can get a “safe” interval

     $$ \delta D_j\in
     \bigl[\max_{i:M_{i j}>0}-\tfrac{x^*_{B,i}}{M_{i j}},
           \;\min_{i:M_{i j}<0}-\tfrac{x^*_{B,i}}{M_{i j}}\bigr].$$

     These give a hyper‐rectangle inside $\mathcal K_B$.
   * **Solver-built sensitivity**: most LP solvers can directly report allowable increases/decreases in each RHS so that the basis stays fixed—they’re just the 1-d projections of $\mathcal K_B$.
   * **Full polyhedral processing**: if you need the extremal rays or vertices of the cone, feed the H-rep into a polyhedral library (e.g. cddlib, polyhedra.jl) to get a V-representation or generate extreme directions.

With those steps you can efficiently test or sample any sequence $\{\delta D^k\}$ against $\mathcal K_B$ and guarantee the linear cost‐response remains valid.
