# Why the next sourced-star move needs a different local integrator

The source-preserving ear operator integrates a polynomial in its local
ratios. That closure does not extend directly to the existing sourced
`SingleMove` chart. A two-loop example already retains a rational
denominator and a nonzero logarithmic interval integral.

This is a **local obstruction to that polynomial interpolation method**.
The complete original kernel in the example is the rational number `-26`.

## Graph, sources, and normalization

Use the ordered oriented graph

```python
edges = ((0,1), (0,2), (0,3), (1,3), (2,3))
marked = (2,3)
sources = ((1,2), (3,-1), (-2,4), (-2,-5))
```

The full source sum is zero. Vertex `0` has contracted degree three and
one surviving field of a quartic superpotential. Vertex `1` has contracted
degree two and one surviving cubic field. Vertices `2,3` are the two marked
insertions. All four nonzero source rows are retained.

Apply `SingleMove(graph, 0, 1, 2, 3)` at the retained times
`(T_12,T_13,T_23)=(2,3,5)`. Its two existing local coordinates are `r1,r2`.
For clarity, in the formulas below `r1,r2` are written as \(r_1,r_2\).
The original times are

\[
\begin{aligned}
t_{01}&=\frac{6(1-r_2)}{2r_1r_2-7r_1-5r_2+10},&
t_{02}&=\frac{10(1-r_1)}{2r_1r_2-7r_1-5r_2+10},\\
t_{03}&=\frac{15}{2r_1r_2-7r_1-5r_2+10},&
t_{13}&=\frac3{r_1},\qquad t_{23}=\frac5{r_2}.
\end{aligned}
\]

The scalar density convention is exactly `htk.CONVENTIONS`, with the
signed `SingleMove.jacobian` in coordinate order
`(r1,r2,T_12,T_13,T_23)`. In this edge order that Jacobian is

\[
-\frac{450}
 {r_1^2r_2^2(2r_1r_2-7r_1-5r_2+10)^3}.
\]

Thus the chart coefficient includes its orientation sign. Replacing the
signed Jacobian by its absolute value reverses the entire local expression;
it does not remove the logarithm.

## Exact local integral

The full transformed density, including every source row, is

\[
\frac{r_1r_2(2r_1r_2-11r_1-9r_2+18)
                  (20r_1r_2-119r_1-71r_2+170)}
 {500(2r_1r_2-7r_1-5r_2+10)^2}.
\]

At \(r_2=1/3\), this becomes

\[
\frac{r_1(31r_1-45)(337r_1-439)}{1500(19r_1-25)^2}.
\]

An exact partial-fraction decomposition gives

\[
\frac{10447r_1}{541500}-\frac{6089}{2572125}
-\frac{30781}{514425(19r_1-25)}
-\frac{112}{6859(19r_1-25)^2}.
\]

Consequently,

\[
\int_0^1
\frac{r_1(31r_1-45)(337r_1-439)}{1500(19r_1-25)^2}\,dr_1
=\frac{147541}{20577000}
+\frac{30781}{9774075}\log\frac{25}{6}.
\]

The logarithmic coefficient is nonzero. The denominator also depends on
the retained core times before their displayed specialization. A general
sourced-star transfer therefore needs an exact rational-function integration
step, or another identity that treats these terms; the existing finite
polynomial/Beta rule alone does not apply in this chart.

## Independent controls and scope

[`check_sourced_star.py`](check_sourced_star.py) reconstructs the complete
transformed density from the original graph numerator, checks it against
the original current determinant at three rational points, and compares
all thirty original-current coordinates at those points. It differentiates
a real exact primitive to check the interval integral. Independently,
generic series-parallel integration and the sourced-ear operator both give
the complete original kernel

\[
-26.
\]

The local logarithm therefore establishes neither irrationality of this
complete kernel nor failure of a full Wick/color cancellation. It identifies
the additional analytic structure a general sourced-star backend must
handle.

Run from the repository root:

```bash
python3 explore/circuit-general/forest-general-reduction/check_sourced_star.py
```
