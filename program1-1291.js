(md,tex) => {
    return (
        md`A dot product is what allow us to compute lengths and angles. We compute the length of a vector ${tex`u`} (say in ${tex`\mathbb R^3`}) we compute its length as ${tex`\Vert u \Vert = \sqrt{\langle u, u \rangle}`}. Given another vector ${tex`v`}, we compute the angle as:
${tex.block`\begin{aligned}
\measuredangle(u,v) = \operatorname{acos} \left( \frac{\langle u, v \rangle}{\Vert u \Vert \cdot \Vert v \Vert} \right)
\end{aligned}`}
Now that we can compute lengths and angles of vectors, we can compute lengths and angles or curves. Given a curve ${tex`\gamma : [0,1] \rightarrow \mathbb R^3`} we can compute its length as:
${tex.block`\begin{aligned}
{\sf length}(\gamma) = \int_0^1 \sqrt{\langle \gamma'(t), \gamma'(t) \rangle} dt
\end{aligned}`}
And if two curves intersect at a given point, we can compute the angles between curves as the angle between tangent vectors at that point.

We usually rely on the usual dot product ${tex`{\langle u, v \rangle}_E = u_x v_x + u_y v_y + u_z v_z `} which gives us the Euclidian geometry. If we change the dot product to something else, we can get other very interesting geometries. The hyperbolic geometry arises when we start measuring lengths and angles using the Minkowski product:
${tex.block`\begin{aligned}
{\langle u, v \rangle}_M = u_x v_x + u_y v_y - u_z v_z\end{aligned}`}
It is very similar to the standard Euclidian product except that the sign of the last coordinate is reversed. This product is very important in relativity (where the ${tex`z`}-coordinate is usually the time and this product measures what is called [proper time](https://en.wikipedia.org/wiki/Proper_time)). But even if you don't care about physics, there is some very interesting geometry that comes out of it (and produces nice pictures, like the one on the top).`
)};