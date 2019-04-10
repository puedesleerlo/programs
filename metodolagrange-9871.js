(md,tex) => {
    return(
        md`El lagrangiano del sistema es
        ${tex.block`{\mathcal  {L}}=T-V={\frac  {1}{2}}ml^{2}{\dot  {\theta }}^{2}+mgl\cos {\theta }`}
        donde ${tex`\theta\,`} es la elongación angular (ángulo que forma el hilo con la vertical) y ${tex`l\,`} es la longitud del hilo. Aplicando las ecuaciones de Lagrange se sigue
        ${tex.block`{\frac {d}{dt}}{\frac {\partial {\mathcal {L}}}{\partial {\dot {\theta }}}}-{\frac {\partial {\mathcal {L}}}{\partial \theta }}=0 \Rightarrow ml^{2}{\ddot {\theta }}+mgl\sin \theta =0`}
        y obtenemos la ecuación del movimiento es
        ${tex.block`l{\ddot  {\theta }}+g\sin {\theta }=0`}
        de modo que la masa no interviene en el movimiento de un péndulo.`
    )
};