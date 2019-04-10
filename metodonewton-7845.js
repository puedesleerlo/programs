(md,tex)=>{
    return(
        md`Consideremos un péndulo simple, como el representado en la Figura. Si desplazamos la partícula desde la posición de equilibrio hasta que el hilo forme un ángulo Θ con la vertical, y luego la abandonamos partiendo del reposo, el péndulo oscilará en un plano vertical bajo la acción de la gravedad. Las oscilaciones tendrán lugar entre las posiciones extremas Θ y -Θ, simétricas respecto a la vertical, a lo largo de un arco de circunferencia cuyo radio es la longitud, ${tex`\ell`} , del hilo. El movimiento es periódico, pero no podemos asegurar que sea armónico.
        Para determinar la naturaleza de las oscilaciones deberemos escribir la ecuación del movimiento de la partícula.
        La partícula se mueve sobre un arco de circunferencia bajo la acción de dos fuerzas: su propio peso (mg) y la tensión del hilo (N), siendo la fuerza motriz la componente tangencial del peso. Aplicando la segunda ley de Newton obtenemos:
        ${tex.block`F_{\text{t}}=-mg\sin{\theta }=ma_{{\text{t}}}\,`}
        siendo at, la aceleración tangencial y donde hemos incluido el signo negativo para manifestar que la fuerza tangencial tiene siempre sentido opuesto al desplazamiento (fuerza recuperadora).
        Al tratarse de un movimiento circular, podemos poner
        ${tex.block`a_{{\text{t}}}=\ell {\ddot  \theta }\,`} 
        siendo ${tex`{\ddot  \theta }\,`} la aceleración angular, de modo que la ec. dif. del movimiento es:
        ${tex.block`-mg\sin \theta =m\ell {\ddot  \theta } \Rightarrow  \ell {\ddot  \theta }+g\sin \theta =0\,`}
        Esta ec. dif. no corresponde a un movimiento armónico simple (m.a.s.) debido a la presencia de la función seno, de modo que podemos asegurar que el movimiento del péndulo simple no es armónico simple, en general.`)
};