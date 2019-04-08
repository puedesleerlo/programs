(svg,poleX,poleY,efectorX,efectorY) => {
    return(
        svg`<svg style="overflow: visible">    
        <line [attr.x1]="${centroX}" [attr.y1]="${centroY}" [attr.x2]="${poleX}" [attr.y2]="${poleY}" style="stroke:rgb(255,0,0);stroke-width:0.5%" />
        <line [attr.x1]="${poleX}" [attr.y1]="${poleY}" [attr.x2]="${efectorX}" [attr.y2]="${efectorY}" style="stroke:rgb(255,0,0);stroke-width:0.5%" />
        <circle [attr.cx]="${efectorX}" [attr.cy]="${efectorY}" r="2%"   fill="black" />
        <circle [attr.cx]="${centroX}" [attr.cy]="${centroY}" r="2%"  fill="red" />
        <circle [attr.cx]="${poleX}" [attr.cy]="${poleY}" r="2%"  fill="black" />
        </svg>
        `
    )
}