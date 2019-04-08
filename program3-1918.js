(d3,DOM,PoincareCircleOutput,HyperbolicTesselation) =>
{
  var width = 640;
  var height = 260; 
  var svg = d3.select(DOM.svg(width, height))  
      .style("overflow", "visible"); 
   
  var mixed_input = PoincareCircleOutput(svg, 300, 130, 120, 10);

  var tesselation = HyperbolicTesselation([3,5,5], 5);
  var vertices = [];
  var faces = [];

  tesselation.faces.forEach(function(f) {    
    mixed_input.AddSegment(f.v[0].pos, f.v[1].pos);
    mixed_input.AddSegment(f.v[1].pos, f.v[2].pos);
    mixed_input.AddSegment(f.v[2].pos, f.v[0].pos);
    
  });
  
  
  tesselation.vertices.forEach(function(v) {
    var nfaces = 2 * tesselation.k[v.i];
    for (var i = 0; i < nfaces; ++i) {
      if (v.f[i] != null & v.f[(i+1)%nfaces] != null) {
        mixed_input.AddSegment(v.f[i].center, v.f[(i+1)%nfaces].center, "red");
      }
    }
  });

  return svg.node();
};