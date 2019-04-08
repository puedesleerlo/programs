(geom) => {return(
    function PoincareCircleOutput(svg, x, y, scale, k=20) {
    //   var circle = ProjectionInput(svg, x, y, scale, 1.0);
    //   var lines = [];
    //   var segments = [];
      
      var shape = svg.append("g");
      shape.append("circle")
        .attr("r", scale)
        .attr("cx", x)
        .attr("cy", y)
        .attr("stroke","red")
        .attr("stroke-width", 1.5)
        .attr("fill", "none");
      
      function NewSegment(p1, p2, color, stroke) {
        shape.append("line")
          .attr("stroke", color)
          .attr("stroke-width", stroke)
          .attr("x1", x + scale * p1[0])
          .attr("y1", y + scale * p1[1])
          .attr("x2", x + scale * p2[0])
          .attr("y2", y + scale * p2[1]);
      }
      
      function PoincareSegmentPoint(l,t) {
        if (l.c != Infinity) {
          var a = l.thetap * (1-t) + l.thetaq * t;
          return geom.Translate(l.c,
                           geom.Scale(geom.Translate(geom.Scale(l.pn, Math.cos(a)),
                                           geom.Scale(l.qn, Math.sin(a))),
                                 l.r));
        } else {
          return geom.Translate(geom.Scale(l.p, 1-t),
                           geom.Scale(l.q, t));
        }
      }
      
      function PoincareLinePoint(l,t) {
        if (l.c != Infinity) {
          var a = -l.theta * (1-t) + l.theta * t;
          return geom.Translate(l.c,
                           geom.Scale(geom.Translate(geom.Scale(l.pn, Math.cos(a)),
                                           geom.Scale(l.qn, Math.sin(a))),
                                 l.r));
        } else {
          return geom.Translate(geom.Scale(l.pn, 1-t),
                           geom.Scale(l.qn, t));
        }
      }
    
      function AddLine(l, color="blue", stroke = 2) {
        if (l.c != Infinity) {
          var k = Math.ceil((2 * l.theta * l.r) / 0.1);
          if (k == 0) console.log("k=0");
          for (var i = 0; i < k; ++i) {
            NewSegment(PoincareLinePoint(l, i/k),
                       PoincareLinePoint(l, (i+1)/k),
                       color, stroke);
          }
        } else {
          NewSegment(PoincareLinePoint(l, 0),
                     PoincareLinePoint(l, 1),
                     color, stroke);
        }
      }
      
      function AddSegmentLine(l, color="blue", stroke = 2 ) {
       if (l.c != Infinity) {
          var k = Math.ceil((Math.abs(l.thetaq - l.thetap) * l.r) / 0.1);
          for (var i = 0; i < k; ++i) {
            NewSegment(PoincareSegmentPoint(l, i/k),
                       PoincareSegmentPoint(l, (i+1)/k),
                       color, stroke);
          }
        } else {
          NewSegment(PoincareSegmentPoint(l, 0),
                     PoincareSegmentPoint(l, 1),
                     color, stroke);
        }
      }
      
      function AddSegment(p,q, color="blue", stroke = 2) {
        var l = geom.PoincareLine(p, q);
        AddSegmentLine(l, color, stroke);
      }
      
      function AddPoint(p) {
        var point = shape.append("circle")
          .attr("r", 3.5)
          .attr("cx", x + scale * p[0])
          .attr("cy", y + scale * p[1]);
      }
      
      return {
        AddPoint: AddPoint,
        AddLine: AddLine,
        AddSegment: AddSegment,
        AddSegmentLine: AddSegmentLine,
      }
    }
    )};