(d3) => {return(
    (function() {
      
      function ApplyToSegment(f,s) {
        return {p1: f(s.p1), p2: f(s.p2)}
      }
      
      function Dot(u,v) {
        var dp = 0;
        for (var i = 0; i < u.length; ++i) {
          dp += u[i] * v[i];
        }
        return dp;
      }
      
      function Norm(v) {
        return Math.sqrt(Dot(v,v));
      }
      
      function Scale(p,s) { 
        var q = [];
        for (var i = 0; i < p.length; ++i) {
          q.push(p[i] * s);
        }
        return q;
      }
      
      function ScaleCoordinates(p,s) { 
        var q = [];
        for (var i = 0; i < p.length; ++i) {
          q.push(p[i] * s[i]);
        }
        return q;
      }
      
      function Translate(p, dp) { 
        var q = [];
        for (var i = 0; i < p.length; ++i) {
          q.push(p[i] + dp[i]);
        }
        return q;
      }
      
      function Orthonormalize(u,v) {
        var un = Scale(u, 1 / Norm(u));
        var p = Translate(v, Scale(un, -Dot(un,v)));
        return Scale(p, 1 / Norm(p));
      }
      
      // Finds a 2x2 determinant
      function Det2(M) {
        return M[1][1] * M[0][0] - M[0][1] * M[1][0];
      }
      
      // Solves a 2x2 linear system. Assumes Det(M) is not too small.
      function LinearSolve2(M, z) {
        var x = [ M[1][1] * z[0] - M[0][1] * z[1],
                  z[1] * M[0][0] - z[0] * M[1][0] ];
        x = Scale(x, 1 / Det2(M));
        return x;
      }
      
      function EuclidianDistance(p,q) {
        return Norm(Translate(p, Scale(q, -1)));
      }
      
      function PoincareCenter(p,q) {
        if (Math.abs(Det2([p,q])) < 0.001) {
          return Infinity;
        } else {
          return LinearSolve2([p,q], [(Dot(p,p) + 1 )/2, (Dot(q,q) + 1 )/2]);
        }
      }
      
      function PoincareCenterDir(p,d) {
        if (Math.abs(Det2([p,d])) < 0.001) {
          return Infinity;
        } else {
          return LinearSolve2([p,d], [(Dot(p,p) + 1 )/2, Dot(p,d)]);
        }
      }
      
      function Acos(x) {
        return Math.acos(Math.max(Math.min(x,1), -1));
      }
      
      function PoincareLine(p,q) {
        var c = PoincareCenter(p,q);
        if (c == Infinity) {
          return {
            c:Infinity,
            cp: Scale(p , 1/Norm(p)),
            cq: Scale(p , -1/Norm(p)),
            p: p,
            q: q
          }
        } else {
          var r = Norm(Translate(p, Scale(c, -1)));
          var cp1 = Translate(
              Scale(c, 1/ ( Norm(c)**2)),
              Scale([-c[1],c[0]], Math.sqrt(1-1/(Norm(c)**2))/Norm(c))
          );
          var cp2 = Translate(
              Scale(c, 1/ ( Norm(c)**2)),
              Scale([c[1],-c[0]], Math.sqrt(1-1/(Norm(c)**2))/Norm(c))
          );
          
          var pn = Translate(cp1, Scale(c, -1));
          var qn = Translate(cp2, Scale(c, -1));
          pn = Scale(pn, 1 / Norm(pn));
          qn = Scale(qn, 1 / Norm(qn));
          var theta = Acos(Dot(pn, qn));
          qn = Orthonormalize(pn, qn);
          var vp = Translate(p, Scale(c, -1));
          var thetap = Acos(Dot(pn, vp) / Norm(vp));
          var vq = Translate(q, Scale(c, -1));
          var thetaq = Acos(Dot(pn, vq) / Norm(vq));
          return {
            c: c,
            pn: pn,
            qn: qn,
            r: r,
            theta: theta,
            thetap: thetap,
            thetaq: thetaq,
          }
        }
      }
      
      function PoincareLineDir(p,d) {
        var c = PoincareCenterDir(p,d);
        if (c == Infinity) {
          return {
            c:Infinity,
            cq: Scale(Translate(p,Scale(d,+1.5)) ,  1/Norm(Translate(p,Scale(d,+1.5)))),
            cp: Scale(Translate(p,Scale(d,-1.5)) , 1/Norm(Translate(p,Scale(d,-1.5)))),
            p: p,
          }
        } else {
          var r = Norm(Translate(p, Scale(c, -1)));
          var cp1 = Translate(
              Scale(c, 1/ ( Norm(c)**2)),
              Scale([-c[1],c[0]], Math.sqrt(1-1/(Norm(c)**2))/Norm(c))
          );
          var cp2 = Translate(
              Scale(c, 1/ ( Norm(c)**2)),
              Scale([c[1],-c[0]], Math.sqrt(1-1/(Norm(c)**2))/Norm(c))
          );
          
          var pn = Translate(cp1, Scale(c, -1));
          var qn = Translate(cp2, Scale(c, -1));
          pn = Scale(pn, 1 / Norm(pn));
          qn = Scale(qn, 1 / Norm(qn));
          if (Dot(d, qn) < 0) {
            var temp = pn;
            pn = qn;
            qn = temp;
          }
          var theta = Acos(Dot(pn, qn));
          qn = Orthonormalize(pn, qn);
          var vp = Translate(p, Scale(c, -1));
          var thetap = Acos(Dot(pn, vp) / Norm(vp));
          return {
            c: c,
            pn: pn,
            qn: qn,
            r: r,
            theta: theta,
            thetap: thetap,
          }
        }
      }
      
      function PoincareDistance(p,q) {
        return Math.acosh(1 + (2 * EuclidianDistance(p,q)**2) / ((1-Dot(p,p))*(1-Dot(q,q))) );
      }
      
      function PoincarePoint(l,t) {
        if (l.c != Infinity) {
          var a = l.thetap * (1-t) + l.theta * t;
          return Translate(l.c,
                           Scale(Translate(Scale(l.pn, Math.cos(a)),
                                           Scale(l.qn, Math.sin(a))),
                                 l.r));
        } else {
          return Translate(Scale(l.p, 1-t),
                           Scale(l.cq, t));
        }
      }
      
      function PoincareSegmentPoint_(l,t) {
        if (l.c != Infinity) {
          var a = l.thetap * (1-t) + l.thetaq * t;
          return Translate(l.c,
                           Scale(Translate(Scale(l.pn, Math.cos(a)),
                                           Scale(l.qn, Math.sin(a))),
                                 l.r));
        } else {
          return Translate(Scale(l.p, 1-t),
                           Scale(l.q, t));
        }
      }
      
      function PoincareLinePoint_(l,t) {
        if (l.c != Infinity) {
          var a = l.theta * t;
          return Translate(l.c,
                           Scale(Translate(Scale(l.pn, Math.cos(a)),
                                           Scale(l.qn, Math.sin(a))),
                                 l.r));
        } else {
          return Translate(Scale(l.cp, 1-t),
                           Scale(l.cq, t));
        }
      }
      
      
      function PoincareSegmentPoint(p,q,t) {
        var l = PoincareLine(p,q);
        var d = PoincareDistance(p, q);
        var ts = BinarySearch(function(s) {
          return PoincareDistance(p, PoincareSegmentPoint_(l,s)) - d*t;
        }, 0.001, 0,1);
        return PoincareSegmentPoint_(l,ts);    
      }
      
      
      function TangentVector(l,t) {
        if (l.c != Infinity) {
          var a = Math.PI/2 + (l.thetap * (1-t) + l.theta * t);
          return Translate(Scale(l.pn, Math.cos(a)),
                           Scale(l.qn, Math.sin(a)));
        } else {
          var v = Translate(Scale(l.p, -1),
                            Scale(l.cq, 1));
          return Scale(v, 1 / Norm(v));
        }
      }
      
      function BinarySearch(f, epsilon, a0, b0) {
        var a = a0, b = b0;
        while (b - a > epsilon) {
          var x = (b+a)/2;
          if (f(x) < 0) {
            a = x;
          } else {
            b = x;
          }
        }
        return (a+b)/2;
      }
      
      
      function PoincareSegmentDir(p, d, dist, eps = 0.001) {
        var l = PoincareLineDir(p,d);
        var t = BinarySearch(s => {return PoincareDistance(p, PoincarePoint(l,s)) - dist; },
                             eps, 0, 1);
        if (l.c == Infinity) {
          l.q = PoincarePoint(l,t);
          l.tq = TangentVector(l,t);
        } else {
          l.thetaq = l.thetap * (1-t) + l.theta * t;
          l.q = PoincarePoint(l,t);
          l.tq = TangentVector(l,t);
        }
        return l;
      }
      
      function Rotate(v, theta) {
        return [ Math.cos(theta) * v[0] - Math.sin(theta) * v[1],
                 Math.sin(theta) * v[0] + Math.cos(theta) * v[1] ];
      }
      
      function CrossProd(p,q) {
        return p[0] * q[1] - p[1] * q[0];
      }
      
      function SegmentIntersect(a,b,c,d) {
        var ab = PoincareLine(a,b);
        var cd = PoincareLine(c,d);
        if (ab.c != Infinity) {
          var sign = (Norm(Translate(c, Scale(ab.c, -1))) < ab.r) ? +1 : -1;
          function pred(t) {
            var p = PoincareSegmentPoint_(cd,t);
            return sign * (Norm(Translate(p, Scale(ab.c, -1))) - ab.r);
          }
          var t = BinarySearch(pred, 0.0001, 0, 1);
          return PoincareSegmentPoint_(cd,t);
        } else {
          var sign = CrossProd(Translate(a, Scale(b,-1)),
                               Translate(c, Scale(b, -1))) < 0 ? +1 : -1;
          function pred(t) {
            var p = PoincareSegmentPoint_(cd,t);
            return sign * CrossProd(Translate(a, Scale(b,-1)),
                                    Translate(p, Scale(b, -1)));
          }
          var t = BinarySearch(pred, 0.0001, 0, 1);
          return PoincareSegmentPoint_(cd,t);
        }
      }
      
        
      function LineIntersect(ab, cd) {
        if ((ab.c == Infinity) && (cd.c == Infinity)) {
          return [0,0];
        }
        if (ab.c != Infinity) {
          
          var sign = (Norm(Translate(PoincareLinePoint_(cd,0), Scale(ab.c, -1))) < ab.r) ? +1 : -1;
          function pred(t) {
            var p = PoincareLinePoint_(cd,t);
            return sign * (Norm(Translate(p, Scale(ab.c, -1))) - ab.r);
          }
          var t = BinarySearch(pred, 0.0001, 0, 1);
          return PoincareLinePoint_(cd,t);
        } else {
          var sign = (Norm(Translate(PoincareLinePoint_(ab,0), Scale(cd.c, -1))) < cd.r) ? +1 : -1;
          function pred(t) {
            var p = PoincareLinePoint_(ab,t);
            return sign * (Norm(Translate(p, Scale(cd.c, -1))) - cd.r);
          }
          var t = BinarySearch(pred, 0.0001, 0, 1);
          return PoincareLinePoint_(ab,t);
        }
      }
      
      function HyperbolicCenter(p,q,r) {
        var pq = PoincareSegmentPoint(p,q,0.5);
        var qr = PoincareSegmentPoint(q,r,0.5);
        return SegmentIntersect(pq,r,p,qr);
      }
    
      
      return {
        // Given p = (\theta, \phi) returns a point in the sphere
        SphericalCoordinate: (p) => {
          var z1 = p[0];
          var z2 = - p[1];
          return [Math.cos(z1) * Math.cos(z2),
                  Math.sin(z1) * Math.cos(z2),
                  Math.sin(z2)];
        },
        HyperbolicalCoordinate: (p) => {
          var z1 = p[0];
          var z2 = p[1];
          return [Math.cos(z1) * Math.sinh(z2),
                  Math.sin(z1) * Math.sinh(z2),
                  Math.cosh(z2)];
        },
        StereographicProjection: (p) => {
          return [p[0] / (1+p[2]), p[1] / (1+p[2])];
        },
        StereographicProjectionInverse: (p) => {
          var n = p[0] * p[0] + p[1] * p[1];
          var t = 2 / (n + 1);
          var np = [p[0] * t, p[1] * t, -1 +t];
          return np;
        },
        StereographicProjectionHyperboloidInverse: (p) => {
          var n = p[0] * p[0] + p[1] * p[1];
          var t = 2 / (-n + 1);
          var np = [p[0] * t, p[1] * t, -1 +t];
          return np;
        },
        // Rotate a (x,y,z) point with respect in plane-yz
        Rotateyz : (p, theta) => {
          var angle = theta * 2 * Math.PI / 360;
          return [p[0], Math.cos(angle) * p[1] + Math.sin(angle) * p[2],
                       - Math.sin(angle) * p[1] + Math.cos(angle) * p[2]];
        },
        // Translates a point. Assumes point and translation are the same size
        Translate: Translate,
        // Scales a point p by a scalar s
        Scale: Scale,
        // Scales a point p by s coordinate-wise
        ScaleCoordinates: ScaleCoordinates,
        // Dot product
        Dot: Dot,
        // 2-Norm of a vector
        Norm: Norm,
        // Return a vector v' that is orthogonal to v such that u and v' span the same space as u and v
        Orthonormalize: Orthonormalize,
        // Projects a point p into the z=1 plane. The scaling is kept as the third point for re of ference
        Project : (p) => { return  [p[0] / (1-p[2]), p[1] / (1-p[2]), p[2] /*1 / (1-p[2])*/]; },
        // Apply a function to a segment
        ApplyToSegment: ApplyToSegment,
        // Apply to a set of segments
        ApplyToSegmentSet: (f,ss) => { 
          return ss.map(s => { return ApplyToSegment(f,s); }) },
        // Grid
        Grid: (k1, k2) => {
          var segments = [];
          for (var i = 0; i < k1; ++i) {
            for (var j = 0; j < k2; ++j) {
              segments.push({p1: [i/k1,j/k2], p2: [(i+1)/k1,j/k2]});
              segments.push({p1: [i/k1,j/k2], p2: [i/k1,(j+1)/k2]});
            }
          }
          return segments;
        },
        // Given a set of pairs of 3d, returns the range of the z-coordinate
        DepthRange: ss => {
          var maxz = Math.max(...ss.map(s => {return Math.max(s.p1[2], s.p2[2]);}));
          var minz = Math.min(...ss.map(s => {return Math.min(s.p1[2], s.p2[2]);}));
          return [minz, maxz];
        },
        // Given a range [minz, maxz] and a vector of 3 components representing a rgb color
        // we return a function that maps a point to its color
        DepthColor: (zrange, rgb) => {
          return p => {
            var c = 200 - 150 * (p[2] - zrange[0]) / (zrange[1] - zrange[0]);
            return d3.rgb(c * rgb[0], c * rgb[1], c * rgb[2]);
          }
        },
        PoincareLine: PoincareLine,
        PoincareLineDir: PoincareLineDir,
        PoincareSegmentDir: PoincareSegmentDir,
        PoincareDistance: PoincareDistance,
        PoincareSegmentPoint: PoincareSegmentPoint,
        BinarySearch: BinarySearch,
        LinearSolve2: LinearSolve2,
        Rotate: Rotate,
        SegmentIntersect: SegmentIntersect,
        HyperbolicCenter: HyperbolicCenter,
        TangentVector: TangentVector,
        LineIntersect: LineIntersect,
      }
    })()
    )};