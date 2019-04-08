(geom) => {return(
    function HyperbolicTesselation(k, n_levels) {
        // We are expecting parameters such that 1/a + 1/b + 1/c = 1;
        var a = k[0], b = k[1], c = k[2];
        var k = k;
        
        var vertices = [];
        var faces = [];
        
        function GetTriangleLengths() {
          var gamma = Math.PI / c;
          var v = [Math.cos(gamma), Math.sin(gamma)];
          var angle = -Math.PI * (1- (1/b) - (1/c));
          var d = [Math.cos(angle), Math.sin(angle)];
          function pred(t) {
            var line = geom.PoincareLineDir(geom.Scale(v,t), d);
            if (line.c[1] - line.r > 0) return 1;
            return - Math.cos(Math.PI / a) + line.c[1] / line.r;
          }
          var t = geom.BinarySearch(pred, 0.000001, 0, 1);
          var p = geom.Scale(v,t);
          var l = geom.PoincareLineDir(p, d);
          var q =  [l.c[0] - Math.sqrt((l.r**2) - (l.c[1]**2)),0];
          var z = [0,0];
          return [geom.PoincareDistance(z,p),
                  geom.PoincareDistance(q,z),
                  geom.PoincareDistance(p,q)];
        }
        
        var lenghts = GetTriangleLengths();
        
        function NewVertex(index, i, pos, base_dir) {
          return {
            index: index,
            base_dir: base_dir,
            i : i,
            pos: pos,
            f: new Array(2 * k[i]).fill(null),
            v: new Array(2 * k[i]).fill(null),
            level: 0,
          }
        }
        
        function NewFace(a,b,c) {
          return {
            v: [a,b,c],
            prev: null,
          }
        }
        
        // We add the first face with its vertices to the
        // tesselation and start expanding from there.
        
        function InitializeTesselation() {
          var base_dir = [1,0];
          var pos = [0,0];
          
          for (var i = 0; i < 3; ++i) {
            var p = NewVertex(vertices.length, i, pos, base_dir);
            var lp = geom.PoincareSegmentDir(pos, base_dir, lenghts[(6-(2*i+1))%3]);
            pos = [lp.q[0], lp.q[1]];
            base_dir = geom.Rotate(lp.tq, Math.PI * (1- (1/k[(i+1)%3])));
            vertices.push(p);
          }
          var f = NewFace(vertices[0], vertices[1], vertices[2]);
          for (var i = 0; i < 3; ++i) {
            vertices[i].f[0] = f;
            vertices[i].v[0] = vertices[(i+1)%3];
            vertices[i].v[1] = vertices[(i+2)%3];
          }
          vertices[1].level = 1;
          vertices[2].level = 1;
          faces.push(f);
        }
        
        InitializeTesselation();
        
        function FindIndex(w,v) {
          var wj = 0;
          for (wj = 0; wj < 2*k[w.i]; ++wj) {
            if (w.v[wj] == v) break;
          }
          return wj;
        }
        
        function ProcessVertex(v, max_faces = Infinity ) {
          var nfaces = 2 * k[v.i];
          for (var j = 0; j < Math.min(nfaces, max_faces); ++j) {
            if (v.f[j] == null && v.f[(j+1)%nfaces] == null) {
              var w = v.v[j];
              var wj = FindIndex(w,v);
              wj = (2*k[w.i] + wj - 1) % (2*k[w.i]);
      
              var u_i = 3 - v.i - w.i;
              var l = geom.PoincareSegmentDir(v.pos,
                                              geom.Rotate(v.base_dir, Math.PI * (j+1) / k[v.i] ),
                                              lenghts[w.i]);
              var u_base_dir = geom.Scale(l.tq, -1);
              var u_pos = [l.q[0], l.q[1]];
              var u = NewVertex(vertices.length, u_i, u_pos, u_base_dir);
              u.level = v.level + 1;
              u.v[0] = v;
              u.v[1] = w;
              v.v[(j+1)%nfaces] = u;
              w.v[wj] = u;
              var f = NewFace(v, w, u);
              f.prev = v.f[j-1];
              u.f[0] = f;
              v.f[j] = f;
              w.f[wj] = f;
              faces.push(f);
              vertices.push(u);
            } else if (v.f[j] == null && v.f[(j+1)%nfaces] != null) {
              var w = v.v[j];
              var u = v.v[(j+1)% nfaces ];
              var wv_j = FindIndex(w,v);
              var f = NewFace(v, w, u);
              f.prev = v.f[j-1];
              v.f[j] = f;
              w.v[(2*k[w.i] + wv_j-1 ) % (2*k[w.i])] = u;
              w.f[(2*k[w.i] + wv_j-1 ) % (2*k[w.i])] = f;
              var uv_j = FindIndex(u,v);
              u.v[(uv_j+1) % (2*k[u.i])] = w;
              u.f[uv_j] = f;
              faces.push(f);
            }
          }
        }
        
        var j = 0;
        for (var l = 0; l < n_levels; ++l) {
          var lvertices = [];
          for (; j < vertices.length && vertices[j].level == l; ++j) {
            lvertices.push(vertices[j]);
          }
          lvertices.sort(function(u,w) { return geom.Norm(u.pos) - geom.Norm(w.pos); });
          lvertices.forEach(function(v) {
            ProcessVertex(v);
          });
        }
        
        faces.forEach(function(f) {
          f.center = geom.HyperbolicCenter(f.v[0].pos, f.v[1].pos, f.v[2].pos);
        });
       
        return {
          vertices: vertices,
          faces: faces,
          k: k,
        }
      }

    )};