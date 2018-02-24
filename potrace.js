/*
potrace.js 私家版．2017/12/11版
文字列の整理，notation修正．even領域とodd領域の判別を可能とする．パス精度の追加．
パス文字列をsvgに準拠．パスの精度桁数の追加．枠線処理を内包．処理速度の向上．バグの修正等．
Node.js環境での動作
by defghi1977
TODO:無駄な部分を削除する
TODO:パフォーマンスのさらなる改善
http://defghi1977-onblog.blogspot.jp/2012/06/svg_28.html
移植元
http://www.libspark.org/svn/as3/PotrAs/src/com/nitoyon/potras/
javascript版製作者
http://d.hatena.ne.jp/project_the_tower2/20110724/1311473645
potrace オリジナル 
Copyright © 2001-2011 Peter Selinger. 
http://potrace.sourceforge.net/
Lisence GPL2
http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
*/
Potrace = (function(){
	"use strict";
	//string literal
	var NULL_STR = "";
	var SPACE = " ";
	var COMMA = ",";
	var DOT = ".";
	var PLUS = "+";
	var MINUS = "-";
	var LF = "\n";
	var CLOSED_PATH = "ClosedPath";
	var CANVAS = "canvas";
	var M = "M";
	var L = "L";
	var C = "C";
	var _m = "m";
	var _l = "l";
	var _c = "c";
	var Z = "Z";
	var SVG_NS = "http://www.w3.org/2000/svg";
	var PATH = "path";
	var G = "g";
	var D = "d";
	var CLASS_ID = "class_id";
	var CLASS = "class";
	var EVEN = "even";
	var ODD = "odd";
	var MARGIN = 1;

	//////////////////////////////////////////////////////////
	//define
	//////////////////////////////////////////////////////////
	var INFTY           = 10000000;
	var POTRACE_CURVETO = 1;
	var POTRACE_CORNER  = 2;
	//////////////////////////////////////////////////////////
	//private
	//////////////////////////////////////////////////////////
	var param = {
		threshold:0.625,
		turdSize:3,
		turnPolicy:1,
		alphamax:0.9,
		precision:-1,
		noFilter:false,//set true when source image already binalized.
		isRelative:false//set true to output svg path string via relative commands.
	};
	var pathList = (function(){
		//XOR　反転　行列修正
		var filter = new Float64Array([
			-1, 0, 0, 0, 255,
			0 ,-1, 0 ,0 ,255,
			0 , 0,-1 ,0 ,255,
			0 , 0, 0, 1, 0
		]);
		return function(oriBmd,param){
			var pathList = [];
			var y = 0;
			var cpyBmd = oriBmd.clone();
			var point = new Point();
			while(findNext(cpyBmd, point)){
				// calculate the sign by looking at the original
				let sign = isBlack(oriBmd.getPixel(point.x, point.y)) ? PLUS : MINUS;
				// calculate the path
				let p = findPath(cpyBmd, new Point(point.x, point.y - 1), sign, param.turnPolicy);
				if(!p){
					pathList = null;
					break;
				}
				//update buffered image
				xorPath(cpyBmd, p, filter);
				// if it's a turd, eliminate it, else append it to the list
				if(p.area > param.turdSize){
					p.isEven = sign == MINUS;//by def
					pathList.push(p);
				}
			}
			return pathList;
		};
	})();
	var findNext = function(bmd,pt){
		for(let y = pt.y, height = bmd.height; y < height; y=(y+1)|0){
			for(let x=0, width = bmd.width; x < width ;x=(x+1)|0){
				if(isBlack(bmd.getPixel(x,y))){
					pt.x = x;
					pt.y = y;
					return true;
				}
			}
		}
		return false;
	};
	var isBlack = function(pixel){
		return (pixel.r + pixel.g + pixel.b) == 0;
	};
	var findPath = (function(){
		var rotateRight = new Float64Array([0, -1, 1, 0]);
		var rotateLeft = new Float64Array([0, 1, -1, 0]);
		return function(bmd, startPoint, sign, turnpolicy){
			var area = 0;
			var pointList = [];
			var pt = startPoint.clone();
			var dir = new Point(0, 1);

			while(true){
				pointList.push(pt.clone());
			
				//move to next point
				pt.offset(dir.x,dir.y);
				area += pt.x * (-1) * dir.y;
			
				//path complete?
				if(pt.equals(startPoint)){
					break;
				}
			
				//determine next direction
				let c =isBlack(bmd.getPixel(pt.x + (dir.x - dir.y - 1)/2, pt.y + (dir.y + dir.x + 1)/2));
				let d =isBlack(bmd.getPixel(pt.x + (dir.x + dir.y - 1)/2, pt.y + (dir.y - dir.x + 1)/2));
			
				//c d
				//X _
				if(c && !d){ //ambiguous turn
					//if(true){ //TODO
						dir.transformPoint(rotateLeft);
					//}else{
					//	dir.transformPoint(rotateRight);
					//}
				}else if(c){ //left turn
					dir.transformPoint(rotateLeft);
				}else if(!d){ // right turn 
					dir.transformPoint(rotateRight);
				}
			}
		
			//allocate new path object
			var path = {};
			path.priv = pointList;
			path.area = area;
			path.sign = sign;

			return path;
		};
	})();


	var xorPath = function(bm, p, filter){
		var priv = p.priv;
		var len = priv.length;
		var min = Math.min;
		var max = Math.max;
		//get minimum x
		var minX = 99999;
		for(let i=0; i<len; i=(i+1)|0){
			minX = +min(minX, priv[i].x);
		}
		var y1 = priv[len - 1].y;
		var pt = new Point();

		var rect = new Rectangle();
		for(let i=0; i<len; i=(i+1)|0){
			let x = priv[i].x;
			let y = priv[i].y;

			if(y != y1){
				let y2 = max(y, y1);
				pt.x = minX; pt.y = y2;
				rect.x = minX;
				rect.y = y2;
				rect.width = x - minX;
				rect.height = 1;
				bm.applyFilter(rect, pt, filter);
				y1 = y;
			}
		}
	};
	var processPath = function(pList, param){
		var closedPathList = new ClosedPathList();
		var outpaths = [];
		for(var i=0, len=pList.length; i<len; i = (i+1)|0){
			let pList_i = pList[i];
			let sums = calcSums(pList_i.priv);
			let lons = calcLon(pList_i.priv);
			let po = bestPolygon(pList_i.priv, lons, sums);
			let vertex = adjustVertices(pList_i.priv, sums, po);
			let closedPath = smooth(vertex, pList_i.sign, param.alphamax);
			adjust(closedPath, param.precision);//by def
			closedPath.isEven = pList_i.isEven;//by def
			closedPathList.$a.push(closedPath);//後々のループを省くためにここで生成してしまう．
			if(closedPath.hasOwnProperty(CLASS_ID) && closedPath.class_id == CLOSED_PATH){
				outpaths.push(closedPath.outpath());
			}
		}
		closedPathList._outpath = outpaths.join(NULL_STR);
		return closedPathList;
	};
	var getSvgElement = function(tagName){
		return document.createElementNS(SVG_NS, tagName);
	};
	var getPathElement = function(){
		return getSvgElement(PATH);
	};
	var getGElement = function(){
		return getSvgElement(G);
	};
	/**
	 * Preparation: fill in the sum* fields of a path (used for later
	 * rapid summing). Return 0 on success, 1 with errno set on
	 * failure.
	 */
	var calcSums = function(pt){
		var n = pt.length;

		var sums = new Array(n + 1);

		// origin
		var x0 = pt[0].x;
		var y0 = pt[0].y;

		// preparatory computation for later fast summing
		sums[0] = new Sums();
		var sums_0 = sums[0];
		sums_0.x2 = sums_0.xy = sums_0.y2 = sums_0.x = sums_0.y = 0;
		for(let i = 0; i < n; i=(i+1)|0){
			let x = pt[i].x - x0;
			let y = pt[i].y - y0;
			sums[i + 1] = new Sums();
			let sums_i = sums[i];
			let sums_i_1 = sums[i + 1];
			sums_i_1.x  = sums_i.x  + x;
			sums_i_1.y  = sums_i.y  + y;
			sums_i_1.x2 = sums_i.x2 + x*x;
			sums_i_1.xy = sums_i.xy + x*y;
			sums_i_1.y2 = sums_i.y2 + y*y;
		}
		return sums;
	};
	/**
	 * Stage 1: determine the straight subpaths (Sec. 2.2.1). Fill in the
	 * "lon" component of a path object (based on pt/len).        For each i,
	 * lon[i] is the furthest index such that a straight line can be drawn
	 * from i to lon[i]. Return 1 on error with errno set, else 0.
	 */
	var calcLon=function(pt){
		var pti,ptk,ptk1,ptmod;
		var lon= [];
		var n = pt.length;
		var j;
		var abs = Math.abs;
		var min = Math.min;
		// initialize the nc data structure. Point from each point to the
		// furthest future point to which it is connected by a vertical or
		// horizontal segment. We take advantage of the fact that there is
		// always a direction change at 0 (due to the path decomposition
		// algorithm). But even if this were not so, there is no harm, as
		// in practice, correctness does not depend on the word "furthest"
		// above.
		var nc = [];
		var k = 0;
		for(let i = n - 1; i >= 0; i=(i-1)|0){
			pti = pt[i];
			ptk = pt[k];
			if(pti.x != ptk.x && pti.y != ptk.y){
				k = i + 1; // necessarily i<n-1 in this case
			}
			nc[i] = k;
		}

		// determine pivot points: for each i, let pivk[i] be the furthest k
		// such that all j with i<j<k lie on a line connecting i,k. */
		var pivk = [];
		var ct;
		var dir;
		var k1;
		var constraint0 = new Point();
		var constraint1 = new Point();
		var cur = new Point();
		var off = new Point();
		var dk = new Point();
		for(let i = n - 1; i >= 0; i=(i-1)|0){
			pti = pt[i];
			ptmod = pt[mod(i + 1, n)];
			ct = new Int8Array([0, 0, 0, 0]);

			// keep track of "directions" that have occurred
			dir = (3 + 3 * (ptmod.x - pti.x) + (ptmod.y - pti.y)) / 2;
			ct[dir]++;
			constraint0.x = constraint0.y = constraint1.x = constraint1.y = 0;

			// find the next k such that no straight line from i to k
			k = nc[i];
			k1 = i;

			let foundk = false;
			while(true){
				ptk = pt[k];
				ptk1 = pt[k1];
				
				dir = (3 + 3 * sign(ptk.x - ptk1.x) + sign(ptk.y - ptk1.y)) / 2;
				ct[dir]++;

				// if all four "directions" have occurred, cut this path
				if(ct[0] && ct[1] && ct[2] && ct[3]){
					pivk[i] = k1;
					foundk = true;
					break;
				}

				cur.x = ptk.x - pti.x;
				cur.y = ptk.y - pti.y;

				// see if current constraint is violated
				if(xprod(constraint0, cur) > 0 || xprod(constraint1, cur) < 0){
					break;
				}

				// else, update constraint
				if (abs(cur.x) <= 1 && abs(cur.y) <= 1){
					// no constraint
				}else{
					off.x = cur.x + ((-cur.y >= 0 && (-cur.y > 0 || cur.x < 0)) ? 1 : -1);
					off.y = cur.y + ((-cur.x <= 0 && (-cur.x < 0 || cur.y < 0)) ? 1 : -1);
					if(xprod(constraint0, off) <= 0){
						constraint0.x = off.x;
						constraint0.y = off.y;
					}

					off.x = cur.x + ((-cur.y <= 0 && (-cur.y < 0 || cur.x < 0)) ? 1 : -1);
					off.y = cur.y + ((-cur.x >= 0 && (-cur.x > 0 || cur.y < 0)) ? 1 : -1);
					if(xprod(constraint1, off) >= 0){
						constraint1.x = off.x;
						constraint1.y = off.y;
					}
				}

				k1 = k;
				k = nc[k1];
				if(!cyclic(k, i, k1)){
					break;
				}
			}

			// constraint
			if(!foundk){
				pti = pt[i];
				ptk = pt[k];
				ptk1 = pt[k1];
				dk.x = sign(ptk.x - ptk1.x);
				dk.y = sign(ptk.y - ptk1.y);
				cur.x = ptk1.x - pti.x;
				cur.y = ptk1.y - pti.y;

				// find largest integer j such that xprod(constraint0, cur+j*dk) >= 0 
				// and xprod(constraint1, cur+j*dk) <= 0. Use bilinearity of xprod.
				var a = xprod(constraint0, cur);
				var b = xprod(constraint0, dk);
				var c = xprod(constraint1, cur);
				var d = xprod(constraint1, dk);
				// find largest integer j such that a+j*b>=0 and c+j*d<=0.
				// This can be solved with integer arithmetic.
				j = INFTY;
				if(b > 0){
					j = floordiv(-a, b);
				}
				if(d < 0){
					j = min(j, floordiv(c, -d));
				}

				pivk[i] = mod(k1 + j, n);
			}

			// foundk
		} // for i

		// clean up: for each i, let lon[i] be the largest k such that for
		// all i' with i<=i'<k, i'<k<=pivk[i'].
		j = pivk[n - 1];
		lon[n - 1] = j;
		for(let i = n - 2; i >= 0; i=(i-1)|0){
			if(cyclic(i + 1, pivk[i], j)){
				j = pivk[i];
			}
			lon[i] = j;
		}

		for(let i = n - 1; cyclic(mod(i + 1, n), j, lon[i]); i=(i-1)|0){
			lon[i] = j;
		}

		return lon;
	};
	/**
	 * Stage 2: calculate the optimal polygon (Sec. 2.2.2-2.2.4).
	 *
	 * find the optimal polygon. Fill in the m and po components. Return 1
	 * on failure with errno set, else 0. Non-cyclic version: assumes i=0
	 * is in the polygon. Fixme: ### implement cyclic version.
	 */
	var bestPolygon = function(pt, lon, sums){
		var i, j, k;
		var n = pt.length;

		var clip0 = new Float64Array(n);		// clip0[n]: longest segment pointer, non-cyclic
		var clip1 = new Float64Array(n + 1);	// clip1[n+1]: backwards segment pointer, non-cyclic

		// calculate clipped paths
		for(i = 0; i < n; i=(i+1)|0){
			var c = mod(lon[mod(i - 1, n)] - 1, n);
			if(c == i){
				c = mod(i + 1, n);
			}
			if(c < i){
				clip0[i] = n;
			}else{
				clip0[i] = c;
			}
		}

		// calculate backwards path clipping, non-cyclic. j <= clip0[i] iff clip1[j] <= i, for i,j=0..n.
		j = 1;
		for(i = 0; i < n; i=(i+1)|0){
			while(j <= clip0[i]){
				clip1[j] = i;
				j=(j+1)|0;
			}
		}
		var seg0 = new Array(n + 1);		// seg0[m+1]: forward segment bounds, m<=n
		var seg1 = new Array(n + 1);	 	// seg1[m+1]: backward segment bounds, m<=n

		// calculate seg0[j] = longest path from 0 with j segments
		i = 0;
		for(j = 0; i < n; j=(j+1)|0){
			seg0[j] = i;
			i = clip0[i];
		}
		seg0[j] = n;
		var m = j;

		// calculate seg1[j] = longest path to n with m-j segments
		i = n;
		for(j = m; j > 0; j=(j-1)|0){
			seg1[j] = i;
			i = clip1[i];
		}
		seg1[0] = 0;

		// now find the shortest path with m segments, based on penalty3
		// note: the outer 2 loops jointly have at most n interations, thus
		// the worst-case behavior here is quadratic. In practice, it is
		// close to linear since the inner loop tends to be short.
		var pen   = new Array(n + 1); // pen[n+1]: penalty vector
		var prev  = new Array(n + 1);	// prev[n+1]: best path pointer vector
		var thispen;
		pen[0] = 0;
		for(j = 1; j <= m; j=(j+1)|0){
			for(i = seg1[j]; i <= seg0[j]; i=(i+1)|0){
				let best = -1;
				for(k = seg0[j - 1]; k >= clip1[i]; k=(k-1)|0){
					thispen = penalty3(pt, k, i, sums) + pen[k];
					if(best < 0 || thispen < best){
						prev[i] = k;
						best = thispen;
					}
				}
				pen[i] = best;
			}
		}

		// read off shortest path
		var po = new Array(m);
		for(i = n, j = m - 1; i > 0; j=(j-1)|0){
			po[j] = i = prev[i];
		}

		return po;
	};
	/**
	 * Stage 3: vertex adjustment (Sec. 2.3.1).
	 *
	 * Adjust vertices of optimal polygon: calculate the intersection of
	 * the two "optimal" line segments, then move it into the unit square
	 * if it lies outside. Return 1 with errno set on error; 0 on
	 * success.
	 */
	var adjustVertices= function(pt, sums, po){
		var m = po.length;
		var n = pt.length;
		var abs = Math.abs;
		var i, j, k, l;

		// represent each line segment as a singular quadratic form; the
		// distance of a point (x,y) from the line segment will be
		// (x,y,1)Q(x,y,1)^t, where Q=q[i].
		var q = [];
		var v = new Point3D(); 
		var ctr = new Point();
		var dir = new Point();
		var $v = v.$v;
		for(i = 0; i < m; i=(i+1)|0){
			// calculate "optimal" point-slope representation for each line segment
			j = po[(i + 1) % m];
			j = mod(j - po[i], n) + po[i];
			ctr.x = ctr.y = dir.x = dir.y = 0;
			pointslope(pt, sums, po[i], j, ctr, dir);

			q[i] = new QuadraticForm();
			var d1 = dir.x * dir.x + dir.y * dir.y;
			if(d1 != 0.0){
				$v[0] = dir.y;
				$v[1] = -dir.x;
				$v[2] = -$v[1] * ctr.y - $v[0] * ctr.x;
				q[i].fromVectorMultiply(v).scalar(1 / d1);
			}
		}

		// now calculate the "intersections" of consecutive segments.
		// Instead of using the actual intersection, we find the point
		// within a given unit square which minimizes the square distance to
		// the two lines.
		var vertex = [];
		var p0 = pt[0].clone();
		var s = new Point();
		var w = new Point();
		var _q = new QuadraticForm();
		var minCoord = new Point;
		for(i = 0; i < m; i=(i+1)|0){
			let z;
			vertex[i] = new Point();
			let vertexi = vertex[i];

			// let s be the vertex, in coordinates relative to x0/y0
			s.x = pt[po[i]].x - p0.x;
			s.y = pt[po[i]].y - p0.y;

			// intersect segments i-1 and i
			j = mod(i - 1, m);

			// add quadratic forms
			let Q = q[j].clone().add(q[i]);
			let Q$m = Q.$m;
			let Q$m0 = Q$m[0], Q$m1 = Q$m[1];
			let Q$m0$v = Q$m0.$v, Q$m1$v = Q$m1.$v;
			let v$v = v.$v;
			while(true){
				// minimize the quadratic form Q on the unit square
				// find intersection
				let det = Q$m0$v[0] * Q$m1$v[1] - Q$m0[1] * Q$m1[0];
				if(det != 0.0){
					w.x = (-Q$m0$v[2] * Q$m1$v[1] + Q$m1$v[2] * Q$m0$v[1]) / det;
					w.y = ( Q$m0$v[2] * Q$m1$v[0] - Q$m1$v[2] * Q$m0$v[0]) / det;
					break;
				}

				// matrix is singular - lines are parallel. Add another,
				// orthogonal axis, through the center of the unit square
				if(Q$m0[0] > Q$m1$v[1]){
					v$v[0] = -Q$m0$v[1];
					v$v[1] = Q$m0$v[0];
				}else if(Q$m1$v[1]){
					v$v[0] = -Q$m1$v[1];
					v$v[1] = Q$m1$v[0];
				}else{
					v$v[0] = 1;
					v$v[1] = 0;
				}

				let d = v$v[0] * v$v[0] + v$v[1] * v$v[1];
				v$v[2] = - v$v[1] * s.y - v$v[0] * s.x;
				Q.add(_q.fromVectorMultiply(v)).scalar(1 / d);
			}

			let dx = abs(w.x - s.x);
			let dy = abs(w.y - s.y);
			if(dx <= 0.5 && dy <= 0.5){
				vertexi.x = w.x + p0.x;
				vertexi.y = w.y + p0.y;
				continue;
			}

			// the minimum was not in the unit square; now minimize quadratic
			// on boundary of square
			let min, cand; // minimum and candidate for minimum of quad. form
			minCoord.x = s.x; minCoord.y = s.y; // coordinates of minimum
			min = Q.apply(s);

			if(Q$m0$v[0] != 0.0){
				for(z = 0; z < 2; z=(z+1)|0){  // value of the y-coordinate
					w.y = s.y - 0.5 + z;
					w.x = - (Q$m0$v[1] * w.y + Q$m0$v[2]) / Q$m0$v[0];
					dx = (w.x - s.x > 0 ? w.x - s.x : s.x - w.x);
					cand = Q.apply(w);
					if(dx <= 0.5 && cand < min){
						min = cand;
						minCoord.x = w.x; minCoord.y = w.y;
					}
				}
			}

			if(Q$m1$v[1] != 0.0){
				for(z = 0; z < 2; z=(z+1)|0){ // value of the x-coordinate
					w.x = s.x - 0.5 + z;
					w.y = - (Q$m1$v[0] * w.x + Q$m1$v[2]) / Q$m1$v[1];
					dy = (w.y - s.y > 0 ? w.y - s.y : s.y - w.y);
					cand = Q.apply(w);
					if(dy <= 0.5 && cand < min){
						min = cand;
						minCoord.x = w.x; minCoord.y = w.y;
					}
				}
			}

			// check four corners
			for(l = 0; l < 2; l=(l+1)|0){
				for(k = 0; k < 2; k=(k+1)|0){
					w.x = s.x - 0.5 + l;
					w.y = s.y - 0.5 + k;
					cand = Q.apply(w);
					if(cand < min){
						min = cand;
						minCoord.x = w.x; minCoord.y = w.y;
					}
				}
			}

			vertexi.x = minCoord.x + p0.x;
			vertexi.y = minCoord.y + p0.y;
			//continue;
		}

		return vertex;
	};
	/**
	 *  Stage 4: smoothing and corner analysis (Sec. 2.3.3)
	 *
	 *  Always succeeds and returns 0
	 */
	var smooth = function(vertex, sign, alphamax){

		if(!alphamax){ alphamax =1.0;}
		var m = vertex.length;
		var closedPath = new ClosedPath();

		if(sign == MINUS){
			// reverse orientation of negative paths
			let tmp;
			for(let i = 0, j = m - 1; i < j; i=(i+1)|0, j=(j-1)|0){
				tmp = vertex[i];
				vertex[i] = vertex[j];
				vertex[j] = tmp;
			}
		}

		// examine each vertex and find its best fit
		var p2 = new Point();
		var p3 = new Point();
		var p4 = new Point();
		for(let i = 0; i < m; i=(i+1)|0){
			let j = (i + 1) % m;
			let k = (i + 2) % m;

			let vertexi = vertex[i];
			let vertexj = vertex[j];
			let vertexk = vertex[k];

			interval(1 / 2.0, vertexk, vertexj, p4);

			let curve = new Curve();
			curve.vertex.x = vertexj.x;
			curve.vertex.y = vertexj.y;

			let denom = ddenom(vertexi, vertexk);
			let alpha;
			if(denom != 0.0){
				let dd = dpara(vertexi, vertexj, vertexk) / denom;
				dd = dd < 0 ? -dd : dd;
				alpha = dd > 1 ? (1 - 1.0 / dd) / 0.75 : 0;
			}else{
				alpha = 4 / 3.0;
			}
			curve.alpha0 = alpha;   // remember "original" value of alpha
			
			let c = curve.c;
			let c0 = c[0], c1=c[1], c2=c[2];
			if(alpha > alphamax){ // pointed corner
				curve.tag  = POTRACE_CORNER;
				c0.x = 0; c0.y = 0;
				c1.x = vertexj.x; c1.y = vertexj.y;
				c2.x = p4.x; c2.y = p4.y;
			}else{
				if(alpha < 0.55){
					alpha = 0.55;
				}else if(alpha > 1){
					alpha = 1;
				}
				interval(0.5 + 0.5 * alpha, vertexi, vertexj, p2);
				interval(0.5 + 0.5 * alpha, vertexk, vertexj, p3);
				curve.tag = POTRACE_CURVETO;
				c0.x = p2.x; c0.y = p2.y;
				c1.x = p3.x; c1.y = p3.y;
				c2.x = p4.x; c2.y = p4.y;
			}
			curve.alpha = alpha; // store the "cropped" value of alpha
			curve.beta  = 0.5;

			closedPath.$a[j] = curve;
		}

		return closedPath;
	};

	/**
	 *  Stage 5: adjusting the coordinates of the points.
	 * 
	 *  round to precision. by defghi1977
	 *  Always succeeds
	 */
	var adjust = function(closedPath, precision){
		
		//path出力の精度桁数
		var round;
		if(precision == undefined || precision < 0){
			round = function(value){return value;};
		}else if(precision == 0){
			round = Math.round;
		}else{
			var dec = Math.pow(10, Math.min(15, precision));
			round = (function(dec, r){
				return function(value){return r(value * dec)/dec;}
			})(dec, Math.round);
		}
		var curves = closedPath.$a;
		for(let i = 0, len = curves.length; i<len; i=(i+1)|0){
			let points = curves[i].c;
			let point;
			point = points[0];
			point.x = round(point.x - MARGIN);
			point.y = round(point.y - MARGIN + 1);
			point = points[1];
			point.x = round(point.x - MARGIN);
			point.y = round(point.y - MARGIN + 1);
			point = points[2];
			point.x = round(point.x - MARGIN);
			point.y = round(point.y - MARGIN + 1);
		}
		return closedPath;
	};

	/////////////////////////////////////////
	// dpara function
	/**
	 *  return (p1-p0)x(p2-p0), the area of the parallelogram
	 */
	/////////////////////////////////////////
	var dpara = function(p0, p1, p2){
		var x1 = p1.x - p0.x;
		var y1 = p1.y - p0.y;
		var x2 = p2.x - p0.x;
		var y2 = p2.y - p0.y;

		return x1 * y2 - x2 * y1;
	};
	/////////////////////////////////////////
	// dorth_infty function
	/**
	 *  return a direction that is 90 degrees counterclockwise from p2-p0,
	 *  but then restricted to one of the major wind directions (n, nw, w, etc)
	 */
	/////////////////////////////////////////
	var dorth_infty=function(p0, p2){
		return new Point(
			sign(p2.x  - p0.x), 
			-sign(p2.y - p0.y)
		);
	};
	/////////////////////////////////////////
	// ddenom function
	/**
	 *  ddenom/dpara have the property that the square of radius 1 centered
	 *  at p1 intersects the line p0p2 iff |dpara(p0,p1,p2)| <= ddenom(p0,p2)
	 */
	/////////////////////////////////////////
	var ddenom = function(p0, p2){
		var r = dorth_infty(p0, p2);
		return r.y * (p2.x - p0.x) - r.x * (p2.y - p0.y);
	};
	/////////////////////////////////////////
	// interval function
	/**
	 *  range over the straight line segment [a,b] when lambda ranges over [0,1]
	 */
	/////////////////////////////////////////
	var interval = function(lambda, a, b, ret){
		ret.x = a.x + lambda * (b.x - a.x);
		ret.y = a.y + lambda * (b.y - a.y);
	};
	/////////////////////////////////////////
	// penalty3 function
	/////////////////////////////////////////
	var pointslope=function(pt, sums, i, j, ctr, dir){
		// assume i<j

		var n = pt.length;

		var x, y, x2, xy, y2;
		var k;
		var a, b, c, lambda2, l;
		var r = 0; // rotations from i to j

		while (j >= n){
			j -= n;
			r += 1;
		}
		while(i >= n){
			i -= n;
			r -= 1;
		}
		while(j < 0){
			j += n;
			r -= 1;
		}
		while(i < 0){
			i += n;
			r += 1;
		}
		
		var sums_j_1 = sums[j + 1];
		var sums_i = sums[i];
		var sums_n = sums[n];
		x  = sums_j_1.x  - sums_i.x  + r * sums_n.x;
		y  = sums_j_1.y  - sums_i.y  + r * sums_n.y;
		x2 = sums_j_1.x2 - sums_i.x2 + r * sums_n.x2;
		xy = sums_j_1.xy - sums_i.xy + r * sums_n.xy;
		y2 = sums_j_1.y2 - sums_i.y2 + r * sums_n.y2;
		k = j + 1 - i + r * n;
		
		ctr.x = x / k;
		ctr.y = y / k;

		a = (x2 - x * x / k) / k;
		b = (xy - x * y / k) / k;
		c = (y2 - y * y / k) / k;
		
		var sqrt = Math.sqrt;
		var abs = Math.abs;
		lambda2 = (a + c + sqrt((a - c)*(a - c) + 4 * b * b)) / 2; // larger e.value 

		// now find e.vector for lambda2
		a -= lambda2;
		c -= lambda2;

		if(abs(a) >= abs(c)){
			l = sqrt(a * a + b * b);
			if(l != 0){
				dir.x = -b / l;
				dir.y =  a / l;
			}
		}else{
			l = sqrt(c * c + b * b);
			if (l != 0){
				dir.x = -c / l;
				dir.y =  b / l;
			}
		}
		if(l == 0){
			dir.x = dir.y = 0;  // sometimes this can happen when k=4:
			                      // the two eigenvalues coincide
		}
	};
	/////////////////////////////////////////
	// penalty3 function
	/////////////////////////////////////////
	var penalty3 = function(pt, i, j, sums){
		var n = pt.length;

		// assume 0 <= i < j <= n
		var r = 0; // rotations from i to j

		if(j >= n){
			j -= n;
			r += 1;
		}
		var sums_j1 = sums[j + 1];
		var sums_i = sums[i];
		var sums_n = sums[n];
		var x  = sums_j1.x  - sums_i.x  + r * sums_n.x;
		var y  = sums_j1.y  - sums_i.y  + r * sums_n.y;
		var x2 = sums_j1.x2 - sums_i.x2 + r * sums_n.x2;
		var xy = sums_j1.xy - sums_i.xy + r * sums_n.xy;
		var y2 = sums_j1.y2 - sums_i.y2 + r * sums_n.y2;
		var k  = j + 1 - i + r * n;
		
		var pt_i = pt[i];
		var pt_j = pt[j];
		var pt_0 = pt[0];
		var px = (pt_i.x + pt_j.x) / 2.0 - pt_0.x;
		var py = (pt_i.y + pt_j.y) / 2.0 - pt_0.y;
		var ey = (pt_j.x - pt_i.x);
		var ex =  -(pt_j.y - pt_i.y);

		var a = ((x2 - 2 * x * px) / k + px * px);
		var b = ((xy - x * py - y * px) / k + px * py);
		var c = ((y2 - 2 * y * py) / k + py * py);
		
		var s = ex * ex * a  +  2 * ex * ey * b  +  ey * ey * c;

		return Math.sqrt(s);
	};
	/////////////////////////////////////////
	// mod function
	// floordiv function
	/**
	 *  some useful macros. Note: the "mod" macro works correctly for
	 *  negative a. Also note that the test for a>=n, while redundant,
	 *  speeds up the mod function by 70% in the average case (significant
	 *  since the program spends about 16% of its time here - or 40%
	 *  without the test). The "floordiv" macro returns the largest integer
	 *  <= a/n, and again this works correctly for negative a, as long as
	 *  a,n are integers and n>0.
	 */
	////////////////////////////////////////
	var mod=function(a,n){
		return a >= n ? a % n : a >= 0 ? a : n - 1 - (-1 - a) % n;
	};
	var floordiv=function(a, n){
		return a >= 0 ? Math.floor(a / n) : Math.floor(-1 - (-1 - a) / n);
	};
	/////////////////////////////////////////
	// sign function
	////////////////////////////////////////
	var sign = function(x){
		return (x > 0 ? 1 : x < 0 ? -1 : 0);
	};
	/////////////////////////////////////////
	// xprod function
	// calculate p1 x p2
	////////////////////////////////////////
	var xprod=function(p1, p2){
		return p1.x * p2.y - p1.y * p2.x;
	};
	/////////////////////////////////////////
	// cyclic  function
	// return 1 if a <= b < c < a, in a cyclic sense (mod n)
	////////////////////////////////////////
	var cyclic = function(a, b, c){
		if (a <= c){
			return (a <= b) && (b < c);
		}else{
			return (a <= b) || (b < c);
		}
	};
	/////////////////////////////////////////
	// Curve class
	////////////////////////////////////////
	var Curve = function(){
		this.initialize.apply(this, arguments);
	};
	Curve.prototype = {
		/**
		 *  Constructor.
		 *
		 *  initialize the members of the given curve structure to size m.
		 *  Return 0 on success, 1 on error with errno set.
		 */
		initialize:function(){
			this.c = new Array(3);
			this.c[0] = new Point();
			this.c[1] = new Point();
			this.c[2] = new Point();
			this.vertex = new Point();
		},
		toString:function(){
			return ["alpha0: " , this.alpha0 , LF
			     , "alpha:  " , this.alpha , LF
			     , "beta:   " , this.beta , LF
			     , "corner: " , (this.tag == POTRACE_CORNER) , LF
			     , "bezier: " , this.c[0] , COMMA , this.c[1] , COMMA , this.c[2]].join(NULL_STR);
		},
		tag:0, /* tag[n]: POTRACE_CORNER or POTRACE_CURVETO */
		c:null,
		/* c[n][i]: control points. c[n][0] is unused for tag[n]=POTRACE_CORNER */
		/* the remainder of this structure is special to privcurve, and is
		   used in EPS debug output and special EPS "short coding". These
		   fields are valid only if "alphacurve" is set. */
		vertex:null, /* for POTRACE_CORNER, this equals c[1] */
		alpha:0,     /* only for POTRACE_CURVETO */
		alpha0:0,    /* "uncropped" alpha parameter - for debug output only */
		beta:0
	};

	/////////////////////////////////////////
	// ClosedPath class
	////////////////////////////////////////
	var ClosedPath = function(){
		this.initialize.apply(this, arguments);
	};
	ClosedPath.prototype = {
		$a:null,
		initialize:function(array){
			if(!array){array = [];}
			this.$a = array;
			this.class_id=CLOSED_PATH;
		},
		/**
		 * Get quad bezier curve point.
		 */
		getBezierPoint:function(p0, p1, p2, p3, t){
			var pow = Math.pow;
			return new Point(
				pow(1 - t, 3)  * p0.x + 3 * t * pow(1 - t, 2) * p1.x + 3 * t * t * (1 - t) * p2.x + t * t * t * p3.x,
				pow(1 - t, 3)  * p0.y + 3 * t * pow(1 - t, 2) * p1.y + 3 * t * t * (1 - t) * p2.y + t * t * t * p3.y
			);
		},
		/**
		 * outpath
		 */
		outpath:function(){
			if(this._outpath != undefined){
				return this._outpath;
			}

			var pathtextbits = [];
			var $a = this.$a;
			var pt = $a[$a.length - 1].c[2];
			pathtextbits.push([M, pt.x, SPACE, pt.y, SPACE].join(NULL_STR));
			var prevTag;
			if(!param.isRelative){
				for(let i = 0, len=$a.length; i < len; i=(i+1)|0){
					let c = $a[i];
					let cc = c.c;
					if(c.tag == POTRACE_CORNER){
						pathtextbits.push([prevTag != c.tag ? L : NULL_STR, cc[1].x, SPACE, cc[1].y, SPACE, /*L,*/ cc[2].x, SPACE, cc[2].y, SPACE].join(NULL_STR));
					}else{
						pathtextbits.push([prevTag != c.tag ? C : NULL_STR, cc[0].x, SPACE, cc[0].y, SPACE, cc[1].x, SPACE, cc[1].y, SPACE, cc[2].x, SPACE, cc[2].y, SPACE].join(NULL_STR));
					}
					pt = cc[2];
					prevTag = c.tag;
				}
			}else{
				//precision
				let base = Math.pow(10, param.precision);
				let r = function(val){
					return (Math.round(val * base))/base;
				}
				for(let i = 0, len=$a.length; i < len; i=(i+1)|0){
					let cx = pt.x;
					let cy = pt.y;
					let c = $a[i];
					let cc = c.c;
					if(c.tag == POTRACE_CORNER){
						pathtextbits.push([prevTag != c.tag ? _l : NULL_STR, r(cc[1].x - cx), SPACE, r(cc[1].y - cy), SPACE, /*_l,*/ r(cc[2].x - cc[1].x), SPACE, r(cc[2].y - cc[1].y), SPACE].join(NULL_STR));
					}else{
						pathtextbits.push([prevTag != c.tag ? _c : NULL_STR, r(cc[0].x - cx), SPACE, r(cc[0].y - cy), SPACE, r(cc[1].x - cx), SPACE, r(cc[1].y - cy), SPACE, r(cc[2].x - cx), SPACE, r(cc[2].y - cy), SPACE].join(NULL_STR));
					}
					pt = cc[2];
					prevTag = c.tag;
				}
			}
			pathtextbits.push(Z);
			this._outpath = pathtextbits
				.join(NULL_STR)
				.replace(/ ([A-z,\-])/g, "$1")
				.replace(/(\s|-)0\./g, "$1.");
			return this._outpath;
		},
		_outpath: undefined
	};
	/////////////////////////////////////////
	// Point3D class
	////////////////////////////////////////
	var Point3D = function(){
	    this.initialize.apply(this, arguments);
	};
	Point3D.prototype = {
		initialize:function(x, y, z){
			if(!x){ x=0;}
			if(!y){ y=0;}
			if(!z){ z=0;}
			this.$v = new Float64Array([x, y, z]);
		},
		toString:function(){
			return ["(", [this.$v[0],this.$v[1],this.$v[2]].join(COMMA), ")"].join(NULL_STR);
		},
		$v:null
	};
	/////////////////////////////////////////
	// QuadraticForm class
	/**
	 *  the type of (affine) quadratic forms, represented as symmetric 3x3
	 *  matrices.  The value of the quadratic form at a vector (x,y) is v^t
	 *  Q v, where v = (x,y,1)^t.
	 */
	////////////////////////////////////////
	var QuadraticForm = function(){
		this.initialize.apply(this, arguments);
	};
	QuadraticForm.prototype = {
		initialize:function(){
			this.$m = [new Point3D(), new Point3D(), new Point3D()];
		},
		/**
		 *  Apply quadratic form Q to vector w = (w.x,w.y)
		 */
		apply:function(w){
			var v = new Float64Array([w.x, w.y, 1]);
			var sum = 0.0;
			var $m = this.$m;
			for(var i = 0; i < 3; i=(i+1)|0){
				var $mi$v = $m[i].$v;
				var vi = v[i];
				for(var j = 0; j < 3; j=(j+1)|0){
					sum += vi * $mi$v[j] * v[j];
				}
			}
			return sum;
		},
		clone:function(){
			var ret = new QuadraticForm();
			var ret_$m = ret.$m;
			var $m = this.$m;
			for(var i = 0; i < 3; i=(i+1)|0){
				var ret_$mi$v = ret_$m[i].$v;
				var $mi$v = $m[i].$v;
				for(var j = 0; j < 3; j=(j+1)|0){
					ret_$mi$v[j] = $mi$v[j];
				}
			}
			return ret;
		},
		add:function(m2){
			var $m = this.$m;
			for(var i= 0; i < 3; i=(i+1)|0){
				var $mi$v = $m[i].$v;
				var m2$mi$v = m2.$m[i].$v;
				for(var j= 0; j < 3; j=(j+1)|0){
					$mi$v[j] += m2$mi$v[j];
				}
			}
			return this;
		},
		scalar:function(s){
			var $m = this.$m;
			for(var i = 0; i < 3; i=(i+1)|0){
				var $mi$v = $m[i].$v;
				for(var j = 0; j < 3; j=(j+1)|0){
					$mi$v[j] *= s;
				}
			}
			return this;
		},
		fromVectorMultiply:function(v){
			var $m = this.$m;
			var v$v = v.$v;
			for(var i= 0; i < 3; i=(i+1)|0){
				var $mi$v = $m[i].$v;
				for(var j = 0; j < 3; j=(j+1)|0){
					$mi$v[j] = v$v[i] * v$v[j];
				}
			}
			return this;
		},
		toString:function(){
			return [
				"[" ,
				[
					this.$m[0].$v[0], this.$m[0].$v[1], this.$m[0].$v[2], 
					this.$m[1].$v[0], this.$m[1].$v[1], this.$m[1].$v[2],
					this.$m[2].$v[0], this.$m[2].$v[1], this.$m[2].$v[2] 
				].join(COMMA),
				"]"
			].join(NULL_STR);
		},
		$m:null
	};

	/////////////////////////////////////////
	// Sums class
	////////////////////////////////////////
	var Sums = function(){
	    this.initialize.apply(this, arguments);
	};
	Sums.prototype = {
		initialize:function(){
			this._values = new Float64Array(5);
			this.x = this.y = this.x2 = this.xy = this.y2 = 0;
		}
	};
	(function(proto){
		Object.defineProperty(proto, "x", {
			get: function(){return this._values[0];},
			set: function(value){this._values[0] = value;}
		});
		Object.defineProperty(proto, "y", {
			get: function(){return this._values[1];},
			set: function(value){this._values[1] = value;}
		});
		Object.defineProperty(proto, "x2", {
			get: function(){return this._values[2];},
			set: function(value){this._values[2] = value;}
		});
		Object.defineProperty(proto, "xy", {
			get: function(){return this._values[3];},
			set: function(value){this._values[3] = value;}
		});
		Object.defineProperty(proto, "y2", {
			get: function(){return this._values[4];},
			set: function(value){this._values[4] = value;}
		});
	})(Sums.prototype);
	/////////////////////////////////////////
	// ClosedPathList class
	////////////////////////////////////////
	var ClosedPathList = function(){
	    this.initialize.apply(this, arguments);
	};
	ClosedPathList.prototype = {
		initialize:function(){
			this.$a =[];
		},
		/**
		 *  trace bitmap
		 */
		trace:function(bmpdata){
			var pList = pathList(bmpdata);
			return processPath(pList);
		},
		/**
		 *  outpath.
		 */
		outpath:function(){
			return this._outpath;
		},
		/**
		 * toPathElements
		 */
		toPathElements:function(){
			var gElem = getGElement();
			var pathElemTmpl = getPathElement();
			var pathElem;
			for(var i=0, len=this.$a.length; i<len; i=(i+1)|0){
				var path = this.$a[i];
				if(path.hasOwnProperty(CLASS_ID) &&  path.class_id == CLOSED_PATH){
					pathElem = pathElemTmpl.cloneNode(false);
					pathElem.setAttribute(D, path.outpath());
					pathElem.setAttribute(CLASS, path.isEven? EVEN: ODD);
					gElem.appendChild(pathElem);
				}
			}
			return gElem;
		},
		/**
		 * toPathElement
		 */
		toPathElement:function(){
			var pathElem = getPathElement();
			pathElem.setAttribute(D, this.toPathString());
			return pathElem;
		},
		/**
		 * toPathString
		 */
		toPathString:function(){
			var outpath = this.outpath();
			return outpath != "" ? outpath: "M0 0";
		},
		$a:null,
		_outpath:undefined
	};
	/////////////////////////////////////////
	// Point class
	////////////////////////////////////////
	var Point = function () {
	    this.initialize.apply(this, arguments);
	};
	Point.prototype = {
		initialize:function(x,y){
			if(!this._values){
				this._values = new Float64Array(2);
			}
			if(!x){ x=0;}
			this.x = x;
			if(!y){ y=0;}
			this.y = y;
		},
		offset:function(x,y){
			this.x +=x;
			this.y +=y;
		},
		transformPoint:function(matrix){
			var x = this.x;
			var y = this.y;
			this.x = x*matrix[0] + y*matrix[2];
			this.y = x*matrix[1] + y*matrix[3];
		},
		clone:function(){
			return new Point(this.x, this.y);
		},
		equals:function(point){
			return (this.x == point.x && this.y == point.y);
		}
	};
	(function(proto){
		Object.defineProperty(proto, "x", {
			get: function(){return this._values[0];},
			set: function(value){this._values[0] = value;}
		});
		Object.defineProperty(proto, "y", {
			get: function(){return this._values[1];},
			set: function(value){this._values[1] = value;}
		});
	})(Point.prototype);

	////////////////////////////////////////
	// BitMap class
	/////////////////////////////////////////
	var BitMap = function () {
		this.initialize.apply(this, arguments);
	};
	BitMap.prototype = {
		initialize:function(canvas){
			if(canvas instanceof Uint8ClampedArray){
				this.__data = canvas;
				this.width = arguments[1];
				this.height = arguments[2];
				return;
			}else{
				this.height = canvas.height;
				this.width = canvas.width;
				this._canvas = canvas;
				this._context = canvas.getContext('2d');
				this._data = this._context.getImageData(0, 0, this.width, this.height);
				this.__data = this._data.data;
			}
		},
		getPixel:function(x, y){
			var data = this.__data;
			var index = x * 4 + y * this.width * 4;
			return new Pixel(data[index], data[index + 1], data[index + 2], data[index + 3]);
		},
		applyFilter:function(rect, pt, f){//NOTE:pt?
			var data = this.__data;
			for(var y = rect.y, yLim = rect.y + rect.height; y < yLim; y=(y+1)|0){
				for(var x = rect.x, xLim = rect.x + rect.width;x < xLim; x=(x+1)|0){
					var index = x * 4 + y * this.width * 4;
					var r = +data[index],
						g = +data[index + 1] ,
						b = +data[index + 2] ,
						a = +data[index + 3] ;
					data[index    ] = r*f[0] + g*f[1] + b*f[2] + a*f[3] + f[4]; //Rnew
					data[index + 1] = r*f[5] + g*f[6] + b*f[7] + a*f[8] + f[9]; //Gnew
					data[index + 2] = r*f[10] + g*f[11] + b*f[12] + a*f[13] + f[14]; //Bnew
					data[index + 3] = r*f[15] + g*f[16] + b*f[17] + a*f[18] + f[19]; //Anew
				}
			}
		},
		binarization:(function(){
			//gray scale
			var filter = new Float64Array([
				0.298912, 0.586611, 0.114478, 0, 0,
				0.298912 ,0.586611, 0.114478 ,0 ,0,
				0.298912 ,0.586611, 0.114478 ,0 ,0,
				0 , 0, 0, 1, 0
			]);
			return function(threshold){
				this.applyFilter(new Rectangle(0, 0, this.width, this.height), null, filter);
				threshold = 255*threshold;
				var data = this.__data;
				for(var i = 0, len =  data.length; i < len; i=(i+4)|0){
					data[i  ] = data[i+1] = data[i+2] = threshold > data[i] ? 0: 255;
				}
			}
		})(),
		clone:function(){
			if(!this._canvas){
				var array = new Uint8ClampedArray(this.__data.length);
				array.set(this.__data);
				return new BitMap(array, this.width, this.height);
			}else{
				var canvas = this._canvas.cloneNode(false);
				var context = canvas.getContext('2d');
				context.putImageData(this._data,0,0);
				var bmd = new BitMap(canvas);
				return bmd;
			}
		},
		reflesh:function(){
			this._context.putImageData(this._data, 0, 0);
		},
		height:0,
		width:0,
		_data:null,
		_context:null,
		_canvas:null
	};
	/////////////////////////////////////////
	// Pixel class
	////////////////////////////////////////
	var Pixel = function(){
		this.initialize.apply(this, arguments);
	}
	Pixel.prototype = {
		initialize: function(r,g,b,a){
			this._values = new Uint8Array(4);
			this.r = r;
			this.g = g;
			this.b = b;
			this.a = a;
		}
	};
	(function(proto){
		Object.defineProperty(proto, "r", {
			get: function(){return this._values[0];},
			set: function(value){this._values[0] = value;}
		});
		Object.defineProperty(proto, "g", {
			get: function(){return this._values[1];},
			set: function(value){this._values[1] = value;}
		});
		Object.defineProperty(proto, "b", {
			get: function(){return this._values[2];},
			set: function(value){this._values[2] = value;}
		});
		Object.defineProperty(proto, "a", {
			get: function(){return this._values[3];},
			set: function(value){this._values[3] = value;}
		});
	})(Pixel.prototype);
	/////////////////////////////////////////
	// Rectangle class
	////////////////////////////////////////
	var Rectangle = function () {
		this.initialize.apply(this, arguments);
	};
	Rectangle.prototype = {
		initialize:function(x,y,w,h){
			if(!this._values){this._values = new Float64Array(4);}
			if(!x){ x=0;}
			this.x = x;
			if(!y){ y=0;}
			this.y = y;
			if(!h){ h=0;}
			this.height =h;
			if(!w){ w=0;}
			this.width = w;
		}
	};
	(function(proto){
		Object.defineProperty(proto, "x", {
			get: function(){return this._values[0];},
			set: function(value){this._values[0] = value;}
		});
		Object.defineProperty(proto, "y", {
			get: function(){return this._values[1];},
			set: function(value){this._values[1] = value;}
		});
		Object.defineProperty(proto, "width", {
			get: function(){return this._values[2];},
			set: function(value){this._values[2] = value;}
		});
		Object.defineProperty(proto, "height", {
			get: function(){return this._values[3];},
			set: function(value){this._values[3] = value;}
		});
	})(Rectangle.prototype);
	//////////////////////////////////////////////////////////
	//public
	//////////////////////////////////////////////////////////
	var trace = function(img){
		//arguments must be HTMLImageElement or HTMLCanvasElement or Uint8ClampedArray(with width and height)
		if(img instanceof Uint8ClampedArray || img instanceof Uint8Array){
			var arr = wrapArray(img, arguments[1], arguments[2]);
			return traceArray(arr, arguments[1]+2, arguments[2]+2);
		}else{
			//NOTE:image must be wrapped by white pixels.
			var canvas = document.createElement(CANVAS);
			var ctx = canvas.getContext('2d');
			canvas.width = img.width + MARGIN * 2;
			canvas.height = img.height + MARGIN * 2;
			ctx.fillStyle = "white";
			ctx.fillRect(0, 0, canvas.width, canvas.height);
			ctx.drawImage(img, MARGIN, MARGIN, img.width, img.height);
			return traceCanvas(canvas);
		}
	};
	var traceCanvas = function(canvas){
		var bmd = new BitMap(canvas);
		bmd.binarization(param.threshold);
		var pList = pathList(bmd, param);
		var curveList = processPath(pList, param);
		return curveList;
	};
	var wrapArray = function(array, width, height){
		var narr = new Uint8ClampedArray((width+MARGIN*2)*(height+MARGIN*2)*4);
		narr.fill(255);
		for(let y = 0; y < height; y = (y+1)|0){
			for(let x = 0; x < width; x = (x+1)|0){
				let pos = (y * width + x)*4;
				let npos = ((y+1) * (width + MARGIN*2) + x + MARGIN)*4;
				let a = array[pos+3]/255;
				narr[npos] = array[pos] * a + 255*(1 - a);
				narr[npos+1] = array[pos+1] * a + 255*(1 - a);
				narr[npos+2] = array[pos+2] * a + 255*(1 - a);
			}
		}
		return narr;
	};
	var traceArray = function(array, width, height){
		var bmd = new BitMap(array, width, height);
		if(param.noFilter){
			bmd.binarization(param.threshold);
		}
		var pList = pathList(bmd, param);
		var curveList = processPath(pList, param);
		return curveList;
	};
	var setParam = function(_param){
		if(_param.threshold !== undefined){param.threshold = +_param.threshold;}
		if(_param.turdSize !== undefined){param.turdSize = +_param.turdSize;}
		if(_param.turnPolicy !== undefined){param.turnPolicy = +_param.turnPolicy;}
		if(_param.alphamax !== undefined){param.alphamax = +_param.alphamax;}
		if(_param.precision !== undefined){param.precision = ~~(_param.precision);}
		if(_param.noFilter !== undefined){param.noFilter = !!_param.noFilter;}
		if(_param.isRelative !== undefined){param.isRelative = !!_param.isRelative;}
	};
	var getParam = function(){
		return param;
	};
	return {
		trace:trace,
		setParam:setParam,
		getParam:getParam
	};
}());

//for Node
if(typeof exports != "undefined"){
	for(var i in Potrace){
		exports[i] = Potrace[i];
	}
	exports.normalizePixel = data => {
		const arr = new Uint32Array(data.buffer);
		return new Uint8Array((arr.map((val, i) => (val & 0xff000000) == 0 ? 0 : val)).buffer);
	};
}
