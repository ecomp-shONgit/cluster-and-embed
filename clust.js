/**

    2022 Prof. Schäfer Uni Trier, Prof. Schubert Uni Leipzig
    JS library for some cluster/dimension reduction functions

    Cluster analysis / SORTING: 
    hierar clust (neighb joining) - done
    Boostrap Consensus Tree (? needed)

    DIMENSION REDUCTION / low dim embedding:
    own MD - done
    PCA
    MDS - done
    tSNE
    Isomap 
    Locally Linear Embedding 
    Spectral Embedding/Laplacian Eigenmaps

**/

"use strict";

/*helper some day need to be seperated*/
function getIndexatMin( arr, i ){ 
	let MinX = Infinity;
    let ail = arr.length;
    let X = 0;
    let mX = -1;
    for( X = 0; X < ail; X+=1 ){
        if( MinX > arr[X] && arr[X] != NaN && i != X ){
            MinX = arr[X];
            mX = X;
		}
	}
    return mX;
}

function getIndexatMax( arr, i ){ 
	let MinX = 0;
    let ail = arr.length;
    let X = 0;
    let mX = -1;
    for( X = 0; X < ail; X+=1 ){
        if( MinX < arr[X] && arr[X] != NaN && i != X ){
            MinX = arr[X];
            mX = X;
		}
	}
    return mX;
}

function getIJatMax( arr, i ){ //arr quadratic
    let MaxX = 0;
    let ail = arr.length;
    let X = 0;
    let Y = 0;
    let mX = -1;
    let mY = -1;
    for( X = 0; X < ail; X += 1 ){
        if( i != X ){
            for( Y = X; Y < ail; Y += 1){
                if( i != Y ){
                    if( MaxX < arr[X][Y] && arr[X][Y] != NaN ){
                        MaxX = arr[X][Y];
                        mX = X;
                        mY = Y;
		            }
                }
            }
        }
	}
    return [mX, mY];

}


function indexnotin(A, i){
	let ail = A.length;
    let Ai = 0;
	for( Ai = 0; Ai < ail; Ai+=1 ){
		if( A[Ai] === i ){
			return false;
		}
	}
	return true;
}


/*
    hierarchical cluster analysis  (bottom up (joining)) 
    clustbez = 0 ( single linkage,shortest linkage of clusters, equal to neighbor jojning )
    clustbez = 1 ( UPGMA, avaerage distance of cluster)
*/
function sortDML(DM, L){
    //avoid wrong clusterings because of Reihenfolge, Eigenvalues???
    let eigval = [];
    
    let ll = DM.length;
    for( let i = 0; i < ll; i += 1 ){
        let sumitup = 0;
        let tt = DM[i].length;
        for( let j = 0; j < tt; j += 1 ){
            sumitup += DM[i][j];
        }
        sumitup /= tt;
        eigval.push( sumitup );
    }
    
    let sorteon = true;
    let sortedDM = [];//new Float64Array(ll); //return ;
    let sortedL = [];
    let neworder = [];
    while( sorteon ){
            let topic = getIndexatMin( eigval, Infinity );
            if( topic != -1 ){ //no more to pic
                neworder.push( topic );
                sortedL.push( L[topic] );
                eigval[topic] = NaN;
            } else {
                sorteon = false;
                break;
            }
    }
    for(let t = 0; t < ll; t += 1 ){
        let i = neworder[t];
        let newrow = [];
        for(let tt = 0; tt < ll; tt += 1 ){
            let j = neworder[tt];
            newrow.push(DM[i][j]);
        }
        sortedDM.push(newrow);
    }
    return [sortedDM, sortedL];
}

function clusthierarch( DM, L, clustbez ){
    
    let temp = sortDML(DM, L);
    DM = temp[0];
    L = temp[1];
    let ll = L.length;
    let clusters = [];
    //build first cluster layer - singele nodes
    let f1clust = [];
    for( let i = 0; i < ll; i += 1 ){
        f1clust.push([i]);
    }
    clusters.push(f1clust);
    //build all other cluster layers 
    let goonbuilding = true;
    while( goonbuilding ){
        let lastclusterlayer = clusters[ clusters.length-1 ];
        let inuse = [];
        let newclusterlayer = [];
        for(let l = 0; l < lastclusterlayer.length; l += 1){
            
            let memlc = 0;
            
            let mindist = Infinity;
            for(let ll = 0; ll < lastclusterlayer[l].length; ll += 1){
                let i = lastclusterlayer[l][ll];
                if(i == NaN){
                    continue
                }
                if( clustbez == 0 ){ //single linkage
                    for(let lc = 0; lc < lastclusterlayer.length; lc += 1){
                        if(lc != l){
                            for(let llc = 0; llc < lastclusterlayer[lc].length; llc += 1){
                                let ii = lastclusterlayer[lc][llc];
                                //console.log(i, ii, lc, llc)
                                if( mindist > DM[i][ii] ){
                                    memlc = lc;
                                    mindist = DM[i][ii];
                                }
                            }
                        }
                    }
                } else { //averade distance of clusters
                    for(let lc = 0; lc < lastclusterlayer.length; lc += 1){
                        if(lc != l){
                            let aver = 0;
                            for(let llc = 0; llc < lastclusterlayer[lc].length; llc += 1){
                                let ii = lastclusterlayer[lc][llc];
                                aver += DM[i][ii]
                                
                            }
                            aver /= lastclusterlayer[lc].length;
                            if( mindist >  aver ){
                                memlc = lc;
                                mindist = aver;
                            }
                        }
                    }
                }

            }
            
            let insuelc = indexnotin(inuse, memlc);
            let insuel = indexnotin(inuse, l);
            if( insuel &&  insuelc ){
                inuse.push(l);
                inuse.push(memlc);
                newclusterlayer.push(lastclusterlayer[l].concat(lastclusterlayer[memlc]));
            } else {
                if( insuel ){
                    inuse.push(l);
                    newclusterlayer.push(lastclusterlayer[l].concat([]));
                } else if( insuelc ){
                    inuse.push(memlc);
                    newclusterlayer.push(lastclusterlayer[memlc].concat([]));
                }
            }
        }
        clusters.push( newclusterlayer );
        if( clusters[ clusters.length-1 ].length == 1 ){
            goonbuilding = false;
        }
    }
    //console.log(clusters);
    return [clusters, L];
}

/*
    dimensional reduction
    based on angles and points
*/
function usemeasure( v1, v2 ){
    return euclideanM( v1, v2 );
}

function genpointoncircle( cx, cy, r, angel ){
    let px = Math.round( cx + ( r*Math.cos( angel * ( Math.PI/180 ) ) ) );
    let py = Math.round( cy + ( r*Math.sin( angel * ( Math.PI/180 ) ) ) );
    return [px,py];
}

function getDM( V ){ //computer distance matrix - should move to vecspacemeasure lib
    let DM = [];
    let ll = V.length;
    for(let i = 0; i < ll; i += 1 ){
        let row = [];
        for( let j = 0; j < ll; j += 1 ){
            row.push( usemeasure(V[i], V[j] ) );
        }
        DM.push( row );
    }
    return DM;
}

function getDisssimofDist( DM1, DM2 ){ //DM1.length == DM2.length !!!!
    let DD = [];
    let ll = DM1.length;
    for( let i = 0; i < ll; i += 1 ){
        let dd = [];
        for( let j = 0; j < ll; j += 1 ){
             dd.push( Math.abs( DM1[i][j] - DM2[i][j]) );
        }
        DD.push( dd );
    }
    return DD;
}   

function avversumM( M ){//avarage sum of matrix
    let aver = 0;
    let ll = M.length;
    for( let i = 0; i < ll; i += 1 ){
        for( let j = 0; j < ll; j += 1 ){
            aver += M[i][j];
        }
    }
    return (aver/(ll*ll));
}
function avversumV( V ){//average sum of vector
    let aver = 0;
    let ll = V.length;
    for( let i = 0; i < ll; i += 1 ){       
            aver += V[i];
    }
    return (aver/ll);
}

function MDsomething2D( DM, acc ){ //prperties unkonw, passes small test
    let ll = DM.length;
    //init result points
    let RP = [];
    let fp = [100, 100];
    RP.push( fp );
    let angs = [0,]
    for( let i = 1; i < ll; i += 1 ){
        let point = genpointoncircle( fp[0], fp[1], DM[0][i], 94 );
        //console.log(fp[0], fp[1], DM[0][i], L[i]);
        RP.push( point );
        fp[0] = point[0];
        fp[1] = point[1];
    }
    //console.log(RP.toString());
    let goon = true;
    let seccount = 0;
    let allerrorlast = NaN;
    let optistep = 4;
    while( goon ){
        let DM2 = getDM( RP );
        //console.log(DM2.join(" "));
        let DD = getDisssimofDist( DM, DM2 );
        //console.log(DD.join(" "));
        let inin = getIJatMax( DD, Infinity ); //get i j - of matrix pairs
        let mintempdist = Infinity;
        let optipoint = undefined; 
        for( let ang = 0; ang < 360; ang += optistep ){
            let temppoint = genpointoncircle( RP[inin[1]][0], RP[inin[1]][1], DM[inin[0]][inin[1]], ang );
            let tempd = 0;
            for(let i = 0; i < ll; i += 1){ //minimale distanz zu allen punkten
                //if(i != inin[0] && i != inin[1]){
                    tempd += Math.abs(usemeasure(temppoint, RP[i]) - DM[i][inin[0]]);
                //}
            }
            tempd /= (ll-2);
            if( tempd < mintempdist ){
                optipoint = temppoint;
                mintempdist = tempd;
            }
        }
        let allerror = avversumM( DD );
        //console.log(seccount, inin[0], inin[1], allerror, allerrorlast);
        if( allerrorlast < allerror && allerror < acc){
            //console.log(DM2);
            goon = false;
        } else {
            
            //console.log(mintempdist, optipoint, RP[inin[0]])
            RP[inin[0]][0] = ((RP[inin[0]][0] + optipoint[0] )/2);
            RP[inin[0]][1] = ((RP[inin[0]][1] + optipoint[1] )/2);
            //RP[inin[0]][0] = (RP[inin[0]][0] + (step*optipoint[0]) );
            //RP[inin[0]][1] = (RP[inin[0]][1] + (step*optipoint[1]) );
        }
        allerrorlast = allerror;
        if( seccount == 1000 ){
            console.log(DM2);
            goon = false;
        }
        seccount += 1;
    }
    //console.log(DM);
    //console.log(RP);
    return RP; //check what need to be returned to build stylo like output
}

/*
    dimensional reduction
    MDS classical
*/
function getEMatrix(dim){ //Einheitsmatrix
    let mr = [];
    for( let i = 0; i < dim; i += 1 ){
        let row = [];
        for( let j = 0; j < dim; j += 1 ){
            row.push(1.0);
        }
        mr.push( row );
    }
    return mr;
}

function getNMatrix(dim){ //Zeromatrix
    let mr = [];
    for( let i = 0; i < dim; i += 1 ){
        let row = [];
        for( let j = 0; j < dim; j += 1 ){
            row.push(0.0);
        }
        mr.push( row );
    }
    return mr;
}

function getIMatrix( dim ){ //diagonal matrix
    let mr = [];
    for( let i = 0; i < dim; i += 1 ){
        let row = [];
        for( let j = 0; j < dim; j += 1 ){
            if( i == j ){
                row.push(1.0);
            } else {
                row.push(0.0);
            }
        }
        mr.push( row );
    }
    return mr;
}

function scalerProdV( v1, v2 ){
    //given len(v1) == len(v2)
    let l = v1.length;
    let prodsu = 0.0;
    for( let i = 0; i < l; i += 1 ){
        prodsu += (v1[i]*v2[i]);
    }
    return prodsu; 
}

function scalrmult(s, m){
    let l = m.length;
    let rm = getNMatrix( l ); 
    for( let i = 0; i < l; i += 1 ){
        for( let j = 0; j < l; j += 1 ){
            rm[i][j] = m[i][j]*s;
        }
    }
    return rm;
}

function innerProddistmatrix( M ){
    let l = M.length;
    let IPM = getNMatrix(l); //init
    for( let i = 0; i < l; i += 1 ){
        for( let j = 0; j < l; j += 1 ){
            IPM[i][j] = scalerProdV(M[i], M[j]);
        }
    }
    return IPM;
}
function innerProdM( m1, m2 ){
    let l = m1.length;
    let m11 = getNMatrix(l); //init
    for( let i = 0; i < l; i += 1 ){
        for( let j = 0; j < l; j += 1 ){
            m11[i][j] = m1[j][i];
        }
    }
    let rm = getNMatrix(l); //init
    for( let i = 0; i < l; i += 1 ){
        for( let j = 0; j < l; j += 1 ){
            rm[i][j] = scalerProdV(m11[i], m2[j]);
        }
    }
    return rm;
}

function matrixsubt( m1, m2 ){
    let l = m1.length;
    let rm = getNMatrix(l); //init
    for( let i = 0; i < l; i += 1 ){
        for( let j = 0; j < l; j += 1 ){
            rm[i][j] = (m1[i][j]-m2[i][j]);
        }
    }
    return rm;
}


function MDS( DM ){
    //MD square shaped
    let ll = DM.length;
    //compute gramsche matrix from DM
    let B = innerProddistmatrix( DM );
    console.log(DM, B);
    //build centering matrix
    let C = matrixsubt( getIMatrix( ll ) , scalrmult( (1/ll), getEMatrix( ll ) ) ); //varianz
    console.log(C);
    //use double centering with centering matrix computed from dimension and the I an E Matrix
    let cB = scalrmult( (-1/2), innerProdM( innerProdM( C, B ), C) ); //covarianz
    
    let eigenvalvec = numeric.eig( cB, 100 );
    console.log(cB, eigenvalvec);
}

/*
    dimensional reduction
    tSNE classical
*/

/*
    tSNE main fkt
    INPUT: 
    howmuch:itteration to call step, perple: perpelxity, di: resulting dimensions, epi: epsilon, X: array of high dim points, D: matix of diffences 
    examp: tSNE(100, 30, 2, 10, [],[] );

*/

//globals to store things
let perplexity = 30;
let dim = 2;
let epsilon = 10;
let iter = 0;

//length of current input
let LL = NaN;
let P = []; //typed needed

//outputs
let Y = null;
let gains = null;
let ystep = null;

// return 0 mean unit standard deviation random number
let return_v = false;
let v_val = 0.0;

function gaussRandom( ){
    if( return_v ) { 
        return_v = false;
        return v_val; 
    }
    let u = 2*Math.random()-1;
    let v = 2*Math.random()-1;
    let r = u*u + v*v;
    if(r == 0 || r > 1) return gaussRandom();
    let c = Math.sqrt(-2*Math.log(r)/r);
    v_val = v*c; // cache this for next function call for efficiency
    return_v = true;
    return u*c;
}

// return random normal number
function randn( mu, std ){ 
    return mu+gaussRandom()*std; 
}

// utilitity that creates contiguous vector of zeros of size n
function zeros( n ){
    if(typeof(n)==='undefined' || isNaN(n)) { return []; }
    if(typeof ArrayBuffer === 'undefined') {
        // lacking browser support
        let arr = new Array(n);
        for(let i=0;i<n;i++) { arr[i]= 0; }
        return arr;
    } else {
        return new Float64Array(n); // typed arrays are faster
    }
}

// utility that returns 2d array filled with random numbers
// or with value s, if provided
function randn2d( n, d, s ){
    let uses = typeof s !== 'undefined';
    let x = [];
    for( let i = 0; i < n; i+= 1 ){
        let xhere = [];
        for( let j = 0; j < d; j+= 1 ){ 
            if( uses ){
                xhere.push(s); 
            } else {
                xhere.push(randn(0.0, 1e-4)); 
            }
        }
        x.push(xhere);
    }
    return x;
}

// compute pairwise distance in all vectors in X
function xtod( X ){
    let N = X.length;
    let dist = zeros(N * N); // allocate contiguous array
    for( let i = 0; i < N; i+=1 ){
        for( let j = i+1; j < N; j+= 1 ){
            let d = L2(X[i], X[j]);
            dist[i*N+j] = d;
            dist[j*N+i] = d;
        }
    }
    return dist;
}

// compute (p_{i|j} + p_{j|i})/(2n)
function d2p(D, perplexity, tol) {
    let Nf = Math.sqrt( D.length ); // this better be an integer
    let N = Math.floor( Nf );
    if( N !== Nf ){ 
        throw "D should have square number of elements.";
    } 
    let Htarget = Math.log( perplexity ); // target entropy of distribution
    let P = zeros( N * N ); // temporary probability matrix

    let prow = zeros( N ); // a temporary storage compartment
    for( let i = 0; i < N; i+=1 ){
        let betamin = -Infinity;
        let betamax = Infinity;
        let beta = 1; // initial value of precision
        let done = false;
        let maxtries = 50;

        // perform binary search to find a suitable precision beta
        // so that the entropy of the distribution is appropriate
        let num = 0;
        while( !done ){
            //debugger;

            // compute entropy and kernel row with beta precision
            let psum = 0.0;
            for( let j=0; j<N; j+=1 ){
                let pj = Math.exp(- D[i*N+j] * beta);
                if( i === j ){ pj = 0; } // we dont care about diagonals
                prow[j] = pj;
                psum += pj;
            }

            // normalize p and compute entropy
            let Hhere = 0.0;
            for( let j=0; j<N; j+=1 ){
                let pj = 0;
                if( psum === 0 ){
                    pj = 0;
                } else {
                    pj = prow[j] / psum;
                }
                prow[j] = pj;
                if( pj > 1e-7 ) Hhere -= pj * Math.log(pj);
            }

            // adjust beta based on result
            if( Hhere > Htarget ){
                // entropy was too high (distribution too diffuse)
                // so we need to increase the precision for more peaky distribution
                betamin = beta; // move up the bounds
                if( betamax === Infinity ){ 
                    beta = beta * 2; 
                } else { 
                    beta = (beta + betamax) / 2; 
                }

            } else {
                // converse case. make distrubtion less peaky
                betamax = beta;
                if( betamin === -Infinity ){ 
                    beta = beta / 2; 
                } else { 
                    beta = (beta + betamin) / 2; 
                }
            }

            // stopping conditions: too many tries or got a good precision
            num++;
            if( Math.abs( Hhere - Htarget ) < tol ){ 
                done = true; 
            }
            if( num >= maxtries ){ 
                done = true; 
            }
        } //while end

        // console.log('data point ' + i + ' gets precision ' + beta + ' after ' + num + ' binary search steps.');
        // copy over the final prow to P at row i
        for( let j=0; j<N; j += 1 ){ P[i*N+j] = prow[j]; }

    } // end loop over examples i

    // symmetrize P and normalize it to sum to 1 over all ij
    let Pout = zeros(N * N);
    let N2 = N*2;
    for( let i=0; i<N; i+=1 ){
        for( let j=0; j<N; j+=1 ){
            Pout[i*N+j] = Math.max((P[i*N+j] + P[j*N+i])/N2, 1e-100);
        }
    }

    return Pout;
}

// helper function
function sign(x) { return x > 0 ? 1 : x < 0 ? -1 : 0; }



function initSolution( ){
    // generate random solution to t-SNE
    Y = randn2d( LL, dim); // the solution
    gains = randn2d( LL, dim, 1.0); // step gains to accelerate progress in unchanging directions
    ystep = randn2d( LL, dim, 0.0); // momentum accumulator
    iter = 0;
}

function step( ){
    iter += 1;
    let cg = costGrad( ); // evaluate gradient
    let cost = cg.cost;
    let grad = cg.grad;

    // perform gradient step
    let ymean = zeros( dim );
    for( let i = 0; i < LL; i += 1 ){
        for( let d = 0; d < dim; d += 1 ){
            let gid = grad[ i ][ d ];
            let sid = ystep[i][d];
            let gainid = gains[i][d];

            // compute gain update
            let newgain = sign(gid) === sign(sid) ? gainid * 0.8 : gainid + 0.2;
            if(newgain < 0.01) newgain = 0.01; // clamp
            gains[i][d] = newgain; // store for next turn

            // compute momentum step direction
            let momval = iter < 250 ? 0.5 : 0.8;
            let newsid = momval * sid - epsilon * newgain * grad[i][d];
            ystep[i][d] = newsid; // remember the step we took

            // step!
            Y[i][d] += newsid; 

            ymean[d] += Y[i][d]; // accumulate mean so that we can center later
        }
    }

    // reproject Y to be zero mean
    for( let i = 0; i < LL; i += 1 ){
        for( let d = 0; d < dim; d += 1 ){
            Y[i][d] -= ymean[d]/LL;
        }
    }

    //if(this.iter%100===0) console.log('iter ' + this.iter + ', cost: ' + cost);
    return cost; // return current cost
}

function debugGrad( ){
    let cg = costGrad( ); // evaluate gradient
    let cost = cg.cost;
    let grad = cg.grad;

    let e = 1e-5;
    for( let i = 0; i < LL; i += 1 ){
        for( let d = 0; d < dim; d += 1 ){
        let yold = Y[i][d];

        Y[i][d] = yold + e;
        let cg0 = costGrad( Y );

        Y[i][d] = yold - e;
        let cg1 = costGrad( Y );

        let analytic = grad[i][d];
        let numerical = (cg0.cost - cg1.cost) / ( 2 * e );
        console.log(i + ',' + d + ': gradcheck analytic: ' + analytic + ' vs. numerical: ' + numerical);

        Y[i][d] = yold;
        }
    }
}

function costGrad( ){
    let pmul = iter < 100 ? 4 : 1; // trick that helps with local optima

    // compute current Q distribution, unnormalized first
    let Qu = zeros(LL * LL);
    let qsum = 0.0;
    for( let i = 0; i < LL; i += 1 ){
        for( let j = i+1; j < LL; j += 1 ){
            let dsum = 0.0;
            for( let d = 0; d < dim; d+=1 ){
                let dhere = Y[i][d] - Y[j][d];
                dsum += dhere * dhere;
            }
            let qu = 1.0 / (1.0 + dsum); // Student t-distribution
            Qu[i*LL+j] = qu;
            Qu[j*LL+i] = qu;
            qsum += 2 * qu;
        }
    }
    // normalize Q distribution to sum to 1
    let NN = LL*LL;
    let Q = zeros(NN);
    for( let q=0; q < NN; q+=1 ){ Q[q] = Math.max(Qu[q] / qsum, 1e-100); }

    let cost = 0.0;
    let grad = [];
    for( let i=0; i < LL; i+=1 ){
        let gsum = new Array(dim); // init grad for point i
        for( let d=0; d<dim; d+=1 ){ gsum[d] = 0.0; }
        for( let j=0; j<LL; j++ ){
            cost += - P[i*LL+j] * Math.log(Q[i*LL+j]); // accumulate cost (the non-constant portion at least...)
            let premult = 4 * (pmul * P[i*LL+j] - Q[i*LL+j]) * Qu[i*LL+j];
            for( let d=0; d < dim; d+=1 ){
                gsum[d] += premult * (Y[i][d] - Y[j][d]);
            }
        }
        grad.push(gsum);
    }

    return {cost: cost, grad: grad};
}

function tSNE( howmuch, perple, di, epi, D ) { 
    perplexity = perple; // effective number of nearest neighbors
    dim = di; // by default 2-D tSNE
    epsilon = epi; // learning rate
    let lD = D.length;
    if(!lD > 0){
        throw " D is empty? You must have some data!";
    }
    LL = lD; // back up the size of the dataset
    // convert D to a (fast) typed array version
    let dists = zeros( lD * lD ); // allocate contiguous array
        for(let i = 0; i < lD; i+=1 ){
        for(let j = i+1; j < lD; j+=1 ){
            let d = D[i][j];
            dists[i*LL+j] = d;
            dists[j*LL+i] = d;
        }
    }
    P = d2p( dists, perplexity, 1e-4); // attach to object
    
    initSolution( ); // refresh this
    
    //steps
    for( let goon = 0; goon < howmuch; goon += 1 ){
        step();
    }
    // return pointer to current solution
    return Y;
}





//mini test case
function runclusttest(){
    let texts = ["a","b", "c", "d", "e"];
    let testdistma = 
    [
    [0, 5, 9, 9, 8],
    [5, 0, 10, 10, 9],
    [9, 10, 0, 8, 7],
    [9, 10, 8, 0, 3],
    [8, 9, 7, 3, 0]
    ];


    let clu = clusthierarch( testdistma, texts, 1 );
    console.log(clu)
    for(let cl = 0; cl < clu[0].length; cl+=1){
        console.log("Clusterlayer", cl);
        for(let ci = 0; ci < clu[0][cl].length; ci+=1){
            console.log("Cluster", ci);
            for(let i = 0; i < clu[0][cl][ci].length; i += 1){
                console.log(clu[1][clu[0][cl][ci][i]]);
            }
        }
    }

    let tw = 2000;
    let th = 2000;
    let addx = Math.round(tw/2);
    let addy = Math.round(th/2);
    let mainsvgelem = getsvgMAINELEM( tw.toString(), th.toString() );
    //let dimred1 = MDsomething2D( testdistma, 2.3 ); // distance matrix, label vector, accuracy - how scaled????
    //let dimred1 = MDS( testdistma ); // distance matrix, label vector, target dimension
    let dimred1 = tSNE(300, 2, 2, 10, testdistma ); 
    for( let i = 0; i < dimred1.length; i += 1) {
        let nn = texts[i];
        let nl = ""+nn;
        let posx = Math.round( (dimred1[i][0]*4)+addx );
        let posy = Math.round( (dimred1[i][1]*4)+addy );
        if(posx > tw || posy > th ){
            console.log(nn, "at", posx, posy );
        }
        console.log(posx, posy, nn);

        //svg
        mainsvgelem.appendChild( getdot( posx, posy, 10, nn, nl,"red", "blue" ) );
    }
    document.getElementById("resultsofclust").appendChild( mainsvgelem );
}

