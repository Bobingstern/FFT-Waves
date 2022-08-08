CameraControls.install( { THREE: THREE } );

let N = 256
let L = 100
let A = 20
let PSIZE = 150
let w
let windspeed = 26
let g = 9.81
let T = 0
let h0k = []
let h0k_minus = []
let h0kt = []
let IM
let t
const concat = (xs, ys) => xs.concat(ys);

function setup(){
    //createCanvas(256,256)
    w = createVector(1,0)
    let H = make0()
    h0k = H[0]
    h0k_minus = H[1]
    for (let i=0;i<plane.geometry.attributes.position.count;i++){
        let x = plane.geometry.attributes.position.getX(i);
        let y = plane.geometry.attributes.position.getY(i);
        let z = plane.geometry.attributes.position.getZ(i);
        let k = createVector(2*PI*x/L, 2*PI*z/L)
        let h = h_kt(k, 0)
        h0kt.push(new Fourier.Complex(h[0], h[1]))
        //h0kt.push([constrain(round(h[0]*300),0,255), constrain(round(h[1]*300),0,255), 0, 255])
    }
    //h0kt = h0kt.reduce(concat)
    //IM = new ImageData(Uint8ClampedArray.from(h0kt), 16,16)
    anim()
}

function gaussRnd() {
    var u = 0, v = 0;
    while(u === 0) u = Math.random(); //Converting [0,1) to (0,1)
    while(v === 0) v = Math.random();
    return Math.sqrt( -2.0 * Math.log( u ) ) * Math.cos( 2.0 * Math.PI * v );
}

function make0(){
    let h0k = []
    let h0kminus = []
    let RES = 1
    for (let y=0;y<N;y++){
      let temp1 = []
      let temp2 = []
      for (let x=0;x<N;x++){
        let k = createVector(2*PI*x/L, 2*PI*y/L)
        let k_minus = createVector(-2*PI*x/L, -2*PI*y/L)
        let L_ = windspeed**2/g
        let ksq = k.magSq()
        
        let wind_dir = p5.Vector.normalize(w)
        let wk = p5.Vector.dot(p5.Vector.normalize(k), wind_dir)
        let wk_minus = p5.Vector.dot(p5.Vector.normalize(k_minus), wind_dir)
        let l = windspeed**2 / g
        let phs = exp((-1/(ksq*l*l)))/(ksq*ksq) * wk * wk
        let phs_minus = exp((-1/(ksq*l*l)))/(ksq*ksq) * wk_minus * wk_minus
        
        let re = gaussRnd()/sqrt(2) * sqrt(phs)
        let im = gaussRnd()/sqrt(2) * sqrt(phs)
        
        let re2 = gaussRnd()/sqrt(2) * sqrt(phs_minus)
        let im2 = gaussRnd()/sqrt(2) * sqrt(phs_minus)
        
        temp1.push([re,im])
        temp2.push([re2, im2])
        
      }
      h0k.push(temp1)
      h0kminus.push(temp2)
    }
    
    return [h0k, h0kminus]
}
  
function h_kt(k, t){
    let w = sqrt(g * k.mag())
    let p = p5.Vector.div(k, 2*PI/L)
    // (cos(w*t) + isin(w*t)) * (a+ib)
    // cos(wt)a -bsin(wt) + ibcos(wt) + iasin(wt) 

    let wt = w*t
    let h0 = h0k[round(p.y)+PSIZE/2][round(p.x)+PSIZE/2]
    let h0_minus = h0k_minus[round(p.y)+PSIZE/2][round(p.x)+PSIZE/2]
    let a = h0[0]
    let b = h0[1]
    let c = h0_minus[0]
    let d = h0_minus[1]

    let IM = (a+c)*sin(wt) + (b-d)*cos(wt)
    let RE = (a+c)*cos(wt) + (d-b)*sin(wt)
    return [RE, IM]

}

function h_kt_bogo(k, t){
    let w = sqrt(g * k.mag())
    let p = p5.Vector.div(k, 2*PI/L)
    // (cos(w*t) + isin(w*t)) * (a+ib)
    // cos(wt)a -bsin(wt) + ibcos(wt) + iasin(wt) 

    let wt = w*t
    let h0 = h0k[round(p.y)][round(p.x)]
    let h0_minus = h0k_minus[round(p.y)][round(p.x)]
    let a = h0[0]
    let b = h0[1]
    let c = h0_minus[0]
    let d = h0_minus[1]

    let IM = (a+c)*sin(wt) + (b-d)*cos(wt)
    let RE = (a+c)*cos(wt) + (d-b)*sin(wt)
    return [RE, IM]

}

const width = 1920;
const height = 1080;
const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera( 60, width / height, 0.1, 1000 );
camera.position.set( 0, 0, 100 );

const renderer = new THREE.WebGLRenderer();
renderer.setSize( width, height );
//render.domElement.id = "waves"
document.body.appendChild( renderer.domElement );


let clock = new THREE.Clock();
const controls = new CameraControls(camera, renderer.domElement);
const geometry = new THREE.PlaneBufferGeometry(PSIZE, PSIZE, 75, 75 );
const geometry2 = new THREE.PlaneBufferGeometry(PSIZE, PSIZE, 75, 75 );
const material = new THREE.MeshBasicMaterial( {color: 0x00a6ff, side: THREE.DoubleSide, wireframe: true} );
const plane = new THREE.Mesh( geometry, material );
const plane2 = new THREE.Mesh( geometry2, material );

plane.quaternion.copy( camera.quaternion );
plane2.quaternion.copy( camera.quaternion );

plane.rotation.x = -Math.PI / 2;
plane2.rotation.x = -Math.PI / 2;
//scene.add( plane );



//Wave stuff
let wave = new MultipleWaves();
wave.addWave(new THREE.Vector3(1,0), 0.25, 60);
wave.addWave(new THREE.Vector3(1,0.5), 0.25, 31);
wave.addWave(new THREE.Vector3(1,-0.5), 0.1, 45);
wave.addWave(new THREE.Vector3(-1,0.5), 0.1, 45);


let points = [];
for (let i=0;i<plane.geometry.attributes.position.count;i++){
    let cube = new THREE.Mesh(
        new THREE.BoxGeometry( 0.25, 0.25, 0.25 ), 
        new THREE.MeshBasicMaterial( { color: 0x00a6ff } )
      );
    let x = plane.geometry.attributes.position.getX(i);
    let y = plane.geometry.attributes.position.getY(i);
    let z = plane.geometry.attributes.position.getZ(i);
    cube.position.x = x
    cube.position.y = z
    cube.position.z = y
    points.push(cube)
    scene.add(cube)

    // let k = createVector(2*PI*x/L, 2*PI*y/L)
    // let h = h_kt(k, 0)
    // h0kt.push([h[0], h[1], 0, 255])
}

//flatten h0kt using reduce
//let h0kt_flat = h0kt.reduce((xs, ys) => xs.concat(ys)).map(h0kt).reduce((xs, ys) => xs.concat(ys))

// function draw() {
//     background(220)
//     let RES = 5
//     for (let y=0;y<N;y+=RES){
//       for (let x=0;x<N;x+=RES){
//         let k = createVector(2*PI*x/L, 2*PI*y/L)
//         let F = h_kt_bogo(k, t)
//         fill(F[0]*5000, F[1]*5000, 0)
//         rect(x,y,RES,RES)
//       }
//     }
    
// }

function anim () {
	requestAnimationFrame( anim );
    let delta = clock.getDelta();
    t = clock.getElapsedTime();
    
    let updated = controls.update( delta );
    let GEOM = plane.geometry.attributes.position;
    let GEOM2 = plane2.geometry.attributes.position;
    let EM = []
    // h0kt = []
    // for (let i=0;i<plane.geometry.attributes.position.count;i++){
    //     let x = plane.geometry.attributes.position.getX(i);
    //     let y = plane.geometry.attributes.position.getY(i);
    //     let z = plane.geometry.attributes.position.getZ(i);
    //     let k = createVector(2*PI*x/L, 2*PI*z/L)
    //     let h = h_kt(k, t)
    //     //h0kt.push(new Fourier.Complex(h[0]*5000, h[1]*5000))
    //     h0kt.push(new Fourier.Complex(constrain(round(h[0]*5000),0,255), constrain(round(h[1]*5000),0,255)))
    // }
    // Fourier.invert(h0kt, EM)
    for (let i=0;i<GEOM.count;i++){
        
        let x = GEOM2.getX(i);
        let y = GEOM2.getY(i);
        let z = GEOM2.getZ(i);
        let p = wave.getWave(x,y,z);
        
        points[i].position.x = p.x
        points[i].position.y = p.z
        points[i].position.z = p.y
        //points[i].position.y = EM[i]
        //GEOM.setXYZ(i,p.x,p.y,p.z);
        

    }
    //GEOM.needsUpdate = true;
    wave.update(delta)
	renderer.render( scene, camera );

	
};

