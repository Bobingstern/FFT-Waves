CameraControls.install( { THREE: THREE } );

document.addEventListener("keydown", onDocumentKeyDown, false);
let skyTog = false
function onDocumentKeyDown(event) {
    var keyCode = event.which;
    if (keyCode == 83) {
        if (skyTog){
            scene.remove(sky)
            for (let p of planes){
                p.material.wireframe = true
            }   
        }
        else{
            scene.add(sky)
            for (let p of planes){
                p.material.wireframe = false
            }  
        }
        skyTog = !skyTog
    }
}
const concat = (xs, ys) => xs.concat(ys);


let T = 0
let A = 0.0002

let L = 64
let wind = new THREE.Vector2(32, 0)

const width = window.innerWidth;
const height = window.innerHeight;
const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera( 60, width / height, 0.1, 5000 );
camera.position.set( 0, 5, 64 );

const renderer = new THREE.WebGLRenderer();
renderer.setSize( width, height );
//render.domElement.id = "waves"
document.body.appendChild( renderer.domElement );

let ocean = new Ocean(N, A, wind, L)

let clock = new THREE.Clock();
const controls = new CameraControls(camera, renderer.domElement);
const geometry = new THREE.PlaneBufferGeometry(L, L, N, N );
const material = new THREE.MeshBasicMaterial( {color: 0x00a6ff, side: THREE.DoubleSide, wireframe: true} );


const oceanShader = new THREE.ShaderMaterial({
    wireframe: true, 
    side : THREE.DoubleSide,
    transparent: true,
    uniforms: {
        t : {value: 0},
    },
    vertexShader : `
    const float PI = 3.141592653589793;
    const float g = 9.81;
    const int N = 32;
    uniform float t;
    void main(){
    	gl_Position = projectionMatrix * modelViewMatrix * vec4(position.x, position.y, position.z, 1.0);
    }
    `,
    fragmentShader : `
    //color blue

    varying vec2 vUv;
    void main() {
        mediump vec3 light = vec3(1.0, 0.4, 0.0);

        // ensure it's normalized
        light = normalize(light);

        // calculate the dot product of
        // the light to the vertex normal
        //mediump float dProd = max(0.0, dot(vNormal, light));

        // feed into our frag colour
        // vec3 c = vec3(0, 132, 200);
        gl_FragColor = vec4(48.0 / 255.0, 182.0 / 255.0, 209.0 / 255.0, 0.9);  // A


        //gl_FragColor = texture2D(tex, vUv);
    }
    `
})

let planes = []
let s = 1
let plane = new THREE.Mesh( geometry.clone(), oceanShader.clone() );
plane.quaternion.copy( camera.quaternion );
plane.rotation.x = -Math.PI / 2;
plane.position.x =0;
plane.position.z = 0;
scene.add(plane)
// scene.add(plane2)

let sky = new THREE.Sky()
sky.scale.setScalar( 1000 );

let sun = new THREE.Vector3()
const phi = THREE.MathUtils.degToRad( 90-2 );
const theta = THREE.MathUtils.degToRad( 180 );

// sun.setFromSphericalCoords( 2, phi, theta );
// sky.material.uniforms['sunPosition'].value.copy(sun)
// sky.name = "sky"
// scene.add(sky)


//-------Ocean stuff

//0.02
let precomputed_pos = []

function precompute(){
    let time_ = 0
    let dt = 0.02
    let logg = 0
    while (time_ < 2.4){
        ocean.evalWavesFFT(time_)
        let temp = []
        for (let p of ocean.verticies){
            temp.push({"x":p.x, "y":p.y, "z":p.z})
        }
        precomputed_pos.push(temp)
        time_ += dt
        console.log(time_)
        if (logg > 62){
            logg = 0
            console.log(JSON.stringify(precomputed_pos))
        }
        logg++
    }
}
//precompute()
let indexx = 0
function anim () {
	setTimeout( function() {

        requestAnimationFrame( anim );

    }, 1000 / 60 );
    let delta = clock.getDelta();
    T = clock.getElapsedTime();
    let updated = controls.update( delta );
    let GEOM = plane.geometry.attributes.position
    ocean.evalWavesFFT(T)
    if (indexx >= precomputed_pos.length){
        indexx = 0
    }
    for (let i=0;i<GEOM.count;i++){
    	let p = ocean.verticies[i]
        GEOM.setXYZ(i, p.x, -p.z, p.y)
    }
    GEOM.needsUpdate = true
	renderer.render( scene, camera );
    indexx++
	
};

anim()