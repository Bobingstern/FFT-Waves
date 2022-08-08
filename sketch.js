CameraControls.install( { THREE: THREE } );


const concat = (xs, ys) => xs.concat(ys);


let T = 0
let A = 0.0005
let N = 32
let L = 64
let wind = new THREE.Vector2(32, 0)

const width = window.innerWidth;
const height = window.innerHeight;
const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera( 60, width / height, 0.1, 5000 );
camera.position.set( 0, 0, 64 );

const renderer = new THREE.WebGLRenderer();
renderer.setSize( width, height );
//render.domElement.id = "waves"
document.body.appendChild( renderer.domElement );

let ocean = new Ocean(N, A, wind, L)

let clock = new THREE.Clock();
const controls = new CameraControls(camera, renderer.domElement);
const geometry = new THREE.PlaneBufferGeometry(L, L, N, N );
const material = new THREE.MeshBasicMaterial( {color: 0x00a6ff, side: THREE.DoubleSide, wireframe: true} );

var dummyRGBA = new Uint8Array(4 * ((N+1)**2));
for(var i=0; i < 4 *(N+1)**2; i+=4){
    dummyRGBA[i] = ocean.h_tilde_0[i/4].x * 100000;
    dummyRGBA[i+1] = ocean.h_tilde_0[i/4].y * 100000;
    dummyRGBA[i+2] = ocean.h_tilde_0[i/4].z * 100000;
    dummyRGBA[i+3] = ocean.h_tilde_0[i/4].w * 100000;
}
//flatten dummyRBGA to a single array of floats


let ida = new ImageData(Uint8ClampedArray.from(dummyRGBA), N+1, N+1);
let dummyTex = new THREE.Texture(ida);


const oceanShader = new THREE.ShaderMaterial({
    wireframe: true, 
    side : THREE.DoubleSide,
    uniforms: {
        t : {value: 0},
        A : {value: A},
        L : {value: L},
        wind: {value: wind},
        h_tilde_0: {value: ocean.h_tilde_0},
        offset : {value: new THREE.Vector2(0,0)},
        test: {value: dummyTex}
    },
    vertexShader : `
    const float PI = 3.141592653589793;
    const float g = 9.81;
    const int N = 32;

    uniform float t;
    uniform float A;
    uniform int L;
    uniform vec2 wind;
    uniform vec4 h_tilde_0[(N+1)*(N+1)];


    // vec2 h_tilde[N];
    // vec2 h_tilde_slopex[N];
    // vec2 h_tilde_slopez[N];
    // vec2 h_tilde_dx[N];
    // vec2 h_tilde_dz[N];

    uniform vec2 offset;
    


    //complex multiply
    vec2 cmult(vec2 a, vec2 b) {
        return vec2(a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x);
    }
    //complex add
    vec2 cadd(vec2 a, vec2 b) {
        return vec2(a.x+b.x, a.y+b.y);
    }

    float dispersion(int n, int m) {
        float w_0 = 2.0 * PI / 200.0;
        float kx = PI * (2.0 * float(n) - float(N)) / float(L);
        float kz = PI * (2.0 * float(m) - float(N)) / float(L);
        return floor(sqrt(g * sqrt(kx*kx + kz*kz)) / w_0) * w_0;
    }
    vec2 hTilde(int n, int m) {
        int index = m * (N+1) + n;
        vec4 htilde0 = h_tilde_0[index];
        float omegat = dispersion(n, m) * t;

        float cos_ = cos(omegat);
        float sin_ = sin(omegat);

        vec2 c0 = vec2(cos_, sin_);
        vec2 c1 = vec2(cos_, -sin_);
        return cadd(cmult(c0, htilde0.xy), cmult(c1, htilde0.zw));
    }
   
    void main() {
        vec2 h;
        vec2 D;
        vec2 c;
        vec2 hTilde_c;
        vec2 k;
        for (int m = 0; m < N; m++) {
            float kz = PI * (2.0 * float(m) - float(N)) / float(L);
            for (int n = 0; n < N; n++) {
                float kx = PI * (2.0 * float(n) - float(N)) / float(L);
                k = vec2(kx, kz);
                float k_len = length(k);
                float k_dot_x = dot(k, vec2(position.x+offset.x, position.y));
                c = vec2(cos(k_dot_x), sin(k_dot_x));
                hTilde_c = cmult(hTilde(n, m), c);
                h = cadd(h, hTilde_c);

                if (k_len > 0.000001){
                    D = cadd(D, vec2(kx * hTilde_c.y / k_len, kz * hTilde_c.y / k_len));
                }
            }
        }
        
        gl_Position = projectionMatrix * modelViewMatrix * vec4(position.x + D.x, position.y + D.y, h.x * -1.0, 1.0);
    }
    `,
    fragmentShader : `
    //color blue
    void main() {
        gl_FragColor = vec4(0.0, 1.0, 1.0, 1.0);
    }
    `
})

const plane = new THREE.Mesh( geometry, oceanShader );
const plane2 = new THREE.Mesh( geometry, oceanShader.clone() );

let GEOM = plane.geometry.attributes.position;

plane.quaternion.copy( camera.quaternion );
plane2.quaternion.copy( camera.quaternion );

plane.rotation.x = -Math.PI / 2;
plane2.rotation.x = -Math.PI / 2;

plane2.position.x += L
plane2.material.uniforms.offset.value = new THREE.Vector2(L,0)

let planes = []
let s = 4
for (let i=0;i<s;i++){
    for (let j=0;j<s;j++){
        let p = new THREE.Mesh( geometry.clone(), oceanShader.clone() );
        
        p.quaternion.copy( camera.quaternion );
        p.rotation.x = -Math.PI / 2;
        p.position.x = L * (i) - L/2 * (s-1);
        p.position.z = L * (j) - L/2 * (s-1);
        p.material.uniforms.offset.value = new THREE.Vector2(L*i,L*j)
        planes.push(p)
        scene.add(p)
    }
}
// scene.add(plane)
// scene.add(plane2)



function anim () {
	requestAnimationFrame( anim );
    let delta = clock.getDelta();
    T = clock.getElapsedTime();
    let updated = controls.update( delta );
    
    // for (let i=0;i<GEOM.count;i++){

    //     let p = ocean.verticies[i]
    //     GEOM.setX(i, p.x)
    //     GEOM.setY(i, p.z)
    //     GEOM.setZ(i, p.y)
    // }
    
    // plane.geometry.attributes.position.needsUpdate = true;
    // ocean.evalWavesFFT(T)
    for (let p of planes){
        p.material.uniforms.t.value = T;
    }
    // plane.material.uniforms.t.value = T
    // plane2.material.uniforms.t.value = T
	renderer.render( scene, camera );

	
};

anim()