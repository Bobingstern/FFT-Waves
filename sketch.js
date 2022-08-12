CameraControls.install( { THREE: THREE } );

document.addEventListener("keydown", onDocumentKeyDown, false);
let skyTog = true
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
let A = 0.00025
let N = 32
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
const geometry = new THREE.PlaneBufferGeometry(L, L, 100, 100 );
const material = new THREE.MeshBasicMaterial( {color: 0x00a6ff, side: THREE.DoubleSide, wireframe: true} );

var dummyRGBA = [];
let dummyRGBA2 = []
for(var i=0; i < (N+1)**2; i++){
    dummyRGBA[i] = [ocean.h_tilde_0[i].x * 150000, ocean.h_tilde_0[i].y * 150000, 0, 255];
    dummyRGBA2[i] = [ocean.h_tilde_0[i].z * 100000, ocean.h_tilde_0[i].w * 100000, 0, 255];
    //dummyRGBA[i] = [255, 0, 255, 255]
}

//flatten dummyRBGA to a single array of floats
//dummyRGBA[N+2] = [0, 0, 0, 255]
dummyRGBA.reverse()
dummyRGBA = [].concat.apply([], dummyRGBA);


let ida = new ImageData(Uint8ClampedArray.from(dummyRGBA), N+1, N+1);
let dummyTex = new THREE.Texture(ida, THREE.UVMapping, THREE.ClampToEdgeWrapping, THREE.ClampToEdgeWrapping, THREE.NearestFilter);
dummyTex.needsUpdate = true


const cvs = document.createElement("canvas");
cvs.width = cvs.height = (N+1)**2;
const ctx = cvs.getContext("2d");

ctx.putImageData(ida, 0, 0);
cvs.id = "monke"

document.body.appendChild(cvs);

let ctT = new THREE.CanvasTexture(cvs)

const oceanShader = new THREE.ShaderMaterial({
    wireframe: false, 
    side : THREE.DoubleSide,
    transparent: true,
    uniforms: {
        t : {value: 0},
        A : {value: A},
        L : {value: L},
        wind: {value: wind},
        h_tilde_0: {value: ocean.h_tilde_0},
        offset : {value: new THREE.Vector2(0,0)},
        tex: {value: dummyTex}
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
    uniform sampler2D tex;


    // vec2 h_tilde[N];
    // vec2 h_tilde_slopex[N];
    // vec2 h_tilde_slopez[N];
    // vec2 h_tilde_dx[N];
    // vec2 h_tilde_dz[N];

    uniform vec2 offset;
    varying vec3 vNormal;


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
        //vec4 htilde0 = texture2D(tex, vec2(float(n)/float(N+1), float(m)/float(N+1))).rgba;
        float omegat = dispersion(n, m) * t;

        float cos_ = cos(omegat);
        float sin_ = sin(omegat);

        vec2 c0 = vec2(cos_, sin_);
        vec2 c1 = vec2(cos_, -sin_);
        return cadd(cmult(c0, htilde0.xy), cmult(c1, htilde0.zw));
    }
    int bitRev(int x, int log2n){
        int n = 0;
        for (int i=0;i<log2n;i++){
            n <<= 1;
            n |= (x & 1);
            x >>= 1;
        }
        return n;
    }
    void fft(inout vec2 a[N], inout vec2 A[N]){
        int log2n = int(log(float(N))/log(2.0));
        for (int i=0;i<N;i++){
            int rev = bitRev(i, log2n);
            A[i] = a[rev];
        }

        vec2 J = vec2(0, 1);
        for (int s=1; s<=log2n;++s){
            int m = 1 << s;
            int m2 = m >> 1;
            vec2 w = vec2(0, 1);
            vec2 wm = vec2(cos(PI / float(m2)), sin(PI / float(m2)));
            for (int j=0;j<m2;j++){
                for (int k=j;k<N;k+=m){
                    vec2 t = vec2(w.x*A[k+m2].x - w.y*A[k+m2].y, w.x*A[k+m2].y - w.y*A[k+m2].x);
                    vec2 u = A[k];
                    A[k] = vec2(u.x+t.x, u.y+t.y);
                    A[k+m2] = vec2(u.x-t.x, u.y-t.y);
                }
                w = vec2(w.x*wm.x - w.y*wm.y, w.x*wm.y - w.y*wm.x);
            }
        }
    }      
    void main() {
        vec2 h;
        vec2 D;
        vec2 c;
        vec2 hTilde_c;
        vec2 k;
        vec3 no;
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

                no += vec3(-kx * hTilde_c.y, 0, -kz * hTilde_c.y);

                if (k_len > 0.000001){
                    D = cadd(D, vec2(kx * hTilde_c.y / k_len, kz * hTilde_c.y / k_len));
                }
            }
        }
        no *= -1.0;
        no += vec3(0,1,0);
        no = normalize(no);
        vNormal = no;
        gl_Position = projectionMatrix * modelViewMatrix * vec4(position.x + D.x, position.y + D.y, h.x * -1.0, 1.0);
    }
    `,
    fragmentShader : `
    //color blue
    varying mediump vec3 vNormal;
    uniform sampler2D tex;
    void main() {
        mediump vec3 light = vec3(1.0, 0.4, 0.0);

        // ensure it's normalized
        light = normalize(light);

        // calculate the dot product of
        // the light to the vertex normal
        mediump float dProd = max(0.0, dot(vNormal, light));

        // feed into our frag colour
        vec3 c = vec3(0, 102, 200);
        gl_FragColor = vec4(dProd * (c.x/255.0), // R
                            dProd * (c.y/255.0), // G
                            dProd * (c.z/255.0), // B
                            0.9);  // A


        //gl_FragColor = vec4(texture2D(tex, vec2(0.0/3.0, 0.0/3.0)));
    }
    `
})

const plane = new THREE.Mesh( geometry, oceanShader );


plane.quaternion.copy( camera.quaternion );

plane.rotation.x = -Math.PI / 2;

let planes = []
let s = 4
for (let i=0;i<s;i++){
    for (let j=0;j<s;j++){
        let p = new THREE.Mesh( geometry.clone(), oceanShader.clone() );
        
        p.quaternion.copy( camera.quaternion );
        p.rotation.x = -Math.PI / 2;
        p.position.x = L * (i) - L/2 * (s-1);
        p.position.z = L * (j) - L/2 * (s-1);
        p.material.uniforms.offset.value = new THREE.Vector2(L/2 * (s-1), L/2 * (s-1))
        planes.push(p)
        scene.add(p)
    }
}
// scene.add(plane)
// scene.add(plane2)
let sky = new THREE.Sky()
sky.scale.setScalar( 1000 );

let sun = new THREE.Vector3()
const phi = THREE.MathUtils.degToRad( 90-2 );
const theta = THREE.MathUtils.degToRad( 180 );

sun.setFromSphericalCoords( 2, phi, theta );
sky.material.uniforms['sunPosition'].value.copy(sun)
sky.name = "sky"
scene.add(sky)

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