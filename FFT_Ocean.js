function gaussRnd() {
    var u = 0, v = 0;
    while(u === 0) u = Math.random(); //Converting [0,1) to (0,1)
    while(v === 0) v = Math.random();
    return Math.sqrt( -2.0 * Math.log( u ) ) * Math.cos( 2.0 * Math.PI * v );
}


class FFT{
    constructor(N){
        this.N = N
        this.which = 0
        this.log_2_N = Math.log(this.N)/Math.log(2)
        this.pi2 = 2*Math.PI
        this.reversed = []
        this.T = []
        this.c = []

        for (let i=0;i<this.N;i++) this.reversed[i] = this.reverse(i)

        let pow2 = 1
        for (let i=0;i<this.log_2_N;i++){
            let temp = []
            for (let j=0;j<pow2;j++){
                temp.push(this.t(j, pow2 * 2))
            }
            pow2 *= 2
            this.T.push(temp)
        }

        let temp = []
        for (let i=0;i<this.N;i++){
            temp.push(new math.complex(0,0))
        }
        this.c.push(temp); this.c.push(temp)
    }
    reverse(i){
        let res = 0
        for (let j=0;j<this.log_2_N;j++){
            res = (res << 1) + (i & 1)
            i >>= 1
        }
        return res
    }
    t(x, N){
        return math.complex(Math.cos(this.pi2 * x / N), Math.sin(this.pi2 * x / N))
    }
    fft(input, output, stride, offset){
        for (let i=0;i<this.N;i++) {
            this.c[this.which][i] = input[this.reversed[i] * stride + offset]
        }

        let loops = this.N >> 1
        let size = 1 << 1
        let size_over_2 = 1
        let w_ = 0
        
        for (let i=1;i<=this.log_2_N;i++){
            this.which ^= 1
            for (let j=0;j<loops;j++){
                for (let k=0;k<size_over_2;k++){
                    this.c[this.which][size * j + k] = math.add(this.c[this.which^1][size * j + k], math.multiply(this.c[this.which ^ 1][size * j + k + size_over_2], this.T[w_][k]))
                }
                for (let k=size_over_2;k<size;k++){
                    this.c[this.which][size * j + k] = math.subtract(this.c[this.which^1][size * j - size_over_2 + k], math.multiply(this.c[this.which ^ 1][size * j + k], this.T[w_][k - size_over_2]))
                }
            }
            loops >>= 1
            size <<= 1
            size_over_2 <<= 1
            w_++
        }
        for (let i=0;i<this.N;i++) {
            output[i * stride + offset] = this.c[this.which][i]
        }
    }
}

class vertex_ocean {
    constructor() {
        this.x = 0 
        this.y = 0
        this.z = 0
        this.nx = 0
        this.ny = 0
        this.nz = 0
        this.a = 0;
        this.b = 0;
        this.c = 0;
        this._a = 0;
        this._b = 0;
        this._c = 0;
        this.ox = 0
        this.oy = 0
        this.oz = 0
    }
}

class complex_vector_normal{
    constructor() {
        this.h = math.complex(0, 0)
        this.D = new THREE.Vector2()
        this.n = new THREE.Vector3()

    }
}

class Ocean{
    constructor(N, A, w, length){
        this.g = 9.8
        this.N = N
        this.A = A
        this.w = w
        //this.windspeed = wsp
        this.length = length

        this.h_tilde = []
        this.h_tilde_slopex = []
        this.h_tilde_dx = []
        this.h_tilde_slopez = []
        this.h_tilde_dz = []
        this.h_tilde_0 = []
        for (let i=0;i<this.N;i++){
            this.h_tilde.push(math.complex(0,0))
            this.h_tilde_slopex.push(math.complex(0,0))
            this.h_tilde_dx.push(math.complex(0,0))
            this.h_tilde_slopez.push(math.complex(0,0))
            this.h_tilde_dz.push(math.complex(0,0))
        }

        this.verticies = []
        this.indicies = []
        this.indicies_count = 0
        this.vbo_verticies;
        this.vbo_indicies;

        this.fft = new FFT(this.N)


        for (let m_prime=0;m_prime<this.N+1;m_prime++){
            for (let n_prime=0;n_prime<this.N+1;n_prime++){
                let index = m_prime * (this.N+1) + n_prime
                let h_tilde_0 = this.hTilde_0(n_prime, m_prime)
                let h_tilde_conj = math.conj(this.hTilde_0(-n_prime, -m_prime))
                this.verticies.push(new vertex_ocean())
                this.verticies[index].a = h_tilde_0.re
                this.verticies[index].b = h_tilde_0.im
                this.verticies[index]._a = h_tilde_conj.re
                this.verticies[index]._b = h_tilde_conj.im

                this.h_tilde_0.push(new THREE.Vector4(h_tilde_0.re, h_tilde_0.im, h_tilde_conj.re, h_tilde_conj.im))

                this.verticies[index].ox  = (n_prime - this.N / 2) * this.length / this.N
                this.verticies[index].oy = 0
                this.verticies[index].oz  = (m_prime - this.N / 2) * this.length / this.N
                this.verticies[index].x = this.verticies[index].ox
                this.verticies[index].y = this.verticies[index].oy
                this.verticies[index].z = this.verticies[index].oz
                this.verticies[index].nx = 0
                this.verticies[index].ny = 1
                this.verticies[index].nz = 0


            }
        }
        // this.indicies_count = 0
        // for (let m_prime=0;m_prime<this.N+1;m_prime++){
        //     for (let n_prime=0;n_prime<this.N+1;n_prime++){
        //         index = m_prime * (this.N+1) + n_prime
        //         this.indicies[this.indicies_count++] = index
        //         this.indicies[this.indicies_count++] = index + this.N+1
        //         this.indicies[this.indicies_count++] = index + this.N+2
        //         this.indicies[this.indicies_count++] = index
        //         this.indicies[this.indicies_count++] = index + this.N+2
        //         this.indicies[this.indicies_count++] = index + this.N+1
        //     }
        // }
    }
    dispersion(n, m){
        let w_0 = 2*Math.PI / 200
        let kx = Math.PI * (2 * n - this.N) / this.length
        let kz = Math.PI * (2 * m - this.N) / this.length
        return Math.floor(Math.sqrt(this.g * Math.sqrt(kx**2 + kz**2)) / w_0) * w_0
    }
    phillips(n, m){
        let k = new THREE.Vector2(Math.PI * (2 * n - this.N) / this.length, Math.PI * (2 * m - this.N) / this.length)
        let k_length = k.length()
        let k_length_2 = k_length ** 2
        let k_length_4 = k_length_2 ** 2
        if (k_length < 0.000001) return 0

        let k_norm = k.clone().normalize()
        let w_norm = this.w.clone().normalize()
        let k_dot_w = k_norm.dot(w_norm)
        let k_dot_w2 = k_dot_w ** 2
        let w_length = this.w.length()
        let L = w_length**2 / this.g
        let L2 = L ** 2

        let damping = 0.001
        let l2 = L2 * damping**2
        return this.A * Math.exp(-1 / (k_length_2 * L2)) / k_length_4 * k_dot_w2 * Math.exp(-k_length_2 * l2)
    }
    hTilde_0(n, m){
        let rnd = math.complex(gaussRnd(), gaussRnd())
        return math.multiply(rnd, Math.sqrt(this.phillips(n, m)) / 2)
    }

    hTilde(t, n_prime, m_prime){
        let index = m_prime * (this.N+1) + n_prime
        let htilde0 = math.complex(this.verticies[index].a, this.verticies[index].b)
        let htilde0mkconj = math.complex(this.verticies[index]._a, this.verticies[index]._b)

        let omegat = this.dispersion(n_prime, m_prime) * t
        
        let cos_ = Math.cos(omegat)
        let sin_ = Math.sin(omegat)

        let c0 = math.complex(cos_, sin_)
        let c1 = math.complex(cos_, -sin_) 
        let res = math.add(math.multiply(c0, htilde0), math.multiply(c1, htilde0mkconj))
        return res;
    }

    h_d_n(x, t){
        let h = math.complex(0,0)
        let D = new THREE.Vector2(0,0)
        let n = new THREE.Vector3(0,0,0)

        let c
        let htilde_c
        let k
        let kx
        let kz
        let k_length
        let k_dot_x

        for (let m_prime=0;m_prime<this.N;m_prime++){
            kz = 2 * Math.PI * (m_prime - this.N / 2) / this.length
            for (let n_prime=0;n_prime<this.N;n_prime++){
                kx = 2 * Math.PI * (n_prime - this.N / 2) / this.length
                k = new THREE.Vector2(kx, kz)
                k_length = k.length()
                k_dot_x = k.dot(x)
                c = math.complex(Math.cos(k_dot_x), Math.sin(k_dot_x))
                htilde_c = math.multiply(this.hTilde(t, n_prime, m_prime), c)
                h = math.add(h, htilde_c)

                n.add(new THREE.Vector3(-kx * htilde_c.im, 0, -kz * htilde_c.im))

                if (k_length < 0.000001) continue
                D.add(new THREE.Vector2(kx * htilde_c.im / k_length, kz * htilde_c.im / k_length))
            }
        }

        let s = new THREE.Vector3(0,1,0)
        n.multiplyScalar(-1)
        n.add(s)
        n.normalize()
        let cvn = new complex_vector_normal()
        cvn.h = h
        cvn.n = n
        cvn.D = D
        return cvn
        
    }

    evalWaves(t){
        let lambda = -1
        for (let m_prime=0;m_prime<this.N;m_prime++){
            for (let n_prime=0;n_prime<this.N;n_prime++){
                let index = m_prime * (this.N+1) + n_prime
                let x = new THREE.Vector2(this.verticies[index].x, this.verticies[index].z)

                let h_d_n = this.h_d_n(x, t)
                let scale = 1
                this.verticies[index].y = h_d_n.h.re/scale
                this.verticies[index].x = this.verticies[index].ox + lambda * h_d_n.D.x/scale
                this.verticies[index].z = this.verticies[index].oz + lambda * h_d_n.D.y/scale

                this.verticies[index].nx = h_d_n.n.x
                this.verticies[index].ny = h_d_n.n.y
                this.verticies[index].nz = h_d_n.n.z

                if (n_prime == 0 && m_prime == 0){
                    this.verticies[index + this.N + this.N*(this.N + 1)].y = h_d_n.a
                    this.verticies[index + this.N + this.N*(this.N + 1)].x = this.verticies[index + this.N + this.N*(this.N + 1)].ox + lambda * h_d_n.D.x
                    this.verticies[index + this.N + this.N*(this.N + 1)].z = this.verticies[index + this.N + this.N*(this.N + 1)].ox + lambda * h_d_n.D.y

                    this.verticies[index + this.N + this.N*(this.N + 1)].nx = h_d_n.n.x
                    this.verticies[index + this.N + this.N*(this.N + 1)].ny = h_d_n.n.y
                    this.verticies[index + this.N + this.N*(this.N + 1)].nz = h_d_n.n.z
                }
                if (n_prime == 0){
                    this.verticies[index + this.N].y = h_d_n.a
                    this.verticies[index + this.N].x = this.verticies[index + this.N].ox + lambda * h_d_n.D.x
                    this.verticies[index + this.N].z = this.verticies[index + this.N].ox + lambda * h_d_n.D.y

                    this.verticies[index + this.N].nx = h_d_n.n.x
                    this.verticies[index + this.N].ny = h_d_n.n.y
                    this.verticies[index + this.N].nz = h_d_n.n.z
                }
                if (n_prime == 0){
                    // this.verticies[index + this.N*(this.N + 1)].y = h_d_n.a
                    // this.verticies[index + this.N*(this.N + 1)].x = this.verticies[index + this.N*(this.N + 1)].ox + lambda * h_d_n.D.x
                    // this.verticies[index + this.N*(this.N + 1)].z = this.verticies[index + this.N*(this.N + 1)].ox + lambda * h_d_n.D.y

                    // this.verticies[index + this.N*(this.N + 1)].nx = h_d_n.n.x
                    // this.verticies[index + this.N*(this.N + 1)].ny = h_d_n.n.y
                    // this.verticies[index + this.N*(this.N + 1)].nz = h_d_n.n.z
                }
            }
        }
    }

    evalWavesFFT(t){
        let kx = -1
        let kz = -1
        let len = -1
        let lambda = -1
        let index = 0
        let index1 = 0

        for (let m_prime=0;m_prime<this.N;m_prime++){
            kz = Math.PI * (2 * m_prime - this.N) / this.length
            for (let n_prime = 0; n_prime < this.N ; n_prime++){
                kx = Math.PI * (2 * n_prime - this.N) / this.length
                len = Math.sqrt(kx * kx + kz * kz)
                index = m_prime * this.N + n_prime
                this.h_tilde[index] = this.hTilde(t, n_prime, m_prime)
                this.h_tilde_slopex[index] = math.multiply(this.h_tilde[index], math.complex(0, kx))
                this.h_tilde_slopez[index] = math.multiply(this.h_tilde[index], math.complex(0, kz))
                if (len < 0.000001){
                    this.h_tilde_dx[index] = math.complex(0,0)
                    this.h_tilde_dz[index] = math.complex(0,0)
                }
                else{
                    this.h_tilde_dx[index] = math.multiply(this.h_tilde[index], math.complex(0, -kx/len))
                    this.h_tilde_dz[index] = math.multiply(this.h_tilde[index], math.complex(0, -kz/len))
                }
            }
        }
        this.h_tilde = fft2(this.h_tilde)
        this.h_tilde_slopex = fft2(this.h_tilde_slopex)
        this.h_tilde_slopez = fft2(this.h_tilde_slopez)
        this.h_tilde_dx = fft2(this.h_tilde_dx)
        this.h_tilde_dz = fft2(this.h_tilde_dz)
        
        
        
        let sign = 0
        let signs = [1, -1]
        let n = new THREE.Vector3()
        for (let m_prime=0;m_prime<this.N;m_prime++){
            for (let n_prime = 0; n_prime < this.N ; n_prime++){
                index = m_prime * this.N + n_prime
                index1 = m_prime * (this.N+1) + n_prime
                sign = signs[(n_prime + m_prime) & 1]
                this.h_tilde[index] = math.multiply(this.h_tilde[index], sign)
                this.verticies[index1].y = this.h_tilde[index].re

                this.h_tilde_dx[index] = math.multiply(this.h_tilde_dx[index], sign)
                this.h_tilde_dz[index] = math.multiply(this.h_tilde_dz[index], sign)
                this.verticies[index1].x = this.verticies[index1].ox + this.h_tilde_dx[index].re * lambda
                this.verticies[index1].z = this.verticies[index1].oz + this.h_tilde_dz[index].re * lambda

                this.h_tilde_slopex[index] = math.multiply(this.h_tilde_slopex[index], sign)
                this.h_tilde_slopez[index] = math.multiply(this.h_tilde_slopez[index], sign)
                n = new THREE.Vector3(0 - this.h_tilde_slopex[index].re, 1, 0 - this.h_tilde_slopez[index].re).normalize()
                this.verticies[index1].nx = n.x
                this.verticies[index1].ny = n.y
                this.verticies[index1].nz = n.z
                if (n_prime == 0 && m_prime == 0) {
                    this.verticies[index1 + this.N + (this.N+1) * this.N].y = this.h_tilde[index].re;
     
                    this.verticies[index1 + this.N + (this.N+1) * this.N].x = this.verticies[index1 + this.N + (this.N+1) * this.N].ox + this.h_tilde_dx[index].re * lambda;
                    this.verticies[index1 + this.N + (this.N+1) * this.N].z = this.verticies[index1 + this.N + (this.N+1) * this.N].oz + this.h_tilde_dz[index].re * lambda;
                 
                    this.verticies[index1 + this.N + (this.N+1) * this.N].nx =  n.x;
                    this.verticies[index1 + this.N + (this.N+1) * this.N].ny =  n.y;
                    this.verticies[index1 + this.N + (this.N+1) * this.N].nz =  n.z;
                }
                if (n_prime == 0) {
                    this.verticies[index1 + this.N].y = this.h_tilde[index].a;
     
                    this.verticies[index1 + this.N].x = this.verticies[index1 + N].ox + this.h_tilde_dx[index].re * lambda;
                    this.verticies[index1 + this.N].z = this.verticies[index1 + N].oz + this.h_tilde_dz[index].re * lambda;
                 
                    this.verticies[index1 + this.N].nx =  n.x;
                    this.verticies[index1 + this.N].ny =  n.y;
                    this.verticies[index1 + this.N].nz =  n.z;
                }
                if (m_prime == 0) {
                    this.verticies[index1 + (this.N+1) * this.N].y = this.h_tilde[index].a;
     
                    this.verticies[index1 + (this.N+1) * this.N].x = this.verticies[index1 + (this.N+1) * this.N].ox + this.h_tilde_dx[index].re * lambda;
                    this.verticies[index1 + (this.N+1) * this.N].z = this.verticies[index1 + (this.N+1) * this.N].oz + this.h_tilde_dz[index].re * lambda;
                 
                    this.verticies[index1 + (this.N+1) * this.N].nx =  n.x;
                    this.verticies[index1 + (this.N+1) * this.N].ny =  n.y;
                    this.verticies[index1 + (this.N+1) * this.N].nz =  n.z;
                }

            }
        }

    }
}