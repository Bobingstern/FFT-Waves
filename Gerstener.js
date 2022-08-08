
class Wave{
    constructor(dir, steepness, wavelength){
        this.t = 0
        this.dir = dir
        this.steepness = steepness
        this.wavelength = wavelength

    }
    getWave(x, y, z){
        let p = new THREE.Vector3(x,y,z);
        let k = 2*Math.PI/this.wavelength;
        let c = Math.sqrt(9.8/k)
        let d = new THREE.Vector2(this.dir.x, this.dir.y)
        d.normalize()
        let f = k * ((d.x*p.x + d.y*p.y) - c * this.t);
        let a = this.steepness / k
        p.z = a * Math.sin(f);
        p.x = d.x * a * Math.cos(f)
        p.y = d.y * a * Math.cos(f)
        return p
    }
    update(delta){
        this.t += delta
    }
}

class MultipleWaves{
    constructor(){
        this.waves = []
    }
    addWave(dir, steepness, wavelength){
        this.waves.push(new Wave(dir, steepness, wavelength))
    }
    getWave(x, y, z){
        let p = new THREE.Vector3(x,y,z)
        for (let w of this.waves){
            let p2 = w.getWave(x,y,z)
            p.x += p2.x
            p.y += p2.y
            p.z += p2.z
        }
        return p
    }
    update(delta){
        for (let w of this.waves){
            w.update(delta)
        }
    }
}