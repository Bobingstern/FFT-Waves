let t = 0
function setup(){
    createCanvas(1920,1080);
}

function troch(x,y,amp,freq, speed){
    let X = amp*(height-y)*cos(freq*(speed*t-x*0.01))+x
    let Y = amp*(height-y)*sin(freq*(speed*t-x*0.01))+y
    return createVector(X,Y)
}
function draw(){
    background(220);
    fill(0,0,255,100)
    stroke(0,0,255)
    beginShape()
    let Y_START = height/2
    for (let y=Y_START;y<height;y+=50){
        for (let x=-100;x<width+100;x+=50){
            
            let amp = 0.2
            let freq = 0.4
            pos = troch(x,y,amp,freq, 0.05)
            //pos.y -= Y_START
            
            push()
            noStroke()
            fill(0,0,255, y == Y_START ? 255 : 50)
            circle(pos.x, pos.y, 10)
            pop()
            if (y == Y_START){
                vertex(pos.x,pos.y)
            }
        }
    }
    vertex(width,height)
    vertex(0,height)
    endShape()
    t+=1
}