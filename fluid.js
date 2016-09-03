var WIDTH = 5;
var HEIGHT = 5;

var Fluid = function(width, height, density) {
  this.width = width;
  this.height = height;
  this.density = density;
  this.prop = new Grid2D(width,   height,     0); // property in each cell

  this.p0 = new Grid2D(width,      height,     0);
  this.p1 = new Grid2D(width,      height,     0);

  this.q0 = new Grid2D(width,      height,     0);
  this.q1 = new Grid2D(width,      height,     0);

  this.u  = new Grid2D(width + 1,  height,     0,  0.5, 0.0);
  this.v  = new Grid2D(width,      height + 1, 0,  0.0, 0.5);

  this.gridSize = 1 / Math.min(this.width, this.height);
}

Fluid.prototype = {
  swap: function() {
    var tmpq = this.q0;
    this.q0 = this.q1;
    this.q1 = tmpq;

    var tmpp = this.p0;
    this.p0 = this.p1;
    this.p1 = tmpp;
  },

  addInFlow: function(sx, sy, w, h, d, u, v) {
    var ssx = sx / this.gridSize;
    var ssy = sy / this.gridSize;

    var sex = (sx + w) / this.gridSize;
    var sey = (sy + h) / this.gridSize;

    this.q0.add(ssx, ssy, sex, sey, d);
    this.u.add(ssx, ssy, sex, sey, u);
    this.v.add(ssx, ssy, sex, sey, v);
  }
}

var advect = function(f, timestep) {
  var scale = timestep / f.gridSize;
  for (let i = 0;i < f.width;i++) {
    for (let j = 0;j < f.height;j++) {
      var u = f.u.sample(i, j) * scale;
      var v = f.v.sample(i, j) * scale;

      var advectedQuantity = f.q0.sample(i - u, j - v);

      f.q1.set(i, j, advectedQuantity);
    }
  }
}

// compute divergence
var computeRHS = function(rhs, f) {
  var scale = 1 / f.gridSize;
  for (let i = 0;i < f.width;i++) {
    for (let j = 0;j < f.height;j++) {
      var rhsVal = -scale * (
          f.u.get(i + 0.5, j) - f.u.get(i - 0.5, j) +
          f.v.get(i, j + 0.5) - f.v.get(i, j - 0.5)
        );
      rhs.set(i, j, rhsVal);
    }
  }
}

// project from incremental-fluid by tunabrain
// solve p from given rhs and matrix A
var pBuff0 = new Grid2D(WIDTH, HEIGHT);
var pBuff1 = new Grid2D(WIDTH, HEIGHT);
var projectJacobi = function(f, rhs, iterCount, timestep) {
  var scale = timestep / (f.density * f.gridSize * f.gridSize);
  var p0 = f.p0;
  var p1 = pBuff1;

  for (let m = 0;m < iterCount;m++) {
    if (m === iterCount - 1) {
      p1 = f.p1;
    }

    for (let x = 0;x < f.width;x++) {
      for (let y = 0;y < f.height;y++) {

        var diag = 0.0;
        var offDiag = 0.0;
        if (x > 0) {
          diag += scale;
          offDiag -= scale * p0.get(x - 1, y);
        }

        if (y > 0) {
          diag += scale;
          offDiag -= scale * p0.get(x, y - 1);
        }

        if (x < f.width - 1) {
          diag += scale;
          offDiag -= scale * p0.get(x + 1, y);
        }

        if (y < f.height - 1) {
          diag += scale;
          offDiag -= scale * p0.get(x, y + 1);
        }

        var newP = (rhs.get(x, y) - offDiag) / diag;
        p1.set(x, y, newP);
      }
    }


    if (m === 0) {
      p0 = pBuff1;
      p1 = pBuff0;
    } else {
      var ptmp = p0;
      p0 = p1;
      p1 = ptmp;
    }

  }
}

var projectGaussSeidel = function(f, rhs, iterCount, timestep) {
  var scale = timestep / (f.density * f.gridSize * f.gridSize);
  var p0 = f.p0;
  var p1 = f.p0;

  for (let m = 0;m < iterCount;m++) {
    if (m == iterCount - 1) {
      p1 = f.p1;
    }

    for (let x = 0;x < f.width;x++) {
      for (let y = 0;y < f.height;y++) {

        var diag = 0.0;
        var offDiag = 0.0;

        if (x > 0) {
          diag += scale;
          offDiag -= scale * p0.get(x - 1, y);
        }

        if (y > 0) {
          diag += scale;
          offDiag -= scale * p0.get(x, y - 1);
        }

        if (x < f.width - 1) {
          diag += scale;
          offDiag -= scale * p0.get(x + 1, y);
        }

        if (y < f.height - 1) {
          diag += scale;
          offDiag -= scale * p0.get(x, y + 1);
        }

        var newP = (rhs.get(x, y) - offDiag) / diag;
        p1.set(x, y, newP);
      }
    }
  }
  //console.log(p1.toString());
}

var applyPressure = function(f, timestep) {
  var scale = timestep / (f.density * f.gridSize);

  // set u
  for (let x = 0;x < f.width;x++) {
    for (let y = 0;y < f.height;y++) {
      var sp = scale * f.p1.get(x, y);
      f.u.set(x - 0.5, y, f.u.get(x - 0.5, y) - sp);
      f.u.set(x + 0.5, y, f.u.get(x + 0.5, y) + sp);

      f.v.set(x, y - 0.5, f.v.get(x, y - 0.5) - sp);
      f.v.set(x, y + 0.5, f.v.get(x, y + 0.5) + sp);
    }
  }

  for (let y = 0;y < f.height;y++) {
    f.u.set(-0.5, y, 0);
    f.u.set(f.width - 0.5, y, 0);
  }

  for (let x = 0;x < f.width;x++) {
    f.v.set(x, -0.5, 0);
    f.v.set(x, f.height - 0.5, 0);
  }

  //console.log(f.v.toString());
}

var fluid = new Fluid(WIDTH, HEIGHT, 0.1);

var fluidCanvas = new FluidCanvas('main', fluid);
var rhs = new Grid2D(WIDTH, HEIGHT, 0);

var simulate = function(timestep) {
  computeRHS(rhs, fluid);
  projectGaussSeidel(fluid, rhs, 300, timestep);
  applyPressure(fluid, timestep);
  advect(fluid, timestep);
  //debugger;
  //console.log(fluid)
}

var draw = function(timestep) {
  fluidCanvas.drawQuantity();
  //fluidCanvas.drawVelocity(timestep);
  //fluidCanvas.drawGrid();
}

var loop = function() {
  fluid.addInFlow(0.45, 0.2, 0.2, 0.2, 1.0, 0.0, 3.0);
  console.log('u');
  console.log(fluid.v.toString());
  console.log('q0');
  console.log(fluid.q0.toString());
  simulate(0.005);
  draw(0.005);
  console.log(fluid.q1.toString());
  fluid.swap();
  debugger;

  window.requestAnimationFrame(loop);
}

loop();
