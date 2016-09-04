var WIDTH = 10;
var HEIGHT = 10;

var Fluid = function(width, height, density) {
  this.width = width;
  this.height = height;
  this.density = density;
  this.prop = new Grid2D(width,   height,     0); // property in each cell

  this.p0 = new Grid2D(width,      height,     0);
  this.p1 = new Grid2D(width,      height,     0);

  this.q0 = new Grid2D(width,      height,     0);
  this.q1 = new Grid2D(width,      height,     0);

  this.u0 = new Grid2D(width + 1,  height,     0,  0.5, 0.0);
  this.u1 = new Grid2D(width + 1,  height,     0,  0.5, 0.0);

  this.v0 = new Grid2D(width,      height + 1, 0,  0.0, 0.5);
  this.v1 = new Grid2D(width,      height + 1, 0,  0.0, 0.5);

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

    var tmpv = this.v0;
    this.v0 = this.v1;
    this.v1 = tmpv;

    var tmpu = this.u0;
    this.u0 = this.u1;
    this.u1 = tmpu;
  },

  addInFlow: function(sx, sy, w, h, d, u, v) {
    var ssx = sx / this.gridSize;
    var ssy = sy / this.gridSize;

    var sex = (sx + w) / this.gridSize;
    var sey = (sy + h) / this.gridSize;

    this.q0.add(ssx, ssy, sex, sey, d);
    this.u0.add(ssx - 0.5, ssy, sex - 0.5, sey, u);
    this.v0.add(ssx, ssy - 0.5, sex, sey - 0.5, v);
  }

}

var advect = function(f, timestep) {

  // advect quantity
  var scale = timestep / f.gridSize;
  for (let i = 0;i < f.width;i++) {
    for (let j = 0;j < f.height;j++) {
      var u = f.u0.sample(i, j) * scale;
      var v = f.v0.sample(i, j) * scale;
      var advectedQuantity = f.q0.sample(i - u, j - v);
      f.q1.set(i, j, advectedQuantity);
    }
  }

  // advect velocity u
  for (let i = 0;i < f.width + 1;i++) {
    for (let j = 0;j < f.height;j++) {
      var u = f.u0.sample(i - 0.5, j) * scale;
      var v = f.v0.sample(i - 0.5, j) * scale;
      var advectedU = f.u0.sample(i - u, j - v);
      f.u1.set(i - 0.5, j, advectedU);
    }
  }

  // advect velocity v
  for (let i = 0;i < f.width;i++) {
    for (let j = 0;j < f.height + 1;j++) {
      var u = f.u0.sample(i, j - 0.5) * scale;
      var v = f.v0.sample(i, j - 0.5) * scale;
      var advectedV = f.v0.sample(i - u, j - v);
      f.v1.set(i, j - 0.5, advectedV);
    }
  }
}

// compute divergence
var computeRHS = function(rhs, f) {
  var scale = 1 / f.gridSize;
  for (let i = 0;i < f.width;i++) {
    for (let j = 0;j < f.height;j++) {
      var rhsVal = -scale * (
          f.u0.get(i + 0.5, j) - f.u0.get(i - 0.5, j) +
          f.v0.get(i, j + 0.5) - f.v0.get(i, j - 0.5)
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
}

var applyPressure = function(f, timestep) {
  var scale = timestep / (f.density * f.gridSize);

  // set u
  for (let x = 0;x < f.width;x++) {
    for (let y = 0;y < f.height;y++) {
      var sp = scale * f.p1.get(x, y);
      f.u0.set(x - 0.5, y, f.u0.get(x - 0.5, y) - sp);
      f.u0.set(x + 0.5, y, f.u0.get(x + 0.5, y) + sp);

      f.v0.set(x, y - 0.5, f.v0.get(x, y - 0.5) - sp);
      f.v0.set(x, y + 0.5, f.v0.get(x, y + 0.5) + sp);
    }
  }

  for (let y = 0;y < f.height;y++) {
    f.u0.set(-0.5, y, 0);
    f.u0.set(f.width - 0.5, y, 0);
  }

  for (let x = 0;x < f.width;x++) {
    f.v0.set(x, -0.5, 0);
    f.v0.set(x, f.height - 0.5, 0);
  }
}

var fluid = new Fluid(WIDTH, HEIGHT, 0.1);

var fluidCanvas = new FluidCanvas('main', fluid);
var rhs = new Grid2D(WIDTH, HEIGHT, 0);

var simulate = function(timestep) {
  computeRHS(rhs, fluid);
  projectGaussSeidel(fluid, rhs, 500, timestep);
  applyPressure(fluid, timestep);
  advect(fluid, timestep);
}

var draw = function(timestep) {
  fluidCanvas.drawQuantity();
  fluidCanvas.drawVelocity(timestep * 100);
  fluidCanvas.drawGrid();
}

var loop = function() {
  fluid.addInFlow(0.45, 0.2, 0.2, 0.2, 1.0, 0.0, 3.0);
  simulate(0.002);
  draw(0.002);
  fluid.swap();
  debugger;

  window.requestAnimationFrame(loop);
}

loop();
