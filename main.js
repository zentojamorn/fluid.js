var initCanvas = function(id) {
  var canvas = document.getElementById(id);
  var context = canvas.getContext("2d");
  canvas.context = context;
  return canvas;
}

var canvas = initCanvas('main');
var WIDTH = 20;
var HEIGHT = 20;
var SCALE_WIDTH = canvas.width / WIDTH;
var SCALE_HEIGHT = canvas.height / HEIGHT;

var vec2 = function(x, y) {
  this.x = (x === undefined) ? 0 : x;
  this.y = (y === undefined) ? 0 : y;
}

vec2.prototype = {
  set: function(x, y) {
    this.x = x;
    this.y = y;
    return this;
  }
}

var Pool = function(con, initialSize) {
  this.con = con;
  this.stack = [];
  this.length = (initialSize === undefined) ? 32 : initialSize;
  for (let i = 0;i < this.length;i++) {
    this.stack.push(new this.con());
  }
}

Pool.prototype = {
  alloc: function() {
    if (this.stack.length == 0) {
      for (let i = 0;i < this.length;i++) {
        this.stack.push(new this.con());
      }
      this.length = this.length * 2;
    }
    return this.stack.pop();
  },
  free: function(obj) {
    this.stack.push(obj);
  }
}

var Vec2Pool = new Pool(vec2);

var createGrid2D = function(width, height, initvalue) {
  var grid = new Array(width);
  for (let i = 0;i < width;i++) {
    grid[i] = new Array(height);
    for (let j = 0;j < height;j++) {
      grid[i][j] = initvalue;
    }
  }
  return grid;
}

var drawLine = function(canvas, beginX, beginY, endX, endY) {
  canvas.context.beginPath();
  canvas.context.moveTo(beginX * SCALE_WIDTH, beginY * SCALE_HEIGHT);
  canvas.context.lineTo(endX * SCALE_WIDTH, endY * SCALE_HEIGHT);
  canvas.context.stroke();
}

var drawCircle = function(canvas, posX, posY, radius) {
  canvas.context.beginPath();
  canvas.context.arc(posX * SCALE_WIDTH, posY * SCALE_HEIGHT, radius, 0, 2 * Math.PI, false);
  canvas.context.stroke();
}

var drawRectangle = function(canvas, startX, startY, width, height) {
  canvas.context.rect(startX * SCALE_WIDTH, startY * SCALE_HEIGHT, width * SCALE_WIDTH, height * SCALE_HEIGHT);
  canvas.context.stroke();
}

var drawRectangleFill = function(canvas, startX, startY, width, height, color) {
  canvas.context.fillStyle = color;
  canvas.context.fillRect(startX * SCALE_WIDTH, startY * SCALE_HEIGHT, width * SCALE_WIDTH, height * SCALE_HEIGHT);
}

var drawGrid = function(canvas) {
  for (var i = 0;i <= WIDTH;i++) {
    drawLine(canvas, i, 0, i, HEIGHT);
  }
  for (var j = 0;j <= HEIGHT;j++) {
    drawLine(canvas, 0, j, WIDTH, j);
  }
}

var drawVelocity = function(canvas, u, v) {
  for (var i = 0;i < WIDTH;i++) {
    for (var j = 0;j < HEIGHT;j++) {
      var ux = (u[i][j] + u[i + 1][j]) / 2;
      var uy = (v[i][j] + v[i][j + 1]) / 2;
      drawRectangle(canvas, i + 0.5, j + 0.5, 0.001, 0.000000001);
      drawLine(canvas, i + 0.5, j + 0.5, ux + i + 0.5, uy + j + 0.5);
    }
  }
}

var drawQuantity = function(canvas, q) {
  for (var i = 0;i < WIDTH;i++) {
    for (var j = 0;j < HEIGHT;j++) {
      drawRectangleFill(canvas, i, j, 1, 1,
        new String('rgb(255,' + ((1 - q[i][j]) * 255) + ',' + ((1 - q[i][j]) * 255) + ')')
      );
    }
  }
}

var oldQ = createGrid2D(WIDTH, HEIGHT, 0);
var oldP = createGrid2D(WIDTH, HEIGHT, 0);
var oldU = createGrid2D(WIDTH + 1, HEIGHT, 0.2);
var oldV = createGrid2D(WIDTH, HEIGHT + 1, 0.2);

var q = createGrid2D(WIDTH, HEIGHT, 0);
var p = createGrid2D(WIDTH, HEIGHT, 0);
var u = createGrid2D(WIDTH + 1, HEIGHT, 0);
var v = createGrid2D(WIDTH, HEIGHT + 1, 0);

drawQuantity(canvas, oldQ);
drawGrid(canvas);
drawVelocity(canvas, oldU, oldV);

var computeVelocity = function(vel, pos) {
  vel.x = (oldU[pos.x][pos.y] + oldU[pos.x + 1][pos.y]) / 2;
  vel.y = (oldV[pos.x][pos.y] + oldV[pos.x][pos.y + 1]) / 2;
  return vel;
}

// RK2 advect
// Xmid = Xg - 1 / 2 * deltaT * u(Xg);
// Xp = Xg - deltaT * u(Xmid);
var advect = function(posP, posG, deltaTime) {
  var uVec = computeVelocity(Vec2Pool.alloc(), posG);

  var posMid = Vec2Pool.alloc();
  posMid.x = posG.x - uVec.x * deltaTime * 0.5;
  posMid.y = posG.y - uVec.y * deltaTime * 0.5;

  var uMid = computeVelocity(Vec2Pool.alloc(), posMid);

  posP.x = posG.x - uMid.x * deltaTime;
  posP.y = posG.y - uMid.y * deltaTime;

  Vec2Pool.free(uVec);
  Vec2Pool.free(uMid);
  Vec2Pool.free(posMid);
  return posP;
}

//
// + q(0, 0)   +q(1, 0)
//
// + q(0, 1)   +q(1, 1)
//
var bilinearInterpolate = function(pos, q) {
  var floorX = Math.floor(pos.x);
  var ceilX = Math.ceil(pos.x);
  var floorY = Math.floor(pos.y);
  var ceilY = Math.ceil(pos.y);

  var q00 = q[floorX][floorY];
  var q10 = q[ceilX][floorY];
  var q01 = q[floorX][ceilY];
  var q11 = q[ceilX][ceilY];

  return q00 * (pos.x - floorX) * (pos.y - floorY) +
    q10 * (ceilX - pos.x) * (pos.y - floorY) +
    q01 * (pos.x - floorX) * (ceilY - pos.y) +
    q11 * (ceilX - pos.x) * (ceilY - pos.y);
}

var advectStep = function() {
  
}
