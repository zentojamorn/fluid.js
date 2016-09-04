// khan academy code
var calcBilinearInterpolate = function(x, y, x1, y1, x2, y2, Q11, Q21, Q12, Q22){
  var ans1 = (((x2-x)*(y2-y))/((x2-x1)*(y2-y1)))*Q11;
  var ans2 = (((x-x1)*(y2-y))/((x2-x1)*(y2-y1)))*Q21;
  var ans3 = (((x2-x)*(y-y1))/((x2-x1)*(y2-y1)))*Q12;
  var ans4 = (((x-x1)*(y-y1))/((x2-x1)*(y2-y1)))*Q22;
  return (ans1+ans2+ans3+ans4);
};

var Grid2D = function(width, height, initialValue, offsetX, offsetY) {
  this.data = new Float32Array(width * height);
  this.width = width;
  this.height = height;

  this.offsetX = (offsetX === undefined) ? 0 : offsetX;
  this.offsetY = (offsetY === undefined) ? 0 : offsetY;

  for (let i = 0;i < width * height;i++) {
    this.data[i] = initialValue;
  }
}

Grid2D.prototype = {
  get: function(x, y) { // getter
    x += this.offsetX;
    y += this.offsetY;
    return this.data[this.height * x + y];
  },

  set: function(x, y, v) { // setter
    x += this.offsetX;
    y += this.offsetY;
    x = Math.round(x);
    y = Math.round(y);
    this.data[this.height * x + y] = v;
  },

  sample: function(x, y) {
    x += this.offsetX;
    y += this.offsetY;
    var xLow = Math.floor(x).clamp(0, this.width - 1);
    var yLow = Math.floor(y).clamp(0, this.height - 1);

    var xHigh = Math.ceil(x).clamp(0, this.width - 1);
    var yHigh = Math.ceil(y).clamp(0, this.height - 1);

    var q00 = this.data[this.height * xLow + yLow];
    var q10 = this.data[this.height * xHigh + yLow];
    var q01 = this.data[this.height * xLow + yHigh];
    var q11 = this.data[this.height * xHigh + yHigh];

    var value = calcBilinearInterpolate(x, y, xLow, yLow, xLow + 1, yLow + 1, q00, q10, q01, q11);
    return value;
  },

  add:function(sx, sy, ex, ey, v) {
    sx = Math.floor(sx + this.offsetX);
    sy = Math.floor(sy + this.offsetY);
    ex = Math.floor(ex + this.offsetX);
    ey = Math.floor(ey + this.offsetY);

    sx = sx.clamp(0, this.width - 1);
    sy = sy.clamp(0, this.height - 1);
    ex = ex.clamp(0, this.width - 1);
    ey = ey.clamp(0, this.height - 1);

    for (let i = sx;i < ex;i++) {
      for (let j = sy;j < ey;j++) {
        this.data[this.height * i + j] = Math.max(v, this.data[this.height * i + j]);
      }
    }
  },

  toString:function() {
    var out = "";
    for (let y = 0;y < this.height;y++) {
      for (let x = 0;x < this.width;x++) {
        out += parseFloat(this.data[x * this.height + y]).toFixed(2) + " ";
      }
      out += "\n";
    }
    return out;
  }
}
