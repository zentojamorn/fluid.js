var FluidCanvas = function(id, f) {
  this.canvas = document.getElementById(id);
  this.context = this.canvas.getContext("2d");
  this.scaleX = this.canvas.width / f.width;
  this.scaleY = this.canvas.height / f.height;
}

FluidCanvas.prototype = {
  drawLine: function(beginX, beginY, endX, endY) {
    this.context.beginPath();
    this.context.moveTo(beginX * this.scaleX, beginY * this.scaleY);
    this.context.lineTo(endX * this.scaleX, endY * this.scaleY);
    this.context.stroke();
  },

  drawCircle: function(posX, posY, radius) {
    this.context.beginPath();
    this.context.arc(posX * this.scaleX, posY * this.scaleY, radius, 0, 2 * Math.PI, false);
    this.context.stroke();
  },

  drawRectangle: function(startX, startY, width, height) {
    this.context.rect(startX * this.scaleX, startY * this.scaleY, width * this.scaleX, height * this.scaleY);
    this.context.stroke();
  },

  drawRectangleFill: function(startX, startY, width, height, color) {
    this.context.fillStyle = color;
    this.context.fillRect(startX * this.scaleX, startY * this.scaleY, width * this.scaleX, height * this.scaleY);
  },

  drawGrid: function(fluid) {
    for (var i = 0;i <= fluid.width;i++) {
      this.drawLine(i, 0, i, fluid.height);
    }
    for (var j = 0;j <= fluid.height;j++) {
      this.drawLine(0, j, fluid.width, j);
    }
  },

  drawVelocity: function(fluidU, fluidV, timestep) {
    for (var i = 0;i < fluidV.width;i++) {
      for (var j = 0;j < fluidU.height;j++) {
        var u = fluidU.sample(i, j);
        var v = fluidV.sample(i, j);
        u *= timestep;
        v *= timestep;
        this.drawCircle(i + 0.5, j + 0.5, 1);
        this.drawLine(i + 0.5, j + 0.5, u + i + 0.5, v + j + 0.5);
      }
    }
  },

  drawQuantity: function(fluidQ) {
    for (var i = 0;i < fluidQ.width;i++) {
      for (var j = 0;j < fluidQ.height;j++) {
        this.drawRectangleFill(i, j, 1, 1,
          new String('rgb(255,' + Math.floor((1.0 - fluidQ.get(i, j)) * 255) + ',' + Math.floor((1.0 - fluidQ.get(i, j)) * 255) + ')')
        );
      }
    }
  }

}
