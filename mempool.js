var MemoryPool = function(con, initialSize) {
  this.con = con;
  this.stack = [];
  this.length = (initialSize === undefined) ? 32 : initialSize;
  for (let i = 0;i < this.length;i++) {
    this.stack.push(new this.con());
  }
}

MemoryPool.prototype = {

  alloc: function() {
    if (this.stack.length == 0) {
      if (this.length >= 5000) {
        console.log("Pool: leak");
      }
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
