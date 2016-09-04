var D = {};

D.debug = true;

D.Log = function(message) {
  if (D.debug) {
    console.log(message);
  }
}

D.Assert = function(truth, message) {
  if (D.debug && !truth) {
    throw new Error(message);
  }
}

D.Pause = function() {
  if (D.debug) {
    debugger;
  }
}
