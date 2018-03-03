const ConfigParser = require('configparser');

const config = new ConfigParser();

config.read('param.cfg');

var src = config.get('path', 'src');
var dst = config.get('path', 'dst');

// console.log(src);
// console.log(dst);

exports.src = src;
exports.dst = dst;
