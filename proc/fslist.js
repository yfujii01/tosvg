var fs = require('fs');
var path = require('path');
const conf = require('./conf');

var list;
try {
    const srcDir = conf.src;
    list = fs.readdirSync(srcDir);
    for (var i = 0; i < list.length; i++) {
        var stats = fs.statSync(path.join(srcDir, list[i]));

    }
}
catch (err) {
    console.error(err);
    process.exit(1);
}

exports.fnames = list;