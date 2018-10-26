'use strict';

var crypto = require('crypto');
var url_util = require('url');

module.exports = {
    scrub_obj_whitespace: function(obj) {
        return JSON.parse(JSON.stringify(obj).replace(/(\s+|\n)/g,' ').trim());
    },
    nvl: function(s1){
        if (s1 === null ){
            return "";
        } else {
            return s1;
        }
    },
    hash_string: function(str){
        var md5sum = crypto.createHash('md5');
        md5sum.update(str.replace(/\s+/g,''),'utf8');
        return md5sum.digest('hex');
    },
    url_parse: function(url_str){
        var url_obj = url_util.parse(url_str, true);
        url_obj.path_array = url_obj.path.split('/');
        return url_obj;
    }
};