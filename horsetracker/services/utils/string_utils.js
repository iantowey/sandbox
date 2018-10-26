'use strict';

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
    }

};