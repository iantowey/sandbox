'use strict';

module.exports = {
    scrub_whitespace: function(str) {
        return str.replace(/(\s+|\n)/g,' ').trim();
    },
    scrub_json_arr_obj_strs:function(arr){
        for(var idx in arr){
            this.scrub_json_obj_strs(arr[idx]);
        }
    },
    scrub_json_obj_strs:function(obj){
        for (var key in obj) {
            if ((typeof obj[key]) === 'string') {
                obj[key] = this.scrub_whitespace(obj[key]);
            } else if(Array.isArray(obj[key])){
                this.scrub_json_arr_obj_strs(obj[key]);
            } else if ((typeof obj[key]) === 'object') {
                this.scrub_json_obj_strs(obj[key]);
            }
        }
    }
};