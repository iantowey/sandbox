'use strict';
var dateFormat = require('dateformat');

module.exports = {
    today_yy_mm_yyyy: function() {
        return dateFormat(new Date(), "dd-mm-yyyy");
    }
};