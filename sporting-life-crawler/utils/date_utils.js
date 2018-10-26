'use strict';
var dateFormat = require('dateformat');

module.exports = {
    formatted_date: function(dt, format) {
        return dateFormat(dt, format);
    }
};