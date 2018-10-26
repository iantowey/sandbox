'use strict';

var request = require('request');

module.exports = {

    get: function(url, headers, callback){
        request({method:'GET', uri:url, headers: headers}, function (error, response, html) {
            if (!error) {
                callback(html);
            }
            else {
                console.error(error);
            }
        });
    }
    ,
    post: function(url, headers, body, callback) {
        request({method: 'POST', uri:url, headers: headers, body: body}, function (error, response, html) {
            if (!error) {
                callback(html);
            }
            else {
                console.error(error);
            }
        });
    }
}

