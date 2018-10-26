'use strict';

var fs = require('fs');
var request = require('request');
var cheerio = require('cheerio');
var dateFormat = require('dateformat');

var today_dt = dateFormat(new Date(), "dd-mm-yyyy");
var url = 'www.sportinglife.com/ajax/racing/profile-search?search=on+the+fringe&type=H';
console.log(url);

// var num_result_str = $('.content-header.alt p').text();
// var total_result_count = num_result_str.match(/[0-9]* (?=results)/g);

var file_name = '../../tmp/on_the_fringle_search_result' + today_dt + ".html"
//TODO: replace with post to queue for processing
fs.writeFile(file_name, html, null, function (err) {
    console.log("results saved " + file_name);
});


console.log('done');

