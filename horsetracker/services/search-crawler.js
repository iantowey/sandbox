'use strict';

var fs = require('fs');
var request = require('request');
var cheerio = require('cheerio');
var dateFormat = require('dateformat');

var today_dt = dateFormat(new Date(), "dd-mm-yyyy");
url = 'http://www.sportinglife.com/racing/racecards/' + today_dt;
console.log(url);

var num_result_str = $('.content-header.alt p').text();
var total_result_count = num_result_str.match(/[0-9]* (?=results)/g);

var file_name = '/tmp/racecard_' + today_dt + ".html"
//TODO: replace with post to queue for processing
fs.writeFile(file_name, html, null, function(err){
    console.log("racecard saved " + file_name);
        });


console.log('done');

