'use strict';

var fs = require('fs');
var request = require('request');
var cheerio = require('cheerio');
var dateFormat = require('dateformat');

var today_dt = dateFormat(new Date(), "dd-mm-yyyy");
var url = 'http://www.sportinglife.com/racing/racecards/' + today_dt;
console.log(url);

request(url, function(error, response, html){
    if(!error) {
        var $ = cheerio.load(html);
        var file_name = '../tmp/racecard_' + today_dt + ".html"
        //TODO: replace with post to queue for processing
        fs.writeFile(file_name, html, null, function(err){
            console.log("racecard saved " + file_name);
        });
    }
    else{
            console.log(error);
        }
    });


console.log('done');

